"""Score fusion: run every Tier 1 check on a VH/VL candidate, collapse the
result into one triage level (1 = straight to the bench, 5 = do not pursue)
with the top 1-2 reasons.

This is the piece that solves the "click fatigue" problem from the design
discussion: all checks still run and all detail is retained, but the default
view is one line per candidate.

Tier 2 (structure-based TAP profiling) is opt-in via run_structure=True. Its
dependencies (torch, ImmuneBuilder, OpenMM) are imported lazily inside that
branch, so a Tier-1-only user never needs them installed.
"""

from dataclasses import dataclass, field

from .checks import Flag, run_all_checks
from .humanness import germline_identity, humanness_flags
from .triage import TriageResult, triage
from .immunogenicity import check_immunogenicity
from .numbering import NumberingError, number_antibody

# Calibrated against 300 human-framework antibodies drawn from PLAbDab:
# 150 clinical-stage therapeutics (pairing == "TheraSAbDab") and 150
# patent-text entries. Full run and analysis in
# docs/calibration-step7-results.md; reproduce with
# `python -m antibody_prescreen.calibration.run_calibration`.
#
# Thresholds are PERCENTILE PLACEMENTS on the therapeutic population, not a
# fitted classifier boundary: T_LOW at its median (half of a clinical-grade
# population passes) and T_HIGH at its p90 (the worst tenth is rejected).
# That makes the scale mean something specific — "where does this candidate
# sit against antibodies that reached the clinic" — which the previous
# 508-structure SAbDab basis could not, being selected for crystallisability.
#
# What this does NOT mean, and the number to keep in view: the score barely
# separates the two populations. Measured AUC (probability a random
# patent-text antibody scores worse than a random therapeutic):
#
#   tier 1 only                  0.578
#   tier 1 + TAP                 0.587   (permutation p = 0.005)
#   tier 1 + TAP + immunogenicity 0.557  <- worse; immunogenicity is now
#                                           report-only, weight 0
#
# 0.587 is statistically real and practically weak. Treat the score as a
# triage ORDERING, not a classifier. The one clean categorical signal found:
# 0/150 therapeutics carry any TAP RED flag versus 5/145 patent-text
# (Fisher exact p = 0.028) — low sensitivity, but no approved-stage antibody
# in this sample tripped one.
#
# The honest limit is the label. "Reached the clinic" versus "was patented"
# is provenance, not an assay, and most patented antibodies are real programs
# rather than developability failures. A stronger threshold needs a
# population with measured outcomes.
#
# These values happen to sit close to the previous SAbDab-derived 31.0/40.0.
# The numbers moved little; what changed is what they are anchored to.
T_LOW = 29.41
T_HIGH = 38.31


@dataclass
class CandidateResult:
    candidate_id: str
    score: float
    top_reasons: list[str]
    all_flags: list[Flag] = field(default_factory=list)
    error: str | None = None
    immunogenicity_available: bool = False
    structure_available: bool = False
    tap_profile: "TapProfile | None" = None  # noqa: F821 - lazy Tier 2 import
    triage_result: TriageResult | None = None
    model_confidence: dict | None = None


def _top_reasons(flags: list[Flag], n: int = 2) -> list[str]:
    ranked = sorted(flags, key=lambda f: f.weight, reverse=True)
    return [f.message for f in ranked[:n]]


def screen_candidate(
    candidate_id: str,
    vh_sequence: str,
    vl_sequence: str,
    run_immunogenicity: bool = True,
    run_structure: bool = False,
    model_cache=None,
    predictor=None,
    iedb_cache=None,
) -> CandidateResult:
    """Screen one VH/VL candidate.

    Args:
        run_immunogenicity: run the IEDB MHC-II scan (network; see README).
        run_structure: run Tier 2 — build an ABodyBuilder2 model and compute
            the TAP profile. Needs the optional structure dependencies. If
            they are missing or modelling fails, the candidate still gets a
            Tier 1 level with structure_available=False, rather than
            erroring out — same contract as immunogenicity.
        model_cache: directory for cached .pdb models. Tier 2 only.
        predictor: a pre-built ABodyBuilder2 instance, to avoid reloading
            weights per candidate. See screen_batch.
    """
    try:
        chains = number_antibody(vh_sequence, vl_sequence)
    except NumberingError as e:
        return CandidateResult(
            candidate_id=candidate_id,
            score=0.0,
            top_reasons=[str(e)],
            error=str(e),
        )

    all_flags: list[Flag] = []
    immuno_available = False

    for chain_name, chain in chains.items():
        result = run_all_checks(chain, chain_name)
        all_flags.extend(result.flags)

        all_flags.extend(humanness_flags(germline_identity(chain), chain_name))

        if run_immunogenicity:
            immuno_result = check_immunogenicity(
                chain, chain_name, cache_dir=iedb_cache
            )
            all_flags.extend(immuno_result.flags)
            immuno_available = immuno_available or immuno_result.available

    structure_available = False
    tap_profile = None
    model_conf = None
    if run_structure:
        # Imported here, not at module scope: Tier 2 pulls in torch /
        # ImmuneBuilder / OpenMM, and Tier-1-only users must not need them.
        from .structure import (
            ModellingError,
            TapError,
            build_model,
            run_tap_profile,
            tap_flags,
        )

        try:
            model_path = build_model(
                vh_sequence, vl_sequence, cache_dir=model_cache, predictor=predictor
            )
            tap_profile = run_tap_profile(model_path)
            all_flags.extend(tap_flags(tap_profile))
            structure_available = True

            # Structure-based disulfide pairing supersedes the sequence-level
            # cysteine count, which can only say "the total is odd" and fires
            # on 5.3% of clinical-stage therapeutics.
            from .structure.disulfide import disulfide_flags, find_disulfides
            from .structure.modelling import model_confidence as _confidence

            all_flags.extend(disulfide_flags(find_disulfides(model_path)))
            model_conf = _confidence(model_path)
        except (ModellingError, TapError) as e:
            # Degrade to Tier 1 for this candidate rather than failing the
            # batch. The flag records *why* there is no structure signal, so
            # a silently-Tier-1 candidate is visible in all_flags.
            all_flags.append(
                Flag(
                    check="tap",
                    severity="soft",
                    region="structure",
                    message=f"Tier 2 unavailable: {e}",
                    weight=0.0,
                )
            )

    # Routing is by triage level, not by the summed score. `score` is kept as
    # a diagnostic only — nothing branches on it. See triage.py for why.
    triage_result = triage(all_flags)
    soft_flags = [f for f in all_flags if f.severity == "soft"]
    score = sum(f.weight for f in soft_flags)

    return CandidateResult(
        candidate_id=candidate_id,
        score=float("inf") if triage_result.level == 5 else round(score, 2),
        top_reasons=triage_result.reasons[:2],
        all_flags=all_flags,
        immunogenicity_available=immuno_available,
        structure_available=structure_available,
        tap_profile=tap_profile,
        triage_result=triage_result,
        model_confidence=model_conf,
    )


def screen_batch(
    candidates: list[dict],
    run_immunogenicity: bool = True,
    run_structure: bool = False,
    model_cache=None,
    iedb_cache=None,
) -> list[CandidateResult]:
    """candidates: list of {"candidate_id", "vh_sequence", "vl_sequence"} dicts.

    With run_structure=True the ABodyBuilder2 predictor is constructed once
    for the whole batch and reused — building it loads ~100MB of weights, so
    doing it per candidate would dominate the runtime of a large batch.
    """
    if run_immunogenicity:
        # One batched request per ~200 chains instead of one per chain. For a
        # 1000-candidate batch that is ~10 requests and ~10 minutes instead of
        # 2000 requests and ~90 minutes, and it removes the need to run IEDB
        # calls concurrently — which is what caused the order-correlated
        # throttling during Step 7 calibration.
        from .immunogenicity import prefetch_mhcii

        prefetch_mhcii(
            [c["vh_sequence"] for c in candidates]
            + [c["vl_sequence"] for c in candidates],
            cache_dir=iedb_cache,
        )

    predictor = None
    if run_structure:
        from .structure.modelling import ModellingError, get_predictor

        try:
            predictor = get_predictor()
        except ModellingError:
            # Leave predictor=None: each candidate then fails its own Tier 2
            # step and records the reason, instead of the batch dying here.
            predictor = None

    return [
        screen_candidate(
            c["candidate_id"],
            c["vh_sequence"],
            c["vl_sequence"],
            run_immunogenicity=run_immunogenicity,
            run_structure=run_structure,
            model_cache=model_cache,
            predictor=predictor,
            iedb_cache=iedb_cache,
        )
        for c in candidates
    ]


def format_report(results: list[CandidateResult]) -> str:
    """One-line-per-candidate markdown table, ranked GO first by score.

    The TAP column appears only when at least one result actually has a
    structure profile — a Tier-1-only run keeps exactly the table it had
    before Tier 2 existed.
    """
    # Lowest level first: Level 1 is the candidate to act on, Level 5 the one
    # not to pursue. Same direction as the Emergency Severity Index.
    ranked = sorted(
        results,
        key=lambda r: (99 if r.triage_result is None else r.triage_result.level),
    )
    lines = [
        "| Candidate | Level | What to do | Repairs | Why |",
        "|---|---|---|---|---|",
    ]
    for r in ranked:
        if r.triage_result is None:
            lines.append(f"| {r.candidate_id} | ERROR | {r.error or '—'} | — | — |")
            continue
        t = r.triage_result
        repairs = (
            "—" if t.level == 5 else f"{t.cdr_repairs} CDR / {t.framework_repairs} FR"
        )
        why = "; ".join(t.reasons[:2]) if t.reasons else "—"
        lines.append(
            f"| {r.candidate_id} | **L{t.level}** | {t.action} | {repairs} | {why} |"
        )
    return "\n".join(lines)
