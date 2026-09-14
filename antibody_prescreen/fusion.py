"""Score fusion: run every Tier 1 check on a VH/VL candidate, collapse the
result into one verdict (GO / CONDITIONAL / NO-GO) with the top 1-2 reasons.

This is the piece that solves the "click fatigue" problem from the design
discussion: all checks still run and all detail is retained, but the default
view is one line per candidate.

Tier 2 (structure-based TAP profiling) is opt-in via run_structure=True. Its
dependencies (torch, ImmuneBuilder, OpenMM) are imported lazily inside that
branch, so a Tier-1-only user never needs them installed.
"""

from dataclasses import dataclass, field

from .checks import Flag, run_all_checks
from .immunogenicity import check_immunogenicity
from .numbering import NumberingError, number_antibody

# Calibrated against 508 real, structurally-solved antibodies from the
# ANTIPASTI-curated SAbDab dataset (see data/sabdab_derived_candidates.json
# and sabdab_source.py) — RCSB/SAbDab/TheraSAbDab were all unreachable from
# the build environment, so this dataset (GitHub-hosted, derived from real
# PDB structures) stood in for the originally planned TheraSAbDab set.
#
# This is a real improvement over the original 3-candidate placeholder, but
# it is NOT the same thing as the plan's original Stage D: this population
# is "antibodies that were solvable/crystallizable" (a structural-biology
# selection), not "antibodies with known real-world clinical/manufacturing
# outcomes." A structure existing says nothing about developability,
# immunogenicity in patients, or manufacturing success — cross-referencing
# against literature-documented liabilities (the original Stage D ask) is
# still outstanding.
#
# Score distribution across 466 non-hard-gated candidates: min=15.59,
# median=31.17, p90=39.72, max=50.91. Thresholds set at the median (roughly
# half of a real antibody population lands GO) and the 90th percentile
# (worst ~10% lands NO-GO), rather than the earlier placeholder values which
# put effectively 100% of real structures below GO.
T_LOW = 31.0
T_HIGH = 40.0


@dataclass
class CandidateResult:
    candidate_id: str
    verdict: str  # "GO", "CONDITIONAL", "NO-GO", or "ERROR"
    score: float
    top_reasons: list[str]
    all_flags: list[Flag] = field(default_factory=list)
    error: str | None = None
    immunogenicity_available: bool = False
    structure_available: bool = False
    tap_profile: "TapProfile | None" = None  # noqa: F821 - lazy Tier 2 import


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
) -> CandidateResult:
    """Screen one VH/VL candidate.

    Args:
        run_immunogenicity: run the IEDB MHC-II scan (network; see README).
        run_structure: run Tier 2 — build an ABodyBuilder2 model and compute
            the TAP profile. Needs the optional structure dependencies. If
            they are missing or modelling fails, the candidate still gets a
            Tier 1 verdict with structure_available=False, rather than
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
            verdict="ERROR",
            score=0.0,
            top_reasons=[str(e)],
            error=str(e),
        )

    all_flags: list[Flag] = []
    immuno_available = False

    for chain_name, chain in chains.items():
        result = run_all_checks(chain, chain_name)
        all_flags.extend(result.flags)

        if run_immunogenicity:
            immuno_result = check_immunogenicity(chain, chain_name)
            all_flags.extend(immuno_result.flags)
            immuno_available = immuno_available or immuno_result.available

    structure_available = False
    tap_profile = None
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

    hard_gates = [f for f in all_flags if f.severity == "hard_gate"]
    soft_flags = [f for f in all_flags if f.severity == "soft"]

    if hard_gates:
        return CandidateResult(
            candidate_id=candidate_id,
            verdict="NO-GO",
            score=float("inf"),
            top_reasons=[hard_gates[0].message],
            all_flags=all_flags,
            immunogenicity_available=immuno_available,
            structure_available=structure_available,
            tap_profile=tap_profile,
        )

    score = sum(f.weight for f in soft_flags)

    if score < T_LOW:
        verdict = "GO"
    elif score < T_HIGH:
        verdict = "CONDITIONAL"
    else:
        verdict = "NO-GO"

    return CandidateResult(
        candidate_id=candidate_id,
        verdict=verdict,
        score=round(score, 2),
        top_reasons=_top_reasons(soft_flags) if verdict != "GO" else [],
        all_flags=all_flags,
        immunogenicity_available=immuno_available,
        structure_available=structure_available,
        tap_profile=tap_profile,
    )


def screen_batch(
    candidates: list[dict],
    run_immunogenicity: bool = True,
    run_structure: bool = False,
    model_cache=None,
) -> list[CandidateResult]:
    """candidates: list of {"candidate_id", "vh_sequence", "vl_sequence"} dicts.

    With run_structure=True the ABodyBuilder2 predictor is constructed once
    for the whole batch and reused — building it loads ~100MB of weights, so
    doing it per candidate would dominate the runtime of a large batch.
    """
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
        )
        for c in candidates
    ]


def format_report(results: list[CandidateResult]) -> str:
    """One-line-per-candidate markdown table, ranked GO first by score.

    The TAP column appears only when at least one result actually has a
    structure profile — a Tier-1-only run keeps exactly the table it had
    before Tier 2 existed.
    """
    order = {"GO": 0, "CONDITIONAL": 1, "NO-GO": 2, "ERROR": 3}
    ranked = sorted(results, key=lambda r: (order[r.verdict], r.score))

    show_tap = any(r.structure_available for r in results)

    if show_tap:
        lines = [
            "| Candidate | Verdict | Score | TAP | Why |",
            "|---|---|---|---|---|",
        ]
    else:
        lines = ["| Candidate | Verdict | Score | Why |", "|---|---|---|---|"]

    for r in ranked:
        score_str = "—" if r.score == float("inf") else str(r.score)
        why = "; ".join(r.top_reasons) if r.top_reasons else "—"
        if show_tap:
            tap = r.tap_profile.summary() if r.tap_profile is not None else "n/a"
            lines.append(
                f"| {r.candidate_id} | **{r.verdict}** | {score_str} | {tap} | {why} |"
            )
        else:
            lines.append(
                f"| {r.candidate_id} | **{r.verdict}** | {score_str} | {why} |"
            )
    return "\n".join(lines)
