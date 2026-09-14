"""Score fusion: run every Tier 1 check on a VH/VL candidate, collapse the
result into one verdict (GO / CONDITIONAL / NO-GO) with the top 1-2 reasons.

This is the piece that solves the "click fatigue" problem from the design
discussion: all checks still run and all detail is retained, but the default
view is one line per candidate.
"""

from dataclasses import dataclass, field

from .checks import Flag, run_all_checks
from .immunogenicity import check_immunogenicity
from .numbering import NumberingError, number_antibody

# Placeholder thresholds — NOT calibrated yet. Stage D (mini calibration on a
# handful of approved mAbs) sets real values; until then these are a rough
# first guess so the pipeline is runnable end-to-end today.
T_LOW = 15.0
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


def _top_reasons(flags: list[Flag], n: int = 2) -> list[str]:
    ranked = sorted(flags, key=lambda f: f.weight, reverse=True)
    return [f.message for f in ranked[:n]]


def screen_candidate(
    candidate_id: str,
    vh_sequence: str,
    vl_sequence: str,
    run_immunogenicity: bool = True,
) -> CandidateResult:
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
    )


def screen_batch(
    candidates: list[dict], run_immunogenicity: bool = True
) -> list[CandidateResult]:
    """candidates: list of {"candidate_id", "vh_sequence", "vl_sequence"} dicts."""
    return [
        screen_candidate(
            c["candidate_id"], c["vh_sequence"], c["vl_sequence"], run_immunogenicity
        )
        for c in candidates
    ]


def format_report(results: list[CandidateResult]) -> str:
    """One-line-per-candidate markdown table, ranked GO first by score."""
    order = {"GO": 0, "CONDITIONAL": 1, "NO-GO": 2, "ERROR": 3}
    ranked = sorted(results, key=lambda r: (order[r.verdict], r.score))

    lines = ["| Candidate | Verdict | Score | Why |", "|---|---|---|---|"]
    for r in ranked:
        score_str = "—" if r.score == float("inf") else str(r.score)
        why = "; ".join(r.top_reasons) if r.top_reasons else "—"
        lines.append(f"| {r.candidate_id} | **{r.verdict}** | {score_str} | {why} |")
    return "\n".join(lines)
