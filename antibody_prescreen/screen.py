"""Public screening workflow: profile -> reference bands -> repair triage."""

from dataclasses import asdict, dataclass
from pathlib import Path

from .bands import CLINICAL_METRICS, BandSet, compose
from .levels import Triage, assign
from .numbering import NumberingError, number_antibody
from .profile import Profile, build_profile

DATA_DIR = Path(__file__).resolve().parent / "data"
DEFAULT_BANDS = DATA_DIR / "bands_human_repertoire.json"
DEFAULT_CLINICAL_BANDS = DATA_DIR / "bands_clinical.json"


@dataclass
class StructureAnnotation:
    """Optional Tier 2 result. It is context only and never changes `triage`."""

    available: bool
    note: str = ""
    model_path: str | None = None
    confidence: dict | None = None
    tap_values: dict[str, float] | None = None
    tap_flags: dict[str, str] | None = None


@dataclass
class RunProvenance:
    """What a reader needs in order to know what a report is worth.

    A level is meaningless without the reference population that defined the
    range, and a batch summary is misleading without the share of candidates
    that could not be numbered at all. Both travel with every report so that
    neither can be quietly dropped when results are copied into a slide.
    """

    band_reference: str
    band_n_derivation: int
    band_n_validation: int
    band_cuts: tuple
    n_candidates: int
    n_numbered: int
    unmeasured_metrics: list[str]

    @property
    def numbering_success_rate(self) -> float:
        return self.n_numbered / self.n_candidates if self.n_candidates else 0.0

    def to_dict(self) -> dict:
        return {
            "band_version": {
                "reference": self.band_reference,
                "n_derivation": self.band_n_derivation,
                "n_validation": self.band_n_validation,
                "cuts": list(self.band_cuts),
            },
            "n_candidates": self.n_candidates,
            "n_numbered": self.n_numbered,
            "numbering_success_rate": round(self.numbering_success_rate, 4),
            "unmeasured_metrics": self.unmeasured_metrics,
        }

    def as_lines(self) -> list[str]:
        pct = 100 * self.numbering_success_rate
        unmeasured = ", ".join(self.unmeasured_metrics) or "none"
        return [
            f"- **Reference bands:** {self.band_reference} "
            f"(fit n={self.band_n_derivation}, validation n={self.band_n_validation}, "
            f"cuts p{'/p'.join(f'{c:g}' for c in self.band_cuts)})",
            f"- **Numbering:** {self.n_numbered}/{self.n_candidates} candidates numbered ({pct:.1f}%)",
            f"- **Metrics with no reference band (unmeasured, NOT typical):** {unmeasured}",
        ]


@dataclass
class ScreenResult:
    candidate_id: str
    triage: Triage
    profile: Profile
    structure: StructureAnnotation | None = None

    def to_dict(self) -> dict:
        return {
            "candidate_id": self.candidate_id,
            "level": self.triage.level,
            "level_name": self.triage.name,
            "action": self.triage.action,
            "drivers": [asdict(finding) for finding in self.triage.drivers],
            "findings": [asdict(finding) for finding in self.triage.findings],
            "unmeasured": self.triage.unmeasured,
            "error": self.triage.error,
            "structure": asdict(self.structure) if self.structure else None,
        }


def provenance(results: list[ScreenResult], band_set: BandSet) -> RunProvenance:
    """Summarise a completed run: which bands, how many numbered, what was not measured."""
    unmeasured = sorted({
        entry.split("=")[0] for result in results for entry in result.triage.unmeasured
    })
    return RunProvenance(
        band_reference=band_set.reference,
        band_n_derivation=band_set.n_derivation,
        band_n_validation=band_set.n_validation,
        band_cuts=band_set.cuts,
        n_candidates=len(results),
        n_numbered=sum(1 for r in results if not r.profile.error),
        unmeasured_metrics=unmeasured,
    )


def load_default_bands() -> BandSet:
    """Load repertoire bands plus clinical humanness bands supplied with evalAB."""
    return compose(
        BandSet.load(DEFAULT_BANDS), BandSet.load(DEFAULT_CLINICAL_BANDS), CLINICAL_METRICS
    )


def screen_candidate(
    candidate_id: str,
    vh_sequence: str,
    vl_sequence: str,
    *,
    band_set: BandSet | None = None,
    run_structure: bool = False,
    model_cache: str | Path | None = None,
    predictor=None,
) -> ScreenResult:
    """Screen one complete paired human VH/VL candidate.

    The declared scope is a single variable heavy and variable light domain.
    A non-human, truncated, non-standard, or wrongly placed sequence is kept
    in the report as a scope/error finding rather than silently discarded.
    """
    band_set = band_set or load_default_bands()
    try:
        profile = build_profile(candidate_id, number_antibody(vh_sequence, vl_sequence))
    except NumberingError as exc:
        profile = Profile(candidate_id=candidate_id, error=str(exc))
    triage = assign(profile, band_set)
    structure = _annotate_structure(vh_sequence, vl_sequence, model_cache, predictor) if run_structure and not profile.error else None
    return ScreenResult(candidate_id, triage, profile, structure)


def screen_batch(
    candidates: list[dict],
    *,
    band_set: BandSet | None = None,
    run_structure: bool = False,
    model_cache: str | Path | None = None,
) -> list[ScreenResult]:
    """Screen a batch with required candidate_id, vh_sequence, vl_sequence keys."""
    _validate_candidates(candidates)
    predictor = None
    if run_structure:
        try:
            from .structure.modelling import get_predictor
            predictor = get_predictor()
        except Exception:
            # Individual results record the Tier 2 availability failure.
            predictor = None
    return [
        screen_candidate(
            row["candidate_id"], row["vh_sequence"], row["vl_sequence"],
            band_set=band_set, run_structure=run_structure,
            model_cache=model_cache, predictor=predictor,
        )
        for row in candidates
    ]


def format_report(
    results: list[ScreenResult],
    *,
    shared_scaffold: bool = False,
    run_provenance: RunProvenance | None = None,
) -> str:
    """Render the action-oriented report; add an introduced-only view on request.

    `shared_scaffold=True` is an explicit declaration by the user, not an
    inference from sequences. It is appropriate only for a congeneric design
    batch whose members share a parent/framework — variants of one parent,
    where every candidate carries the same framework findings. Declaring it
    for an unrelated batch produces an introduced-only view with no meaning,
    which is why it is never inferred.

    `run_provenance` is rendered as a header. Pass it: a level without the
    reference population that defined it, or without the share of candidates
    that failed to number, is not an interpretable result.
    """
    lines = []
    if run_provenance is not None:
        lines += ["### Run provenance", ""] + run_provenance.as_lines() + [""]
    lines += [
        "| Candidate | Level | What to do | Driver |",
        "|---|---|---|---|",
    ]
    for result in sorted(results, key=lambda r: (r.triage.level, r.candidate_id)):
        t = result.triage
        driver = t.drivers[0].message if t.drivers else (t.error or "nothing outside the reference range")
        lines.append(f"| {result.candidate_id} | **L{t.level} {t.name}** | {t.action} | {driver} |")

    if shared_scaffold and results:
        common = set.intersection(*[
            {(f.source, f.message) for f in result.triage.findings} for result in results
        ])
        lines.extend([
            "", "### Introduced-only view", "",
            "Findings shared by every candidate are scaffold context; they do not distinguish variants.",
            "", "| Candidate | Level from introduced findings | Introduced driver |",
            "|---|---|---|",
        ])
        for result in sorted(results, key=lambda r: r.candidate_id):
            introduced = [f for f in result.triage.findings if (f.source, f.message) not in common]
            level = max((f.level for f in introduced), default=1)
            drivers = [f.message for f in introduced if f.level == level]
            lines.append(f"| {result.candidate_id} | L{level} | {'; '.join(drivers[:2]) or 'nothing beyond the shared scaffold'} |")
    return "\n".join(lines)


def _annotate_structure(vh, vl, model_cache, predictor) -> StructureAnnotation:
    try:
        from .structure.modelling import build_model, model_confidence
        from .structure.tap import run_tap_profile
        model = build_model(vh, vl, cache_dir=model_cache, predictor=predictor)
        tap = run_tap_profile(model)
        return StructureAnnotation(True, model_path=str(model), confidence=model_confidence(model), tap_values=tap.values, tap_flags=tap.flags)
    except Exception as exc:
        return StructureAnnotation(False, note=f"Tier 2 unavailable: {exc}")


def _validate_candidates(candidates: list[dict]) -> None:
    required = {"candidate_id", "vh_sequence", "vl_sequence"}
    seen = set()
    for index, row in enumerate(candidates, start=1):
        missing = required - row.keys()
        if missing:
            raise ValueError(f"row {index}: missing required column(s): {', '.join(sorted(missing))}")
        if not str(row["candidate_id"]).strip():
            raise ValueError(f"row {index}: candidate_id is empty")
        if row["candidate_id"] in seen:
            raise ValueError(f"duplicate candidate_id: {row['candidate_id']}")
        seen.add(row["candidate_id"])
