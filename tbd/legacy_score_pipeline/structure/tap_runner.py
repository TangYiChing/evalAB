"""Adapter: an ABodyBuilder2 model -> TAP metrics -> Tier 1-style `Flag`s.

The vendored TAP code speaks in `MetricResult(metric_name, calculated_value,
flag)`. `fusion.py` speaks in `Flag(check, severity, region, message, weight)`.
This module is the only place the two meet, so re-tuning how TAP feeds the
triage level is a change here and nowhere else.
"""

from dataclasses import dataclass, field
from pathlib import Path

from ..checks import Flag

# --- Provisional weights -------------------------------------------------
#
# NOT CALIBRATED. These are placeholders chosen for sane relative behaviour,
# to be replaced by values fitted against a real population (see
# docs/tier2-tap-plabdab-plan.md, Step 7). Do not read meaning into them yet.
#
# Rationale for the placeholder values, so the behaviour is at least
# predictable: with fusion's T_LOW = 31.0, a single RED (6.0) sits below the
# existing odd-cysteine flag (8.0 in checks.py) and cannot on its own push an
# otherwise-clean candidate out of GO — which is deliberate, because several
# approved therapeutics are amber or red on at least one TAP axis. Three or
# more REDs (18.0) does move a candidate materially.
AMBER_WEIGHT = 2.0
RED_WEIGHT = 6.0

# TAP metric name -> (short label, which residues it is computed on).
#
# The residue set matters for triage, because it decides whether a finding has
# a CDR/framework locus at all. Verified against the vendored calculators:
#
#   Total IMGT CDR Length   res.is_cdr            -> CDR
#   Hydrophobic Patch (PSH) res.in_cdr_vicinity   -> CDR region
#   Positive Patch (PPC)    res.in_cdr_vicinity   -> CDR region
#   Negative Patch (PNC)    res.in_cdr_vicinity   -> CDR region
#   SFvCSP                  res.is_surface (H+L)  -> WHOLE MOLECULE, no locus
#
# So four of the five are CDR-region properties and one is a global surface
# property. SFvCSP is therefore classified on the immunogenicity axis rather
# than the repair-locus axis: an Fv-wide charge asymmetry is not something you
# fix with a substitution.
METRIC_LABELS = {
    "Total IMGT CDR Length": "total CDR length",
    "Hydrophobic Patch Score": "PSH (CDR-vicinity hydrophobic patch)",
    "Positive Patch Score": "PPC (CDR-vicinity positive patch)",
    "Negative Patch Score": "PNC (CDR-vicinity negative patch)",
    "SFvCSP": "SFvCSP (Fv-wide charge symmetry)",
}

NO_LOCUS_METRICS = {"SFvCSP"}

# Green/amber boundaries, copied from the vendored calculators so a flag can
# report WHAT was crossed and BY HOW MUCH. "crossed a boundary" on its own is
# not an actionable thing to tell someone.
METRIC_BOUNDS = {
    "Total IMGT CDR Length": ((43, 55), (37, 63)),
    "Hydrophobic Patch Score": ((137.61, 200.71), (106.44, 225.85)),
    "Positive Patch Score": ((0, 1.19), (0, 3.58)),
    "Negative Patch Score": ((0, 1.67), (0, 3.50)),
    "SFvCSP": ((-4.20, 1e5), (-20.50, 1e5)),
}


def describe_excursion(metric: str, value: float, colour: str) -> str:
    """Say which side of which boundary the value fell, and by how much."""
    bounds = METRIC_BOUNDS.get(metric)
    if bounds is None:
        return f"{value:.2f}"
    (green_lo, green_hi), (amber_lo, amber_hi) = bounds
    lo, hi = (amber_lo, amber_hi) if colour == "RED" else (green_lo, green_hi)
    limit_name = "all known therapeutics" if colour == "RED" else "the green range"
    if value > hi:
        return f"{value:.2f}, above the {hi:g} upper bound of {limit_name}"
    if value < lo:
        return f"{value:.2f}, below the {lo:g} lower bound of {limit_name}"
    return f"{value:.2f}"


class TapError(Exception):
    """Raised when TAP cannot be computed for a model (bad PDB, psa failure)."""


@dataclass
class TapProfile:
    """The five TAP metrics for one antibody."""

    values: dict[str, float] = field(default_factory=dict)
    flags: dict[str, str] = field(default_factory=dict)  # metric -> GREEN/AMBER/RED

    @property
    def n_red(self) -> int:
        return sum(1 for f in self.flags.values() if f == "RED")

    @property
    def n_amber(self) -> int:
        return sum(1 for f in self.flags.values() if f == "AMBER")

    @property
    def worst_flag(self) -> str:
        if self.n_red:
            return "RED"
        if self.n_amber:
            return "AMBER"
        return "GREEN"

    def summary(self) -> str:
        """Compact one-cell summary for the report table, e.g. '2R/1A'."""
        if self.n_red or self.n_amber:
            parts = []
            if self.n_red:
                parts.append(f"{self.n_red}R")
            if self.n_amber:
                parts.append(f"{self.n_amber}A")
            return "/".join(parts)
        return "all green"


def run_tap_profile(model_path: Path | str) -> TapProfile:
    """Run the five TAP metrics on an IMGT-numbered ABodyBuilder2 model.

    Args:
        model_path: path to a .pdb with chains named H and L, IMGT numbered.

    Raises:
        TapError: if the model cannot be read or psa fails.
    """
    from .vendor.tap.main import run_tap

    try:
        results = run_tap(str(model_path), outfile=None, quiet=True)
    except Exception as e:
        raise TapError(f"TAP failed on {model_path}: {e}") from e

    return TapProfile(
        values={r.metric_name: float(r.calculated_value) for r in results},
        flags={r.metric_name: r.flag for r in results},
    )


def tap_flags(profile: TapProfile) -> list[Flag]:
    """Convert a TapProfile into weighted soft Flags for fusion.

    Every TAP flag is a *soft* flag, never a hard gate. A RED means "outside
    the range spanned by known therapeutics on this one axis" — real approved
    antibodies do sit outside it, so treating it as an unconditional reject
    would repeat the mistake the odd-cysteine hard gate already made (see the
    package README). Keep it recoverable and weighted.
    """
    flags: list[Flag] = []
    for metric, colour in profile.flags.items():
        if colour == "GREEN":
            continue
        label = METRIC_LABELS.get(metric, metric)
        value = profile.values.get(metric, float("nan"))
        # One check name per metric, so each can be measured and routed on its
        # own. "any TAP RED" was a single gate whose likelihood ratio of 11 was
        # carried entirely by PSH and PNC — measured on 300 antibodies, CDR
        # length, PPC and SFvCSP produced zero RED events in either group.
        flags.append(
            Flag(
                check=f"tap_{_slug(metric)}",
                severity="soft",
                region="structure" if metric in NO_LOCUS_METRICS else "CDR-region",
                message=(
                    f"{label} is {colour.lower()}: "
                    f"{describe_excursion(metric, value, colour)}"
                ),
                weight=RED_WEIGHT if colour == "RED" else AMBER_WEIGHT,
            )
        )
    return flags


def _slug(metric: str) -> str:
    return {
        "Total IMGT CDR Length": "cdr_length",
        "Hydrophobic Patch Score": "psh",
        "Positive Patch Score": "ppc",
        "Negative Patch Score": "pnc",
        "SFvCSP": "sfvcsp",
    }.get(metric, metric.lower().replace(" ", "_"))
