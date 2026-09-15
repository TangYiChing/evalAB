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

# TAP metric name -> short label used in flag messages.
METRIC_LABELS = {
    "Total IMGT CDR Length": "total CDR length",
    "Hydrophobic Patch Score": "PSH (CDR-vicinity hydrophobic patch)",
    "Positive Patch Score": "PPC (CDR-vicinity positive patch)",
    "Negative Patch Score": "PNC (CDR-vicinity negative patch)",
    "SFvCSP": "SFvCSP (Fv charge symmetry)",
}


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
        value = profile.values.get(metric)
        flags.append(
            Flag(
                check="tap",
                severity="soft",
                region="structure",
                message=f"TAP {colour.lower()} flag: {label} = {value:.2f}",
                weight=RED_WEIGHT if colour == "RED" else AMBER_WEIGHT,
            )
        )
    return flags
