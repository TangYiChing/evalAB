"""TAP annotation without a score, threshold, or routing rule.

TAP's published GREEN/AMBER/RED bands are retained as an optional structural
annotation. They are deliberately not fed into evalAB's range-and-repair
matrix: those bands were themselves fitted on therapeutic antibodies and have
not been recalibrated for this workflow.
"""

from dataclasses import dataclass, field
from pathlib import Path


class TapError(Exception):
    """Raised when TAP cannot be calculated for a predicted structure."""


@dataclass
class TapProfile:
    values: dict[str, float] = field(default_factory=dict)
    flags: dict[str, str] = field(default_factory=dict)

    def summary(self) -> str:
        red = sum(flag == "RED" for flag in self.flags.values())
        amber = sum(flag == "AMBER" for flag in self.flags.values())
        if not (red or amber):
            return "all green"
        return "/".join(part for part in (f"{red}R" if red else "", f"{amber}A" if amber else "") if part)


def run_tap_profile(model_path: Path | str) -> TapProfile:
    """Calculate the five TAP metrics for an IMGT-numbered H/L PDB model."""
    from .vendor.tap.main import run_tap

    try:
        results = run_tap(str(model_path), outfile=None, quiet=True)
    except Exception as exc:
        raise TapError(f"TAP failed on {model_path}: {exc}") from exc
    return TapProfile(
        values={result.metric_name: float(result.calculated_value) for result in results},
        flags={result.metric_name: result.flag for result in results},
    )
