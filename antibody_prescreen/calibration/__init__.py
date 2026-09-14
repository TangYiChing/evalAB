"""Calibration support: population sources for setting fusion thresholds.

Nothing here runs as part of screening. These are the tools for answering
"where should T_LOW / T_HIGH / the TAP weights actually sit?" — see
docs/tier2-tap-plabdab-plan.md, Step 7.
"""

from .plabdab_source import (
    PLABDAB_PAIRED_URL,
    fetch_paired_sequences,
    load_paired_sequences,
    sample_pilot,
)

__all__ = [
    "PLABDAB_PAIRED_URL",
    "fetch_paired_sequences",
    "load_paired_sequences",
    "sample_pilot",
]
