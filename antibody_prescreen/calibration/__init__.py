"""Calibration support: the reference populations the bands are fitted to.

Nothing here runs as part of screening. These are the tools for answering
"where should the range boundaries actually sit?", and for measuring what the
resulting bands do on antibodies they were not fitted to. Run
`python -m antibody_prescreen.calibration.fit_bands` to rebuild a band file;
see docs/calibration.md for the procedure and its pre-registration rule.
"""

from .plabdab_source import (
    PLABDAB_PAIRED_URL,
    default_cache_path,
    fetch_paired_sequences,
    load_paired_sequences,
    sample_pilot,
)

__all__ = [
    "PLABDAB_PAIRED_URL",
    "default_cache_path",
    "fetch_paired_sequences",
    "load_paired_sequences",
    "sample_pilot",
]
