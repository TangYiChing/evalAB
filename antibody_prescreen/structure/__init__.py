"""Tier 2 (structure-based) screening: ABodyBuilder2 model -> TAP profile.

Everything in this subpackage is OPTIONAL. Importing `antibody_prescreen`
must not import anything from here, because Tier 2 needs heavy dependencies
(torch, ImmuneBuilder, OpenMM) that a Tier-1-only user should not have to
install. `fusion.screen_candidate` imports this lazily, inside the
`run_structure=True` branch.

See docs/tier2-tap-plabdab-plan.md for the design.
"""

from .modelling import ModellingError, build_model
from .tap_runner import TapError, TapProfile, run_tap_profile, tap_flags

__all__ = [
    "build_model",
    "ModellingError",
    "run_tap_profile",
    "tap_flags",
    "TapProfile",
    "TapError",
]
