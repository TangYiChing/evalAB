"""Tier 2 (structure/TAP) tests.

Split into two halves deliberately:

  - Pure unit tests (flag mapping, cache keying, graceful degradation) run
    with no heavy dependencies and no model building. These are fast and
    always run.
  - Integration tests build a real ABodyBuilder2 model and run the real psa
    binary. They are marked `slow` and skipped when ImmuneBuilder is absent,
    so Tier-1-only checkouts still get a green suite.

The integration assertions follow the same style as the existing Tier 1
suite: take a known-clean reference pair, inject a specific liability by
mutation, and assert the check catches it.
"""

import pytest

from antibody_prescreen import format_report, screen_candidate
from antibody_prescreen.structure.modelling import ModellingError, model_key
from antibody_prescreen.structure.tap_runner import (
    AMBER_WEIGHT,
    RED_WEIGHT,
    TapProfile,
    tap_flags,
)
from tests.test_pipeline import VH_REF, VL_REF

try:
    import ImmuneBuilder  # noqa: F401

    HAVE_IMMUNEBUILDER = True
except ImportError:
    HAVE_IMMUNEBUILDER = False

needs_immunebuilder = pytest.mark.skipif(
    not HAVE_IMMUNEBUILDER, reason="ImmuneBuilder not installed (Tier 2 is optional)"
)

CDR_H3_REF = "SRWGGDGFYAMDY"


def _mutate_cdr_h3(replacement: str) -> str:
    assert CDR_H3_REF in VH_REF, "reference CDR-H3 drifted; update CDR_H3_REF"
    return VH_REF.replace(CDR_H3_REF, replacement)


# --- flag mapping (no dependencies) --------------------------------------


def test_green_metrics_produce_no_flags():
    profile = TapProfile(
        values={"SFvCSP": 9.0, "Positive Patch Score": 0.5},
        flags={"SFvCSP": "GREEN", "Positive Patch Score": "GREEN"},
    )
    assert tap_flags(profile) == []
    assert profile.worst_flag == "GREEN"
    assert profile.summary() == "all green"


def test_amber_and_red_map_to_weighted_soft_flags():
    profile = TapProfile(
        values={"SFvCSP": -10.5, "Negative Patch Score": 10.87},
        flags={"SFvCSP": "AMBER", "Negative Patch Score": "RED"},
    )
    flags = tap_flags(profile)
    weights = sorted(f.weight for f in flags)
    assert weights == sorted([AMBER_WEIGHT, RED_WEIGHT])
    assert profile.worst_flag == "RED"
    assert profile.summary() == "1R/1A"


def test_tap_never_emits_a_hard_gate():
    """A TAP RED must stay recoverable.

    Several approved therapeutics sit outside the green region on at least one
    TAP axis, so a red flag cannot mean an unconditional reject — this is the
    same lesson the odd-cysteine hard gate already taught (see README).
    """
    profile = TapProfile(
        values={m: 99.0 for m in ("SFvCSP", "Positive Patch Score")},
        flags={m: "RED" for m in ("SFvCSP", "Positive Patch Score")},
    )
    assert all(f.severity == "soft" for f in tap_flags(profile))


# --- model cache keying ---------------------------------------------------


def test_model_key_is_content_addressed_and_normalised():
    assert model_key(VH_REF, VL_REF) == model_key(
        VH_REF.lower(), f"  {VL_REF}  "
    ), "cache key should ignore case and surrounding whitespace"
    assert model_key(VH_REF, VL_REF) != model_key(VL_REF, VH_REF)


# --- graceful degradation -------------------------------------------------


def test_structure_failure_degrades_to_tier1(monkeypatch):
    """A Tier 2 failure must not fail the candidate or the batch."""

    def boom(*args, **kwargs):
        raise ModellingError("simulated ABodyBuilder2 failure")

    monkeypatch.setattr("antibody_prescreen.structure.build_model", boom)

    result = screen_candidate(
        "degraded", VH_REF, VL_REF, run_immunogenicity=False, run_structure=True
    )

    assert result.verdict != "ERROR"
    assert result.structure_available is False
    assert result.tap_profile is None
    reasons = [f.message for f in result.all_flags if f.check == "tap"]
    assert any("Tier 2 unavailable" in m for m in reasons), (
        "the reason a candidate silently fell back to Tier 1 must be recorded"
    )


def test_report_omits_tap_column_without_structure():
    result = screen_candidate(
        "tier1only", VH_REF, VL_REF, run_immunogenicity=False, run_structure=False
    )
    report = format_report([result])
    assert "| TAP |" not in report
    assert result.structure_available is False


# --- integration: real model, real psa ------------------------------------


@pytest.mark.slow
@needs_immunebuilder
def test_reference_pair_is_all_green(tmp_path):
    """The clean human-framework reference pair should pass all five metrics."""
    result = screen_candidate(
        "ref",
        VH_REF,
        VL_REF,
        run_immunogenicity=False,
        run_structure=True,
        model_cache=tmp_path,
    )
    assert result.structure_available is True
    assert set(result.tap_profile.flags) == {
        "Hydrophobic Patch Score",
        "Negative Patch Score",
        "Positive Patch Score",
        "SFvCSP",
        "Total IMGT CDR Length",
    }
    assert result.tap_profile.worst_flag == "GREEN"
    assert [f for f in result.all_flags if f.check == "tap"] == []


@pytest.mark.slow
@needs_immunebuilder
def test_injected_positive_patch_is_caught(tmp_path):
    """Replacing CDR-H3 with poly-arginine must drive PPC red.

    Measured: PPC goes 0.59 (green) -> 9.47, against an amber ceiling of 3.58.
    """
    result = screen_candidate(
        "polyR",
        _mutate_cdr_h3("R" * len(CDR_H3_REF)),
        VL_REF,
        run_immunogenicity=False,
        run_structure=True,
        model_cache=tmp_path,
    )
    assert result.structure_available is True
    assert result.tap_profile.flags["Positive Patch Score"] == "RED"
    assert result.tap_profile.values["Positive Patch Score"] > 3.58


@pytest.mark.slow
@needs_immunebuilder
def test_injected_negative_patch_is_caught(tmp_path):
    """Replacing CDR-H3 with poly-aspartate must drive PNC red.

    Measured: PNC goes 0.09 (green) -> 10.87, against an amber ceiling of 3.50.
    It also drags SFvCSP from 9.0 to -10.5 (amber) — charge imbalance is
    exactly what SFvCSP is for, so both moving is the expected behaviour.
    """
    result = screen_candidate(
        "polyD",
        _mutate_cdr_h3("D" * len(CDR_H3_REF)),
        VL_REF,
        run_immunogenicity=False,
        run_structure=True,
        model_cache=tmp_path,
    )
    assert result.tap_profile.flags["Negative Patch Score"] == "RED"
    assert result.tap_profile.values["Negative Patch Score"] > 3.50
    assert result.tap_profile.flags["SFvCSP"] == "AMBER"

    tap_messages = [f.message for f in result.all_flags if f.check == "tap"]
    assert len(tap_messages) == 2
    assert format_report([result]).count("| TAP |") == 1


@pytest.mark.slow
@needs_immunebuilder
def test_model_cache_is_reused(tmp_path):
    from antibody_prescreen.structure import build_model

    first = build_model(VH_REF, VL_REF, cache_dir=tmp_path)
    mtime = first.stat().st_mtime_ns
    second = build_model(VH_REF, VL_REF, cache_dir=tmp_path)
    assert first == second
    assert second.stat().st_mtime_ns == mtime, "cached model was rebuilt"
