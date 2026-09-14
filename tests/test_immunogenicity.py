"""Immunogenicity and germline tests.

The point of these tests is the germline correction, not the API plumbing.
Without it the check fires on every antibody: measured on the reference VH,
all 8 strong predicted MHC-II binders are germline framework, and the
strongest carries a motif present in 34% of clinical-stage therapeutics.

Network tests are marked `network` and skipped by default so the suite stays
offline and fast. Run them with `-m network`.
"""

import pytest

from antibody_prescreen.checks import check_v_domain_integrity
from antibody_prescreen.germline import (
    BEYOND_V,
    GERMLINE,
    MUTATED,
    germline_profile,
)
from antibody_prescreen.immunogenicity import check_immunogenicity
from antibody_prescreen.numbering import number_antibody
from tests.test_pipeline import VH_REF, VL_REF

# Reference VH with influenza HA306-318 (PKYVKQNTLKLAT) grafted into CDR-H3.
# The textbook universal MHC-II epitope, used as the positive control: a
# filter that drops everything would pass the negative test but fail this one.
VH_FLU_EPITOPE = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTIS"
    "ADTSKNTAYLQMNSLRAEDTAVYYCSRPKYVKQNTLKLATWGQGTLVTVSS"
)

needs_network = pytest.mark.network


# --- germline assignment --------------------------------------------------


def test_germline_assignment_is_available():
    chains = number_antibody(VH_REF, VL_REF)
    assert chains["VH"].v_gene == "IGHV3-66*01"
    assert chains["VH"].v_species == "human"
    assert chains["VL"].v_gene == "IGKV1-39*01"


def test_germline_profile_has_all_three_states():
    """germline / mutated / beyond_v are three distinct things, not two.

    beyond_v is the one that is easy to miss and matters most: CDR3 is built
    by V-D-J recombination, so no germline V reference contains it — which is
    exactly why it carries no tolerance protection.
    """
    chain = number_antibody(VH_REF, VL_REF)["VH"]
    profile = germline_profile(chain)

    assert profile.available
    assert len(profile.status) == len(chain.full_sequence())
    assert set(profile.status) == {GERMLINE, MUTATED, BEYOND_V}

    # CDR3 is only PARTLY junctional. IMGT 105-106 are the C-A-R stem, which
    # the V gene does encode; everything from 107 on is V-D-J junction and
    # exists in no germline V reference. Measured on the reference VH:
    #   105 S mutated / 106 R germline / 107-117 all beyond_v
    # So the meaningful assertion is about the junction, not all of CDR3.
    residues = [r for r in chain.residues if r.aa != "-"]
    junction = [
        i for i, r in enumerate(residues) if r.region == "CDR3" and r.position >= 107
    ]
    assert junction
    assert all(profile.status[i] == BEYOND_V for i in junction)

    stem = [i for i, r in enumerate(residues) if r.region == "CDR3" and r.position < 107]
    assert any(profile.status[i] != BEYOND_V for i in stem), (
        "the C-A-R stem is V-gene encoded and must not be treated as junctional"
    )


def test_framework_3_core_is_recognised_as_self():
    """YLQMNSLRA is the germline core behind every false positive we measured."""
    chain = number_antibody(VH_REF, VL_REF)["VH"]
    profile = germline_profile(chain)
    seq = chain.full_sequence()

    start = seq.find("YLQMNSLRA")
    assert start != -1
    assert profile.is_self(start, start + len("YLQMNSLRA"))


def test_grafted_epitope_core_is_not_self():
    chain = number_antibody(VH_FLU_EPITOPE, VL_REF)["VH"]
    profile = germline_profile(chain)
    seq = chain.full_sequence()

    start = seq.find("YVKQNTLKL")
    assert start != -1
    assert not profile.is_self(start, start + len("YVKQNTLKL"))


def test_profile_degrades_without_a_germline_call():
    chain = number_antibody(VH_REF, VL_REF)["VH"]
    chain.v_gene = None
    profile = germline_profile(chain)

    assert not profile.available
    # is_self must be False when germline is unknown, so an unknown chain
    # fails toward reporting binders rather than silently dropping them.
    assert not profile.is_self(0, 9)


# --- V-domain integrity (the AA-level answer to "productivity") ------------


def test_reference_pair_has_intact_anchors():
    for name, chain in number_antibody(VH_REF, VL_REF).items():
        assert check_v_domain_integrity(chain, name) == []


@pytest.mark.parametrize(
    "imgt_position,replacement",
    [(23, "S"), (41, "R"), (104, "A")],
)
def test_broken_anchor_is_caught(imgt_position, replacement):
    chain = number_antibody(VH_REF, VL_REF)["VH"]
    residues = [r for r in chain.residues if r.aa != "-"]
    index = next(i for i, r in enumerate(residues) if r.position == imgt_position)

    mutated = VH_REF[:index] + replacement + VH_REF[index + 1 :]
    flags = check_v_domain_integrity(number_antibody(mutated, VL_REF)["VH"], "VH")

    assert len(flags) == 1
    assert str(imgt_position) in flags[0].message
    assert flags[0].severity == "soft"  # never a hard gate, per the cysteine precedent


def test_truncated_v_domain_is_caught():
    """Truncation before the J-PHE/J-TRP is the most common real failure."""
    flags = check_v_domain_integrity(number_antibody(VH_REF[:-12], VL_REF)["VH"], "VH")
    assert len(flags) == 1
    assert "118" in flags[0].message
    assert "truncated" in flags[0].message


# --- the check itself (network) -------------------------------------------


@needs_network
def test_germline_framework_binders_are_all_dropped(tmp_path):
    """The negative control: a clean human antibody must produce no flags.

    Measured: 8 strong binders predicted, all sharing the germline core
    YLQMNSLRA, all dropped.
    """
    chain = number_antibody(VH_REF, VL_REF)["VH"]
    result = check_immunogenicity(chain, "VH", cache_dir=tmp_path)

    assert result.available
    assert result.n_predicted_total > 0, "expected the predictor to return binders"
    assert result.n_dropped_germline == result.n_predicted_total
    assert result.binders == []
    assert result.flags == []


@needs_network
def test_grafted_epitope_survives_the_filter(tmp_path):
    """The positive control: a filter that drops everything fails this."""
    chain = number_antibody(VH_FLU_EPITOPE, VL_REF)["VH"]
    result = check_immunogenicity(chain, "VH", cache_dir=tmp_path)

    assert result.available
    assert result.n_dropped_germline > 0, "germline framework should still be dropped"
    assert result.binders, "the grafted epitope must survive"
    assert {b.core for b in result.binders} == {"YVKQNTLKL"}

    predicted = [f for f in result.flags if f.check == "immunogenicity"]
    # One flag per distinct core, not per allele or per overlapping window —
    # the same core appears in several rows and must not be weighted twice.
    assert len(predicted) == 1
    assert "overlaps CDR" in predicted[0].message


@needs_network
def test_observed_evidence_is_reported_separately(tmp_path):
    """"Predicted to bind" and "observed to provoke" are different claims."""
    chain = number_antibody(VH_FLU_EPITOPE, VL_REF)["VH"]
    result = check_immunogenicity(chain, "VH", cache_dir=tmp_path)

    observed = [f for f in result.flags if f.check == "immunogenicity_observed"]
    assert len(observed) == 1
    assert "YVKQNTLKL" in observed[0].message
    assert "observed, not predicted" in observed[0].message


@needs_network
def test_unreachable_api_degrades_gracefully(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "antibody_prescreen.immunogenicity.IEDB_MHCII_URL",
        "https://invalid.iedb.example/tools_api/mhcii/",
    )
    chain = number_antibody(VH_REF, VL_REF)["VH"]
    result = check_immunogenicity(chain, "VH", cache_dir=tmp_path)

    assert result.available is False
    assert result.flags == []
    assert "unreachable" in result.note
