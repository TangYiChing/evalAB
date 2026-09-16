"""Tests for the range-based triage: profile -> bands -> level matrix.

These are unit tests of the routing logic, deliberately not of the biology.
A test that asserts "this antibody should be Level 2" would encode a judgement
nobody in this repository is able to make; a test that asserts "deviation 2 at
redesign cost routes to Level 4" encodes the table, which is a decision that
was written down and can be checked.
"""

import pytest

from antibody_prescreen.bands import (
    CLINICAL_METRICS, Band, BandSet, compose, fit_bands, supported_tier,
)
from antibody_prescreen.levels import MATRIX, Finding, assign
from antibody_prescreen.profile import (
    COST_CDR_POINT, COST_FR_POINT, COST_NO_PATH, COST_REDESIGN,
    BinaryObservation, Profile, RangeMetric,
)

REFERENCE = "test"


def _bandset(**bands):
    return BandSet(REFERENCE, 1000, 1000, (5, 95, 99, 99.9), bands)


def _band(metric, lo, hi, n=10000, tail="both"):
    """A band with evenly spaced cut points, so each tier is reachable."""
    span = hi - lo
    return Band(metric, tail,
                low_extreme=lo - 2 * span, low_strong=lo - span, low_mild=lo,
                high_mild=hi, high_strong=hi + span, high_extreme=hi + 2 * span,
                n=n, max_tier=supported_tier(n))


def test_tier_boundaries_are_two_sided():
    b = _band("m", 10, 20)
    assert b.tier(15) == 0
    assert b.tier(25) == 1 and b.tier(5) == 1
    assert b.tier(35) == 2 and b.tier(-5) == 2
    assert b.tier(45) == 3 and b.tier(-15) == 3


def test_high_tail_metric_ignores_the_low_side():
    b = _band("count", 0, 3, tail="high")
    assert b.tier(-100) == 0
    assert b.tier(4) == 1


def test_sample_size_caps_the_tier():
    """The bug that put 23/23 designs at Level 5 on a degenerate lower tail."""
    assert supported_tier(10_000) == 3
    assert supported_tier(534) == 2      # TheraSAbDab: cannot assert 1-in-1000
    assert supported_tier(120) == 1
    assert supported_tier(50) == 0
    small = _band("m", 10, 20, n=534)
    assert small.raw_tier(100) == 3      # it did cross the p99.9 cut
    assert small.tier(100) == 2          # but 534 antibodies cannot say so


def test_degenerate_tail_is_announced_not_hidden():
    b = Band("cdrl2_len", "both", 3, 3, 3, 3, 3, 3, n=534, max_tier=2)
    assert b.degenerate_side(2) and b.degenerate_side(4)
    assert "never seen" in b.describe(4)


def test_fit_refuses_a_population_too_small_to_fit():
    bs = fit_bands({"m": list(range(199))}, {}, REFERENCE, 199, 0)
    assert "m" not in bs.bands
    bs = fit_bands({"m": list(range(500))}, {}, REFERENCE, 500, 0)
    assert "m" in bs.bands


def test_unbanded_metric_is_unmeasured_not_typical():
    """None must never be read as tier 0 — see BandSet.tier."""
    p = Profile("c", metrics=[RangeMetric("nobody_fitted_this", 42.0, COST_REDESIGN)])
    t = assign(p, _bandset())
    assert t.level == 1
    assert t.findings == []
    assert t.unmeasured and "nobody_fitted_this" in t.unmeasured[0]


@pytest.mark.parametrize("deviation,cost,expected", [
    (0, COST_FR_POINT, 1),     # ordinary framework PTM motif — free to fix
    (0, COST_CDR_POINT, 2),    # ordinary CDR motif — costs an affinity re-measure
    (1, COST_REDESIGN, 3),
    (2, COST_REDESIGN, 4),
    (2, COST_NO_PATH, 5),
    (3, COST_FR_POINT, 3),     # far out of range but one substitution back
    (3, COST_REDESIGN, 5),
])
def test_matrix_routes_by_the_two_axes(deviation, cost, expected):
    p = Profile("c", observations=[
        BinaryObservation("o", deviation, cost, "chain", "msg")])
    assert assign(p, _bandset()).level == expected


def test_matrix_is_monotone_in_both_axes():
    """A finding cannot become less serious by getting further out, or dearer."""
    for d in range(4):
        row = [MATRIX[d][c] for c in range(5)]
        assert row == sorted(row), f"deviation {d} is not monotone in cost"
    for c in range(5):
        col = [MATRIX[d][c] for d in range(4)]
        assert col == sorted(col), f"cost {c} is not monotone in deviation"


def test_level_is_the_worst_finding_not_the_sum():
    two_twos = Profile("c", observations=[
        BinaryObservation("a", 0, COST_CDR_POINT, "CDR1", "m"),
        BinaryObservation("b", 0, COST_CDR_POINT, "CDR2", "m"),
        BinaryObservation("c", 0, COST_CDR_POINT, "CDR3", "m"),
    ])
    assert assign(two_twos, _bandset()).level == 2

    one_four = Profile("c", observations=[
        BinaryObservation("a", 2, COST_REDESIGN, "chain", "m")])
    assert assign(one_four, _bandset()).level == 4


def test_drivers_are_only_the_findings_that_set_the_level():
    p = Profile("c", observations=[
        BinaryObservation("minor", 0, COST_FR_POINT, "FR1", "m"),
        BinaryObservation("major", 2, COST_REDESIGN, "chain", "m"),
    ])
    t = assign(p, _bandset())
    assert t.level == 4
    assert [f.source for f in t.drivers] == ["major"]


def test_repair_counts_are_reported_separately_from_the_level():
    p = Profile("c", observations=[
        BinaryObservation("a", 0, COST_FR_POINT, "FR1", "m"),
        BinaryObservation("b", 0, COST_FR_POINT, "FR2", "m"),
        BinaryObservation("c", 0, COST_CDR_POINT, "CDR3", "m"),
    ])
    counts = assign(p, _bandset()).repair_counts()
    assert counts["framework_point"] == 2
    assert counts["cdr_point"] == 1
    assert counts["redesign"] == 0


def test_a_clean_candidate_is_level_one():
    bs = _bandset(m=_band("m", 10, 20))
    p = Profile("c", metrics=[RangeMetric("m", 15.0, COST_REDESIGN)])
    t = assign(p, bs)
    assert t.level == 1 and t.findings == []


def test_error_profile_short_circuits():
    t = assign(Profile("c", error="ANARCI could not number this"), _bandset())
    assert t.level == 0 and t.error


def test_compose_takes_named_metrics_and_keeps_provenance():
    base = _bandset(vh_germline_identity=_band("vh_germline_identity", 0.72, 0.99, n=30000),
                    fv_pi=_band("fv_pi", 4.8, 9.2, n=30000))
    base.bands["vh_germline_identity"].source = "repertoire"
    base.bands["fv_pi"].source = "repertoire"
    clinical = _bandset(vh_germline_identity=_band("vh_germline_identity", 0.76, 0.99, n=393))
    clinical.bands["vh_germline_identity"].source = "clinical"

    out = compose(base, clinical, CLINICAL_METRICS)
    assert out.bands["vh_germline_identity"].source == "clinical"
    assert out.bands["vh_germline_identity"].low_mild == 0.76   # the tighter one
    assert out.bands["fv_pi"].source == "repertoire"            # untouched
    assert "vh_germline_identity" in out.reference              # the mix is named


def test_the_repertoire_reference_widens_a_clinical_excursion_away():
    """A CDR-H3 of 20 is an excursion clinically and ordinary in the repertoire."""
    clinical = _bandset(cdrh3_len=_band("cdrh3_len", 8, 18, n=534))
    p = Profile("c", metrics=[RangeMetric("cdrh3_len", 20.0, COST_REDESIGN)])
    assert assign(p, clinical).level > 1

    repertoire = _bandset(cdrh3_len=_band("cdrh3_len", 8, 21, n=30000))
    assert assign(p, repertoire).level == 1


def test_only_germline_identity_is_pinned_to_the_clinical_population():
    """Everything else moved to the larger, cleaned repertoire reference."""
    assert set(CLINICAL_METRICS) == {"vh_germline_identity", "vl_germline_identity"}


def test_non_human_candidate_is_flagged_as_out_of_scope():
    """Declared scope is human paired VH/VL. A murine chain is not 'unusual',
    it is a candidate nothing in the band set was fitted to describe."""
    from antibody_prescreen.numbering import number_chain
    from antibody_prescreen.profile import build_profile

    # A murine VH (anti-lysozyme D1.3) and a human VL.
    murine_vh = ("QVQLQESGPGLVAPSQSLSITCTVSGFSLTGYGVNWVRQPPGKGLEWLGMIWGDGNTDYNSA"
                 "LKSRLSISKDNSKSQVFLKMNSLHTDDTARYYCARERDYRLDYWGQGTTLTVSS")
    human_vl = ("QSVLTQPPSASGTPGQSVTISCSGSRSNIGGNTVNWYQHLPGMAPKLLIYSSNQRSSGVPDRF"
                "SGSKSGTSASLAISGLQSEDDADYYCASWDDSLNGVVFGGGTKLTVL")
    p = build_profile("murine", {"VH": number_chain(murine_vh, "H"),
                                 "VL": number_chain(human_vl, "L")})
    scope = [o for o in p.observations if o.name == "out_of_scope_species"]
    assert len(scope) == 1 and "not human" in scope[0].message
    assert assign(p, _bandset()).level == 4
