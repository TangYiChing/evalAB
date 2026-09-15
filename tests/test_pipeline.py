"""Tests for the Tier 1 antibody pre-screening pipeline.

Reference sequences below are NOT verified against a canonical database
in this session (RCSB/IMGT/TheraSAbDab were unreachable from this sandbox's
network). VH_REF/VL_REF are used purely as a structurally valid human-
framework IgG Fv pair for exercising the pipeline; ABNUMBER_README_VH is a
real published example sequence (from the AbNumber project's own README).
Do not treat VH_REF/VL_REF as a confirmed match to any specific approved
therapeutic without independent verification.
"""

import pytest

from antibody_prescreen.checks import run_all_checks
from antibody_prescreen.fusion import screen_batch, screen_candidate
from antibody_prescreen.numbering import NumberingError, number_antibody, number_chain

VH_REF = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTIS"
    "ADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
VL_REF = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTD"
    "FTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)

# Real published example VH from AbNumber's README (murine).
ABNUMBER_README_VH = (
    "QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATL"
    "TTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSSAKTTAPSVYPLA"
)


class TestNumbering:
    def test_numbers_real_vh(self):
        chain = number_chain(VH_REF, expect_chain_type="H")
        assert chain.chain_type == "H"
        assert chain.species == "human"
        assert len(chain.full_sequence()) == len(VH_REF)

    def test_numbers_real_vl(self):
        chain = number_chain(VL_REF, expect_chain_type="L")
        assert chain.chain_type in ("K", "L")
        assert len(chain.full_sequence()) == len(VL_REF)

    def test_numbers_published_abnumber_vh(self):
        chain = number_chain(ABNUMBER_README_VH, expect_chain_type="H")
        assert chain.chain_type == "H"

    def test_rejects_non_antibody_sequence(self):
        with pytest.raises(NumberingError):
            number_chain("ACDEFGHIKLMNPQRSTVWY" * 3)

    def test_rejects_wrong_chain_type_expectation(self):
        with pytest.raises(NumberingError):
            number_chain(VH_REF, expect_chain_type="L")

    def test_cdr_extraction_matches_expected_region(self):
        chain = number_chain(VH_REF, expect_chain_type="H")
        # CDR3 should be a substring of the original sequence, non-empty.
        cdr3 = chain.region_sequence("CDR3")
        assert cdr3 and cdr3 in VH_REF


class TestChecks:
    def test_clean_pair_has_no_hard_gates(self):
        chains = number_antibody(VH_REF, VL_REF)
        for name, chain in chains.items():
            result = run_all_checks(chain, name)
            assert result.hard_gates == [], f"{name} unexpectedly hard-gated: {result.hard_gates}"

    def test_odd_cysteine_count_is_a_heavily_weighted_soft_flag_not_a_hard_gate(self):
        # GFNIKDTY -> GFNIKDTC: swaps a Y for an extra, unpaired C in CDR1.
        # Real solved structures (e.g. 1sy6) can have a genuine odd count
        # from a non-canonical disulfide - sequence alone can't tell that
        # apart from a real defect, so this must not auto-reject.
        broken_vh = VH_REF.replace("GFNIKDTY", "GFNIKDTC")
        chains = number_antibody(broken_vh, VL_REF)
        result = run_all_checks(chains["VH"], "VH")
        assert not any(f.check == "cysteine_pairing" for f in result.hard_gates)
        soft_cys_flags = [f for f in result.soft_flags if f.check == "cysteine_pairing"]
        assert any("odd cysteine count" in f.message for f in soft_cys_flags)

    def test_cdr_n_glycosylation_is_caught_but_does_not_hard_gate(self):
        """A CDR N-glyc motif is flagged and repairable, not a hard gate.

        This deliberately changed. It used to assert a hard gate. Measured on
        150 clinical-stage therapeutics, a CDR N-glyc motif fires on 2.7% of
        them — Ispectamab, Zanolimumab and Puxitatug all carry one. The point
        estimate clears the 5% false-rejection budget but the 95% upper bound
        is 6.7%, and the rule is that the upper bound must clear it.

        It is also one substitution away from fixed, which makes it a repair
        cost rather than a rejection.
        """
        # Insert an N-x-S/T motif into CDR-H1 (GFNIKDTY -> GFNISDTY: N-I-S).
        broken_vh = VH_REF.replace("GFNIKDTY", "GFNISDTY")
        chains = number_antibody(broken_vh, VL_REF)
        result = run_all_checks(chains["VH"], "VH")

        motif_flags = [
            f
            for f in result.flags
            if f.check == "n_glycosylation" and f.region.startswith("CDR")
        ]
        assert motif_flags, "the motif must still be detected"
        assert all(f.severity == "soft" for f in motif_flags)
        assert all(f.repairable for f in motif_flags)
        assert not result.hard_gates

    def test_ptm_liability_detects_deamidation_motif(self):
        chains = number_antibody(VH_REF, VL_REF)
        result = run_all_checks(chains["VH"], "VH")
        deamidation_flags = [
            f for f in result.soft_flags if f.check == "ptm_liability" and "deamidation" in f.message
        ]
        # VH_REF's own CDR2 naturally contains an N-G motif (position 62) —
        # this is exercising real biology, not an injected mutation.
        assert deamidation_flags

    def test_non_standard_amino_acid_hard_gates(self):
        broken_vh = VH_REF.replace("EVQL", "EVQB")  # B is not a standard AA
        chains = number_antibody(broken_vh, VL_REF)
        result = run_all_checks(chains["VH"], "VH")
        assert any(f.check == "sequence_sanity" for f in result.hard_gates)


class TestFusion:
    def test_clean_pair_is_not_hard_gated(self):
        """The reference pair must be SOFT-scored, never hard-gated.

        Asserts the hard-gate property directly rather than a routed outcome,
        which would couple this test to wherever the routing rules happen to
        sit. It did once: when thresholds still existed, moving T_HIGH from
        40.0 to 38.31 broke this test over a change that had nothing to do
        with hard gates.

        The pair landing in the worst decile is not a bug. It carries several
        CDR deamidation motifs, and (see the package README) it is not a
        verified match to any approved drug — a structurally valid human
        framework test pair, nothing more.
        """
        result = screen_candidate("test", VH_REF, VL_REF, run_immunogenicity=False)
        assert result.error is None
        assert result.score != float("inf"), "score is inf only when hard-gated"
        assert not [f for f in result.all_flags if f.severity == "hard_gate"]

    def test_odd_cysteine_does_not_auto_reject(self):
        """An odd cysteine count must NOT reject on its own.

        This deliberately changed. It used to assert NO-GO. Measured on 150
        clinical-stage therapeutics, `cysteine_pairing` fires on 5.3% of
        them — above the 5% false-rejection budget — so it is not eligible
        to be a rejecting (Level 5) gate. Sequence alone can only say "the total
        is odd",
        which is true of real approved molecules carrying a non-canonical
        disulfide.

        The question it cannot answer — *which* cysteine is unpaired — needs a
        structure, and is tested in test_structure.py.
        """
        broken_vh = VH_REF.replace("GFNIKDTY", "GFNIKDTC")
        result = screen_candidate("test", broken_vh, VL_REF, run_immunogenicity=False)

        assert result.triage_result.level != 5, "must not auto-reject"
        assert result.score != float("inf")
        assert any(f.check == "cysteine_pairing" for f in result.all_flags)

    def test_garbage_input_is_error_not_crash(self):
        result = screen_candidate("test", "NOTANANTIBODY", VL_REF, run_immunogenicity=False)
        assert result.error is not None
        # An un-numberable sequence gets no triage level at all, rather than
        # being assigned one — it was never screened, which is a different
        # thing from having been screened and rejected.
        assert result.triage_result is None

    def test_batch_screening_regression_fixture(self):
        """One candidate per triage outcome — a regression fixture: re-run
        this after any routing change and confirm the same outcomes."""
        broken_cys_vh = VH_REF.replace("GFNIKDTY", "GFNIKDTC")
        extra_cdr_liability_vh = VH_REF.replace("SRWGGDGFYAMDY", "SRWGGDGFNGMDY")

        candidates = [
            {"candidate_id": "clean", "vh_sequence": VH_REF, "vl_sequence": VL_REF},
            {
                "candidate_id": "broken_cys",
                "vh_sequence": broken_cys_vh,
                "vl_sequence": VL_REF,
            },
            {
                "candidate_id": "extra_cdr_liability",
                "vh_sequence": extra_cdr_liability_vh,
                "vl_sequence": VL_REF,
            },
        ]
        results = screen_batch(candidates, run_immunogenicity=False)
        by_id = {r.candidate_id: r for r in results}

        # No candidate here crosses a Level 5 boundary: none has a broken fold,
        # a non-standard residue, or a poly-residue CDR run. An odd cysteine
        # count is a review observation, not a rejection — see
        # test_odd_cysteine_does_not_auto_reject.
        assert all(r.triage_result.level != 5 for r in results)

        # The ordering property is still the thing to guard: the variant with
        # an extra CDR liability must carry more repair cost than the clean
        # baseline. Asserted on repair count now rather than on the summed
        # score, because nothing routes on the score any more.
        assert (
            by_id["extra_cdr_liability"].triage_result.cdr_repairs
            > by_id["clean"].triage_result.cdr_repairs
        )
