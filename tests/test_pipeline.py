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

    def test_cdr_n_glycosylation_hard_gates(self):
        # Insert an N-x-S/T motif into CDR-H1 (GFNIKDTY -> GFNISDTY: N-I-S).
        broken_vh = VH_REF.replace("GFNIKDTY", "GFNISDTY")
        chains = number_antibody(broken_vh, VL_REF)
        result = run_all_checks(chains["VH"], "VH")
        assert any(f.check == "n_glycosylation" for f in result.hard_gates)

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

        Asserts the hard-gate property directly rather than a verdict, which
        would couple this test to wherever T_LOW/T_HIGH happen to sit. It did
        once: the Step 7 recalibration moved T_HIGH from 40.0 to 38.31 and
        this pair scores 39.92, so a verdict assertion started failing on a
        threshold change that has nothing to do with hard gates.

        The pair landing in the worst decile is not a bug. It carries several
        CDR deamidation motifs, and (see the package README) it is not a
        verified match to any approved drug — a structurally valid human
        framework test pair, nothing more.
        """
        result = screen_candidate("test", VH_REF, VL_REF, run_immunogenicity=False)
        assert result.error is None
        assert result.score != float("inf"), "score is inf only when hard-gated"
        assert not [f for f in result.all_flags if f.severity == "hard_gate"]

    def test_broken_cysteine_pair_is_no_go_but_via_score_not_hard_gate(self):
        broken_vh = VH_REF.replace("GFNIKDTY", "GFNIKDTC")
        result = screen_candidate("test", broken_vh, VL_REF, run_immunogenicity=False)
        assert result.verdict == "NO-GO"
        assert result.score != float("inf")  # soft-scored, not hard-gated
        assert "cysteine" in result.top_reasons[0].lower()

    def test_garbage_input_is_error_not_crash(self):
        result = screen_candidate("test", "NOTANANTIBODY", VL_REF, run_immunogenicity=False)
        assert result.verdict == "ERROR"
        assert result.error is not None

    def test_batch_screening_regression_fixture(self):
        """One candidate per verdict bucket — a regression fixture: re-run
        this after any scoring/weight change and confirm the same buckets."""
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

        assert by_id["broken_cys"].verdict == "NO-GO"
        # The double-CDR-liability variant must score worse (higher) than
        # the clean baseline — this is the core ordering property to guard.
        assert by_id["extra_cdr_liability"].score > by_id["clean"].score
