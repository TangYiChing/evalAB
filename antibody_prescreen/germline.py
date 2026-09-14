"""Per-residue germline identity: is this residue self, or is it engineered?

This module exists because of a gap in what IEDB can tell us. The chain of
events that makes an antibody drug immunogenic is:

    APC engulfs the drug -> chops it into ~15-mers -> loads them onto MHC-II
    -> presents to CD4+ T cells -> a T cell recognises one -> B cells make
    anti-drug antibodies -> the drug is neutralised and cleared

IEDB's prediction API covers step 3 (will this peptide bind MHC-II) and its
query API records step 5 outcomes (was a T-cell response ever observed). What
neither knows is **central tolerance**: T cells reactive to self peptides are
deleted in the thymus, so a germline-encoded framework peptide can be
presented beautifully and still provoke nothing, because nobody is left to
recognise it.

Without that correction the check fires on every antibody. Measured on the
reference VH: all 8 strong predicted MHC-II binders were germline framework,
and the strongest one carries a motif present in 34% of clinical-stage
therapeutics in PLAbDab. See docs/iedb-immunogenicity-evaluation.md.

ANARCI already ships the germline reference — 249 human IGHV alleles plus the
light-chain loci — as IMGT-position-aligned strings, in the same coordinate
system numbering.py produces. So this is an index lookup, not a new data
source, a nucleotide round-trip, or a web service.
"""

from dataclasses import dataclass

from .numbering import NumberedChain

# Per-residue germline states. These are three distinct things, not two:
GERMLINE = "germline"  # matches the assigned V gene -> self, tolerised
MUTATED = "mutated"  # differs from the assigned V gene -> engineered/somatic
BEYOND_V = "beyond_v"  # CDR3/J region, no V-gene germline exists for it
UNKNOWN = "unknown"  # no germline assignment available at all

# BEYOND_V is easy to overlook and matters most: CDR3 is built by V-D-J
# recombination, so it exists in no germline V reference by construction.
# That is exactly why it is both the main antigen-binding region and the
# highest immunogenicity risk — same cause.

# ANARCI's germline tables are keyed by chain class, not by chain_type.
_CHAIN_CLASS = {"H": "H", "K": "K", "L": "L"}


@dataclass
class GermlineProfile:
    v_gene: str | None
    species: str | None
    identity: float | None
    # One entry per residue, aligned to chain.full_sequence() by index.
    status: list[str]

    @property
    def available(self) -> bool:
        return self.v_gene is not None

    def window_status(self, start: int, end: int) -> list[str]:
        """Status slice for a half-open [start, end) window of full_sequence()."""
        return self.status[start:end]

    def is_self(self, start: int, end: int) -> bool:
        """True if every residue in [start, end) is germline-encoded.

        A window that is entirely germline is a self peptide: central
        tolerance should have removed the T cells that would react to it.
        """
        window = self.window_status(start, end)
        return bool(window) and all(s == GERMLINE for s in window)


def germline_profile(chain: NumberedChain) -> GermlineProfile:
    """Classify every residue of a numbered chain as germline / mutated / beyond-V.

    Degrades to all-UNKNOWN rather than raising if ANARCI could not assign a
    germline — the caller then has to treat every peptide as potentially
    non-self, which is the conservative direction to fail in.
    """
    residues = [r for r in chain.residues if r.aa != "-"]

    if not chain.v_gene or not chain.v_species:
        return GermlineProfile(None, None, None, [UNKNOWN] * len(residues))

    try:
        from anarci import germlines as anarci_germlines

        chain_class = _CHAIN_CLASS.get(chain.chain_type, chain.chain_type)
        reference = anarci_germlines.all_germlines["V"][chain_class][chain.v_species][
            chain.v_gene
        ]
    except (ImportError, KeyError):
        return GermlineProfile(
            chain.v_gene, chain.v_species, chain.v_identity, [UNKNOWN] * len(residues)
        )

    status = []
    for residue in residues:
        index = residue.position - 1
        if index < 0 or index >= len(reference):
            # Past the end of the V gene: J/CDR3 territory.
            status.append(BEYOND_V)
            continue
        germline_aa = reference[index]
        if germline_aa == "-":
            # A gap in the germline at this IMGT position means the V gene
            # simply has no residue here — the query's residue is insertional,
            # so it is not covered by tolerance either.
            status.append(BEYOND_V)
        elif germline_aa == residue.aa:
            status.append(GERMLINE)
        else:
            status.append(MUTATED)
    return GermlineProfile(chain.v_gene, chain.v_species, chain.v_identity, status)
