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

# Tolerance belongs to the PATIENT, not to the molecule.
#
# This is the subtle part, and getting it wrong inverts the check. ANARCI
# assigns each chain its own closest germline — for a murine antibody that is
# a mouse V gene. Comparing a murine framework against mouse germline marks it
# "self" and drops its peptides, but a murine framework is emphatically NOT
# self to a human patient: that is the classic HAMA failure mode.
#
# Measured on murine antibodies from PLAbDab, strong binders dropped as self:
#
#   ID                VH v_gene        strong   vs mouse ref   vs human ref
#   UCY74699_UCY74691 IGHV5-9-3*01         10              6              0
#   AVY28417_AVY28469 IGHV2-2*03           11              9              0
#   QRJ72240_QRJ72247 IGHV1-71-16*01        7              5              0
#
# So the reference species is a property of who is being treated, and it
# defaults to human. A non-human framework then correctly shows up as almost
# entirely non-self, which is the signal we want, not noise to suppress.
REFERENCE_SPECIES = "human"


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


def _closest_v_gene(
    chain: NumberedChain, species: str, table: dict
) -> tuple[str, float] | None:
    """Closest V gene to this chain WITHIN `species`, by IMGT-aligned identity.

    Needed when the chain's own assigned germline is a different species from
    the patient's: a murine antibody has no mouse-independent germline call,
    so we ask "what human V gene is this closest to" and measure divergence
    from that.
    """
    residues = [r for r in chain.residues if r.aa != "-"]
    best_gene, best_identity = None, -1.0
    for gene, reference in table.items():
        matched = compared = 0
        for residue in residues:
            index = residue.position - 1
            if index < len(reference) and reference[index] != "-":
                compared += 1
                matched += reference[index] == residue.aa
        if compared:
            identity = matched / compared
            if identity > best_identity:
                best_gene, best_identity = gene, identity
    return (best_gene, best_identity) if best_gene else None


def germline_profile(
    chain: NumberedChain, reference_species: str = REFERENCE_SPECIES
) -> GermlineProfile:
    """Classify every residue as germline / mutated / beyond-V.

    Args:
        reference_species: whose tolerance is being modelled — the PATIENT's
            species, not the antibody's. Defaults to human. See the comment on
            REFERENCE_SPECIES for why this distinction inverts the check if
            got wrong.

    Degrades to all-UNKNOWN rather than raising if no germline can be
    assigned — the caller then has to treat every peptide as potentially
    non-self, which is the conservative direction to fail in.
    """
    residues = [r for r in chain.residues if r.aa != "-"]

    if not chain.v_gene or not chain.v_species:
        return GermlineProfile(None, None, None, [UNKNOWN] * len(residues))

    try:
        from anarci import germlines as anarci_germlines

        chain_class = _CHAIN_CLASS.get(chain.chain_type, chain.chain_type)
        species_table = anarci_germlines.all_germlines["V"][chain_class][
            reference_species
        ]
    except (ImportError, KeyError):
        return GermlineProfile(
            chain.v_gene, chain.v_species, chain.v_identity, [UNKNOWN] * len(residues)
        )

    if chain.v_species == reference_species:
        v_gene, v_identity = chain.v_gene, chain.v_identity
    else:
        # Non-matching species: re-anchor onto the closest gene of the
        # patient's species. Everything that differs from it is then correctly
        # treated as non-self.
        closest = _closest_v_gene(chain, reference_species, species_table)
        if closest is None:
            return GermlineProfile(
                chain.v_gene, chain.v_species, chain.v_identity,
                [UNKNOWN] * len(residues),
            )
        v_gene, v_identity = closest

    reference = species_table.get(v_gene)
    if reference is None:
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
    return GermlineProfile(v_gene, reference_species, v_identity, status)
