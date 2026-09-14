"""IMGT numbering and CDR/framework boundary assignment, via ANARCI.

Every other Tier 1 check (PTM liability, N-glyc, immunogenicity, aggregation
region weighting) needs to know which residues sit in a CDR vs a framework
region. This module is the single place that answers that question, so a
future swap to AbNumber or a different numbering scheme only touches this
file.
"""

from dataclasses import dataclass

from anarci import run_anarci

# Standard IMGT V-domain region boundaries (IMGT unique numbering, 1-128).
# Same boundaries apply to heavy and light chains under the IMGT scheme.
IMGT_REGIONS = [
    ("FR1", 1, 26),
    ("CDR1", 27, 38),
    ("FR2", 39, 55),
    ("CDR2", 56, 65),
    ("FR3", 66, 104),
    ("CDR3", 105, 117),
    ("FR4", 118, 128),
]


class NumberingError(Exception):
    """Raised when ANARCI can't number a sequence (not an antibody V-domain, wrong species, etc)."""


@dataclass
class Residue:
    position: int
    insertion: str
    aa: str
    region: str  # one of IMGT_REGIONS names
    is_cdr: bool


@dataclass
class NumberedChain:
    chain_type: str  # "H", "K", or "L"
    species: str
    evalue: float
    residues: list[Residue]
    # Closest germline V gene, from ANARCI's own assignment. None if ANARCI
    # could not assign one (e.g. an unusual species). Used by the germline
    # module to decide which residues are self and which are engineered.
    v_gene: str | None = None
    v_species: str | None = None
    v_identity: float | None = None

    def region_sequence(self, region: str) -> str:
        return "".join(r.aa for r in self.residues if r.region == region and r.aa != "-")

    def full_sequence(self) -> str:
        return "".join(r.aa for r in self.residues if r.aa != "-")

    def cdr_positions(self) -> set[int]:
        """0-based indices into full_sequence() that fall in a CDR loop."""
        positions = set()
        idx = 0
        for r in self.residues:
            if r.aa == "-":
                continue
            if r.is_cdr:
                positions.add(idx)
            idx += 1
        return positions


def _region_for_position(position: int) -> str:
    for name, start, end in IMGT_REGIONS:
        if start <= position <= end:
            return name
    # ANARCI can number a few residues past 128 in rare long C-terminal
    # stretches; treat anything past the defined range as FR4 rather than
    # erroring out, since it's after CDR3 in the sequence.
    return "FR4"


def number_chain(sequence: str, expect_chain_type: str | None = None) -> NumberedChain:
    """Number one VH or VL sequence and assign IMGT CDR/framework regions.

    expect_chain_type: optional "H" or "L" (K and L are both "light" for our
    purposes) to sanity-check ANARCI's own chain-type call against what the
    caller expected (e.g. caller submitted this as a VL sequence).
    """
    if not sequence or not sequence.isalpha():
        raise NumberingError(f"not a valid amino acid sequence: {sequence!r}")

    # assign_germline=True makes ANARCI also report the closest germline V/J
    # gene. It costs one extra alignment inside the same call — cheaper than
    # running ANARCI twice, and immunogenicity.py needs it to tell a
    # self/tolerised residue apart from an engineered one.
    _, numbered, alignment_details, _ = run_anarci(
        [("query", sequence)], scheme="imgt", assign_germline=True
    )

    domains = numbered[0]
    if not domains:
        raise NumberingError(
            "ANARCI found no antibody V-domain in this sequence — check it's a real VH/VL "
            "sequence and not truncated or containing non-standard residues"
        )
    if len(domains) > 1:
        raise NumberingError(
            f"ANARCI found {len(domains)} V-domains in one sequence — expected exactly one "
            "VH or VL per input; check for a scFv/tandem construct that needs splitting first"
        )

    numbered_residues, _domain_start, _domain_end = domains[0]
    details = alignment_details[0][0]
    chain_type = details["chain_type"]

    if expect_chain_type is not None:
        got_is_heavy = chain_type == "H"
        expect_is_heavy = expect_chain_type == "H"
        if got_is_heavy != expect_is_heavy:
            raise NumberingError(
                f"expected a {'heavy' if expect_is_heavy else 'light'} chain but ANARCI "
                f"identified this as chain_type={chain_type!r} — check the sequence is in "
                "the VH/VL slot you think it is"
            )

    residues = []
    for (position, insertion), aa in numbered_residues:
        region = _region_for_position(position)
        residues.append(
            Residue(
                position=position,
                insertion=insertion.strip(),
                aa=aa,
                region=region,
                is_cdr=region.startswith("CDR"),
            )
        )

    v_gene = v_species = v_identity = None
    germline_call = (details.get("germlines") or {}).get("v_gene")
    if germline_call:
        (v_species, v_gene), v_identity = germline_call

    return NumberedChain(
        chain_type=chain_type,
        species=details["species"],
        evalue=details["evalue"],
        residues=residues,
        v_gene=v_gene,
        v_species=v_species,
        v_identity=v_identity,
    )


def number_antibody(vh_sequence: str, vl_sequence: str) -> dict[str, NumberedChain]:
    """Number a VH/VL pair. Raises NumberingError on either failing."""
    return {
        "VH": number_chain(vh_sequence, expect_chain_type="H"),
        "VL": number_chain(vl_sequence, expect_chain_type="L"),
    }
