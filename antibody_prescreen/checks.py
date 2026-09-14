"""Tier 1 sequence-based liability checks.

Each check function takes a NumberedChain (or a VH/VL pair of them) and
returns a Flag list. A Flag is either a hard gate (short-circuits the whole
candidate to NO-GO, see fusion.py) or a soft, weighted signal.

Nothing here does network calls or needs a GPU — that's the point of Tier 1.
"""

from dataclasses import dataclass, field

from .numbering import NumberedChain

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

MIN_CDR_H3_LEN = 4
MAX_CDR_H3_LEN = 30


@dataclass
class Flag:
    check: str
    severity: str  # "hard_gate" or "soft"
    region: str  # e.g. "CDR3", "FR2", or "chain" for whole-chain issues
    message: str
    weight: float = 0.0  # only meaningful for soft flags; hard gates ignore this


@dataclass
class ChainCheckResult:
    chain_name: str  # "VH" or "VL"
    flags: list[Flag] = field(default_factory=list)

    @property
    def hard_gates(self) -> list[Flag]:
        return [f for f in self.flags if f.severity == "hard_gate"]

    @property
    def soft_flags(self) -> list[Flag]:
        return [f for f in self.flags if f.severity == "soft"]


def check_sequence_sanity(chain: NumberedChain, chain_name: str) -> list[Flag]:
    """Non-standard residues, and CDR-H3 length sanity (heavy chain only)."""
    flags = []
    seq = chain.full_sequence()

    bad_residues = sorted(set(seq) - STANDARD_AA)
    if bad_residues:
        flags.append(
            Flag(
                check="sequence_sanity",
                severity="hard_gate",
                region="chain",
                message=f"non-standard amino acid(s) found: {bad_residues}",
            )
        )

    if chain_name == "VH":
        cdr3_len = len(chain.region_sequence("CDR3"))
        if cdr3_len < MIN_CDR_H3_LEN or cdr3_len > MAX_CDR_H3_LEN:
            flags.append(
                Flag(
                    check="sequence_sanity",
                    severity="hard_gate",
                    region="CDR3",
                    message=(
                        f"CDR-H3 length {cdr3_len} outside expected range "
                        f"[{MIN_CDR_H3_LEN}, {MAX_CDR_H3_LEN}] — likely a design/generation artifact"
                    ),
                )
            )

    return flags


CYSTEINE_ODD_COUNT_WEIGHT = 8.0


def check_cysteine_pairing(chain: NumberedChain, chain_name: str) -> list[Flag]:
    """Flag odd cysteine counts. Real pairing (which Cys bonds to which) needs
    a structure — sequence alone can only say "an odd count means at least
    one cysteine cannot be canonically paired," not whether that's a real
    defect. This was a hard gate initially, but real solved structures (e.g.
    a genuine non-canonical CDR-H3 disulfide) showed that's too blunt: an
    odd count is real biology often enough that auto-rejecting on it alone
    throws away legitimate candidates. Downgraded to a heavily-weighted soft
    flag — real pairing confirmation is deferred to a future Tier 2
    structure-based check, which can actually resolve which cysteines are
    close enough in 3D to bond.
    """
    flags = []
    seq = chain.full_sequence()
    cys_count = seq.count("C")

    if cys_count % 2 != 0:
        flags.append(
            Flag(
                check="cysteine_pairing",
                severity="soft",
                region="chain",
                weight=CYSTEINE_ODD_COUNT_WEIGHT,
                message=(
                    f"{chain_name} has an odd cysteine count ({cys_count}) — at least one "
                    "cysteine cannot be canonically disulfide-paired; sequence alone can't "
                    "confirm whether this is a real defect or a non-canonical disulfide, "
                    "needs structural confirmation before this is disqualifying"
                ),
            )
        )

    # Flag any cysteine sitting inside a CDR loop as a soft flag even when the
    # total count is even — an unexpected CDR cysteine is unusual and worth a
    # human's attention even if it happens to pair with something.
    for r in chain.residues:
        if r.aa == "C" and r.is_cdr:
            flags.append(
                Flag(
                    check="cysteine_pairing",
                    severity="soft",
                    region=r.region,
                    weight=2.0,
                    message=f"cysteine at {chain_name} {r.region} position {r.position} — unusual, verify intentional",
                )
            )

    return flags


def check_n_glycosylation(chain: NumberedChain, chain_name: str) -> list[Flag]:
    """N-x-[S/T] motif, x != P. Hard gate if the N is in a CDR (will be
    glycosylated in mammalian expression with unpredictable effect on
    binding), soft flag if in framework (cheap to fix, lower urgency).
    """
    flags = []
    seq = chain.full_sequence()
    residues = [r for r in chain.residues if r.aa != "-"]

    for i in range(len(seq) - 2):
        n, x, s_t = seq[i], seq[i + 1], seq[i + 2]
        if n == "N" and x != "P" and s_t in ("S", "T"):
            region = residues[i].region
            position = residues[i].position
            motif = seq[i : i + 3]
            if region.startswith("CDR"):
                flags.append(
                    Flag(
                        check="n_glycosylation",
                        severity="hard_gate",
                        region=region,
                        message=(
                            f"N-glycosylation motif '{motif}' at {chain_name} {region} "
                            f"position {position} — inside a CDR loop, will be glycosylated "
                            "in mammalian expression with unpredictable effect on binding"
                        ),
                    )
                )
            else:
                flags.append(
                    Flag(
                        check="n_glycosylation",
                        severity="soft",
                        region=region,
                        weight=3.0,
                        message=(
                            f"N-glycosylation motif '{motif}' at {chain_name} {region} "
                            f"position {position} — framework, cheap to fix later"
                        ),
                    )
                )

    return flags


# PTM liability motifs. Each entry: (name, regex-free 2-char motif check via
# tuple pairs, base weight). Kept as simple substring/pair scans rather than
# a regex engine since these are all fixed-length adjacent-residue motifs.
DEAMIDATION_PAIRS = {("N", "G"), ("N", "S"), ("N", "T"), ("N", "H")}
ISOMERIZATION_PAIRS = {("D", "G"), ("D", "S")}
OXIDATION_PRONE = {"M", "W"}


def check_ptm_liability(chain: NumberedChain, chain_name: str) -> list[Flag]:
    """Deamidation (N-G/N-S/N-T/N-H), isomerization (D-G/D-S), and
    oxidation-prone residues (M, W). All soft flags, CDR-weighted higher
    than framework. Not implemented anywhere in the ToolUniverse
    antibody-engineering skill (confirmed gap) — built here from scratch.
    """
    flags = []
    seq = chain.full_sequence()
    residues = [r for r in chain.residues if r.aa != "-"]

    def region_weight(region: str, cdr_weight: float, fr_weight: float) -> float:
        return cdr_weight if region.startswith("CDR") else fr_weight

    for i in range(len(seq) - 1):
        pair = (seq[i], seq[i + 1])
        region = residues[i].region
        position = residues[i].position

        if pair in DEAMIDATION_PAIRS:
            flags.append(
                Flag(
                    check="ptm_liability",
                    severity="soft",
                    region=region,
                    weight=region_weight(region, cdr_weight=5.0, fr_weight=1.0),
                    message=f"deamidation motif '{''.join(pair)}' at {chain_name} {region} position {position}",
                )
            )
        elif pair in ISOMERIZATION_PAIRS:
            flags.append(
                Flag(
                    check="ptm_liability",
                    severity="soft",
                    region=region,
                    weight=region_weight(region, cdr_weight=4.0, fr_weight=1.0),
                    message=f"isomerization motif '{''.join(pair)}' at {chain_name} {region} position {position}",
                )
            )

    for i, r in enumerate(residues):
        if r.aa in OXIDATION_PRONE:
            flags.append(
                Flag(
                    check="ptm_liability",
                    severity="soft",
                    region=r.region,
                    weight=region_weight(r.region, cdr_weight=2.0, fr_weight=0.3),
                    message=f"oxidation-prone residue '{r.aa}' at {chain_name} {r.region} position {r.position}",
                )
            )

    return flags


# --- V-domain integrity ---------------------------------------------------
#
# This is the amino-acid-level answer to "is this a complete, translatable
# V domain" — the question IMGT/V-QUEST answers under the name "productivity".
#
# V-QUEST's own definition (in-frame, no premature stop codon) is a property
# of a NUCLEOTIDE sequence and cannot be assessed here: this pipeline takes
# amino acids. Reverse-translating them first, as the V-QUEST tutorial
# workflow does, makes the question tautological — no amino acid maps to a
# stop codon, and every indel is a whole codon, so a reverse-translated
# sequence is in-frame and stop-free by construction. Confirmed empirically:
# zero of PLAbDab's 176,894 heavy chains contain a '*'. (If real
# sequencing-derived nucleotide data ever enters the pipeline, V-QUEST
# productivity becomes a meaningful check again — it is not one for AA input.)
#
# What IS meaningful at the amino acid level, and is what actually breaks in
# practice, is whether the conserved IMGT anchor residues that hold the
# immunoglobulin fold together are present and correct:
#
#   position 23  1st-CYS          forms the intradomain disulfide
#   position 41  CONSERVED-TRP    packs the hydrophobic core
#   position 104 2nd-CYS          the other half of that disulfide
#   position 118 J-PHE or J-TRP   marks the end of a complete V domain
#
# Measured on 400 random PLAbDab pairs: 16 (4%) fail this check — mostly
# truncated before position 118, plus a handful of genuinely substituted
# anchors (C23->S, W41->R, C104 missing). ANARCI numbers all 400 without
# complaint, so it does not catch these on its own.
IMGT_ANCHORS = {
    23: ("C", "1st-CYS"),
    41: ("W", "CONSERVED-TRP"),
    104: ("C", "2nd-CYS"),
    118: ("FW", "J-PHE/J-TRP"),
}

# A missing or substituted anchor is a soft flag, not a hard gate. The
# odd-cysteine episode is the precedent: sequence alone cannot distinguish a
# genuinely broken domain from an unusual-but-real one, and a truncated input
# (someone pasted a partial sequence) is a data-entry problem rather than a
# property of the molecule. Position 118 is weighted lower because
# truncation at the J end is the single most common cause and is usually the
# former.
ANCHOR_MISSING_WEIGHT = 6.0
ANCHOR_WRONG_WEIGHT = 8.0
ANCHOR_118_WEIGHT = 3.0


def check_v_domain_integrity(chain: NumberedChain, chain_name: str) -> list[Flag]:
    """Conserved IMGT anchor residues present and correct (see IMGT_ANCHORS)."""
    flags = []
    observed = {r.position: r.aa for r in chain.residues if r.aa != "-"}

    for position, (expected, name) in IMGT_ANCHORS.items():
        base_weight = ANCHOR_118_WEIGHT if position == 118 else None
        aa = observed.get(position)
        if aa is None:
            flags.append(
                Flag(
                    check="v_domain_integrity",
                    severity="soft",
                    region=_region_for(chain, position),
                    weight=base_weight or ANCHOR_MISSING_WEIGHT,
                    message=(
                        f"conserved {name} at IMGT position {position} is absent "
                        f"({chain_name}) — V domain looks truncated or incomplete"
                    ),
                )
            )
        elif aa not in expected:
            flags.append(
                Flag(
                    check="v_domain_integrity",
                    severity="soft",
                    region=_region_for(chain, position),
                    weight=base_weight or ANCHOR_WRONG_WEIGHT,
                    message=(
                        f"conserved {name} at IMGT position {position} is {aa!r}, "
                        f"expected {'/'.join(expected)} ({chain_name}) — "
                        "the immunoglobulin fold depends on this residue"
                    ),
                )
            )
    return flags


def _region_for(chain: NumberedChain, position: int) -> str:
    for r in chain.residues:
        if r.position == position:
            return r.region
    return "chain"


ALL_CHECKS = [
    check_sequence_sanity,
    check_cysteine_pairing,
    check_n_glycosylation,
    check_ptm_liability,
    check_v_domain_integrity,
]


def run_all_checks(chain: NumberedChain, chain_name: str) -> ChainCheckResult:
    from .developability import check_developability  # local import: avoids a cycle, developability.py imports Flag from here

    result = ChainCheckResult(chain_name=chain_name)
    for check_fn in ALL_CHECKS + [check_developability]:
        result.flags.extend(check_fn(chain, chain_name))
    return result
