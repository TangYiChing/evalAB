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


def check_cysteine_pairing(chain: NumberedChain, chain_name: str) -> list[Flag]:
    """Flag odd cysteine counts. Real pairing (which Cys bonds to which) needs
    a structure — this is the sequence-only proxy: an odd count means at
    least one cysteine can't be disulfide-paired at all.
    """
    flags = []
    seq = chain.full_sequence()
    cys_count = seq.count("C")

    if cys_count % 2 != 0:
        flags.append(
            Flag(
                check="cysteine_pairing",
                severity="hard_gate",
                region="chain",
                message=(
                    f"{chain_name} has an odd cysteine count ({cys_count}) — at least one "
                    "cysteine cannot be disulfide-paired, high risk of misfolding/aggregation"
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


ALL_CHECKS = [
    check_sequence_sanity,
    check_cysteine_pairing,
    check_n_glycosylation,
    check_ptm_liability,
]


def run_all_checks(chain: NumberedChain, chain_name: str) -> ChainCheckResult:
    result = ChainCheckResult(chain_name=chain_name)
    for check_fn in ALL_CHECKS:
        result.flags.extend(check_fn(chain, chain_name))
    return result
