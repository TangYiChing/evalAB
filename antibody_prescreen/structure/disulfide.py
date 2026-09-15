"""Disulfide pairing from a predicted structure.

## Why this exists

`checks.check_cysteine_pairing` counts cysteines and flags an odd total. That
is all sequence alone can say, and it is not much: an odd count can be a real
defect or a genuine non-canonical disulfide, and the count says nothing about
*which* cysteine bonds to which. The package README has promised a structural
version since Tier 2 was planned.

Measured consequence of the sequence-only version: it fires on 5.3% of
clinical-stage therapeutics, which is above the 5% false-rejection budget and
therefore disqualifies it as a hard gate. The structural question — "is every
cysteine actually paired in the predicted fold?" — is the one worth asking.

## Implementation note on provenance

The parameter choices below follow the approach in
https://github.com/delvaros/cysteine-bond-prediction (SG-SG distance window,
chi3 dihedral filter, greedy shortest-first pairing). **That repository carries
no licence, so no code was copied from it** — this is an independent
implementation of what is in any case standard structural biology. The
secondary-structure and hydrogen-bond filters in that script are not
reproduced here: they need DSSP, and on a single-conformer predicted model they
would reject on evidence the model does not really carry.

## What the numbers mean

A disulfide bond is an SG-SG distance of ~2.05 A. The window below is
deliberately generous because these are predicted models, not crystal
structures, and ABodyBuilder2's refinement is not geometry-exact.
"""

from dataclasses import dataclass, field
from itertools import combinations
from math import atan2, degrees
from pathlib import Path

SG_DISTANCE_MIN = 1.6
SG_DISTANCE_MAX = 2.8

# chi3 (SG-CB-CB'-SG') clusters near +/-90 degrees in real disulfides. The
# window is wide because a predicted model's side-chain torsions are the least
# reliable part of it; this filter is here to reject coincidental proximity,
# not to grade geometry quality.
CHI3_MIN, CHI3_MAX = -140.0, 140.0

# Minimum separation along the CHAIN so a cysteine cannot pair with its
# immediate neighbour, which is geometrically impossible anyway.
#
# Counted in residues, NOT in IMGT position number. That distinction is not
# pedantic: IMGT numbers CDR-H3 outward from both ends with insertion codes
# (111, 111A...111F, 112F...112A, 112), so two cysteines at IMGT 111D and 112F
# differ by 1 in number while sitting 3 residues apart in the chain. An earlier
# version compared the numbers and silently refused to even evaluate that pair.
# Found on a real design batch where CDR-H3 read C-L-D-C.
MIN_SEQUENCE_GAP = 3

UNPAIRED_WEIGHT = 6.0


@dataclass
class DisulfideResult:
    available: bool
    pairs: list[tuple[str, str]] = field(default_factory=list)  # ("H23", "H104")
    unpaired: list[str] = field(default_factory=list)
    note: str = ""


def _dihedral(p1, p2, p3, p4) -> float | None:
    """Torsion angle p1-p2-p3-p4 in degrees. None if degenerate."""
    def sub(a, b):
        return [a[i] - b[i] for i in range(3)]

    def cross(a, b):
        return [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]

    def dot(a, b):
        return sum(a[i] * b[i] for i in range(3))

    def norm(a):
        length = dot(a, a) ** 0.5
        return [x / length for x in a] if length else None

    b1, b2, b3 = sub(p2, p1), sub(p3, p2), sub(p4, p3)
    n1, n2, b2n = norm(cross(b1, b2)), norm(cross(b2, b3)), norm(b2)
    if n1 is None or n2 is None or b2n is None:
        return None
    m1 = cross(n1, b2n)
    return degrees(atan2(dot(m1, n2), dot(n1, n2)))


def find_disulfides(model_path: Path | str) -> DisulfideResult:
    """Pair the cysteines in an ABodyBuilder2 model.

    Pairs are found across chains as well as within them — the inter-chain
    disulfide is real antibody biology and a within-chain-only search would
    report both partners as unpaired.
    """
    try:
        from Bio import PDB
    except ImportError:
        return DisulfideResult(available=False, note="biopython not installed")

    try:
        structure = PDB.PDBParser(QUIET=True).get_structure("model", str(model_path))
    except Exception as e:
        return DisulfideResult(available=False, note=f"could not parse model: {e}")

    cysteines = []
    for chain in structure[0]:
        # Ordinal position along the chain, which is what "three residues
        # apart" actually means. residue.id[1] is the IMGT number and is not
        # monotonic across insertion codes in CDR-H3.
        for index, residue in enumerate(chain):
            if residue.get_resname() != "CYS" or "SG" not in residue:
                continue
            label = f"{chain.id}{residue.id[1]}{residue.id[2].strip()}"
            cysteines.append((label, chain.id, index, residue))

    if not cysteines:
        return DisulfideResult(available=True, pairs=[], unpaired=[])

    candidates = []
    for (l1, c1, i1, r1), (l2, c2, i2, r2) in combinations(cysteines, 2):
        if c1 == c2 and abs(i1 - i2) < MIN_SEQUENCE_GAP:
            continue
        distance = r1["SG"] - r2["SG"]
        if not (SG_DISTANCE_MIN <= distance <= SG_DISTANCE_MAX):
            continue
        if "CB" in r1 and "CB" in r2:
            chi3 = _dihedral(
                r1["CB"].get_coord(), r1["SG"].get_coord(),
                r2["SG"].get_coord(), r2["CB"].get_coord(),
            )
            if chi3 is not None and not (CHI3_MIN <= chi3 <= CHI3_MAX):
                continue
        candidates.append((distance, l1, l2))

    # Greedy shortest-first: a sulfur can only take part in one bond, so the
    # closest unambiguous pair wins and both partners are then unavailable.
    candidates.sort()
    taken: set[str] = set()
    pairs = []
    for _distance, l1, l2 in candidates:
        if l1 in taken or l2 in taken:
            continue
        pairs.append((l1, l2))
        taken.update((l1, l2))

    unpaired = [label for label, *_ in cysteines if label not in taken]
    return DisulfideResult(available=True, pairs=pairs, unpaired=unpaired)


def disulfide_flags(result: DisulfideResult) -> list:
    """Unpaired cysteines as weighted soft flags.

    Soft, not a hard gate: a free cysteine is a real manufacturability concern
    (it drives covalent aggregation and disulfide scrambling) but it is also
    real biology in some approved molecules, and this is a *predicted* model.
    Same reasoning that downgraded the sequence-level cysteine gate.
    """
    from ..checks import Flag

    if not result.available or not result.unpaired:
        return []
    return [
        Flag(
            check="disulfide_pairing",
            severity="soft",
            region="structure",
            weight=UNPAIRED_WEIGHT,
            message=(
                f"{len(result.unpaired)} unpaired cysteine(s) in the predicted "
                f"structure: {', '.join(result.unpaired)} "
                f"({len(result.pairs)} disulfide(s) found)"
            ),
        )
    ]
