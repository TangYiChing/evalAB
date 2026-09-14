"""Sequence-based developability checks: aggregation propensity (AGGRESCAN),
isoelectric point, and hydrophobic patches.

Core algorithm ported from the ToolUniverse `tooluniverse-antibody-engineering`
skill's `scripts/developability.py` (pure-Python AGGRESCAN implementation,
Conchillo-Sole et al. 2007 values) rather than reimplemented, per the Tier 1
plan's "reuse, don't rebuild" note. This module adds the Flag-based wrapper
and CDR-region weighting that the original script didn't have (it returns raw
JSON with no region awareness).
"""

from .checks import Flag
from .numbering import NumberedChain

# --- AGGRESCAN intrinsic aggregation-propensity values (a3v), Conchillo-Sole 2007 ---
A3V = {
    "A": -0.036, "R": -1.240, "N": -1.302, "D": -1.836, "C": 0.604,
    "Q": -1.231, "E": -1.412, "G": -0.535, "H": -1.033, "I": 1.822,
    "L": 1.380, "K": -0.931, "M": 0.910, "F": 1.754, "P": -1.236,
    "S": -0.294, "T": -0.159, "W": 1.037, "Y": 1.159, "V": 1.594,
}
HOT_SPOT_THRESHOLD = -0.02  # AGGRESCAN HST

PKA = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1, "H": 6.5, "K": 10.8, "R": 12.5}
N_TERM, C_TERM = 8.6, 3.6

# pI range considered "normal" for a human IgG Fv; outside this is a soft flag.
PI_LOW, PI_HIGH = 5.5, 9.5


def _window(values, w=5):
    half = w // 2
    out = []
    n = len(values)
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        seg = values[lo:hi]
        out.append(sum(seg) / len(seg))
    return out


def aggrescan_hot_spots(seq: str) -> list[dict]:
    a3v = [A3V.get(aa, 0.0) for aa in seq]
    a4v = _window(a3v, 5)
    spots = []
    start = None
    for i, v in enumerate(a4v):
        if v > HOT_SPOT_THRESHOLD and start is None:
            start = i
        elif v <= HOT_SPOT_THRESHOLD and start is not None:
            if i - start >= 5:
                spots.append((start, i - 1))
            start = None
    if start is not None and len(a4v) - start >= 5:
        spots.append((start, len(a4v) - 1))
    return [
        {
            "start_idx": s,  # 0-based index into seq, inclusive
            "end_idx": e,  # 0-based index into seq, inclusive
            "peptide": seq[s : e + 1],
            "mean_a4v": round(sum(a4v[s : e + 1]) / (e - s + 1), 3),
        }
        for s, e in spots
    ]


def isoelectric_point(seq: str) -> float:
    counts = {aa: seq.count(aa) for aa in set(seq)}

    def charge(ph):
        pos = 10 ** (N_TERM - ph) / (1 + 10 ** (N_TERM - ph))
        pos += counts.get("K", 0) * 10 ** (PKA["K"] - ph) / (1 + 10 ** (PKA["K"] - ph))
        pos += counts.get("R", 0) * 10 ** (PKA["R"] - ph) / (1 + 10 ** (PKA["R"] - ph))
        pos += counts.get("H", 0) * 10 ** (PKA["H"] - ph) / (1 + 10 ** (PKA["H"] - ph))
        neg = 1 / (1 + 10 ** (C_TERM - ph))
        neg += counts.get("D", 0) / (1 + 10 ** (PKA["D"] - ph))
        neg += counts.get("E", 0) / (1 + 10 ** (PKA["E"] - ph))
        neg += counts.get("C", 0) / (1 + 10 ** (PKA["C"] - ph))
        neg += counts.get("Y", 0) / (1 + 10 ** (PKA["Y"] - ph))
        return pos - neg

    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def check_developability(chain: NumberedChain, chain_name: str) -> list[Flag]:
    seq = chain.full_sequence()
    residues = [r for r in chain.residues if r.aa != "-"]

    flags: list[Flag] = []

    for spot in aggrescan_hot_spots(seq):
        spot_region_residues = residues[spot["start_idx"] : spot["end_idx"] + 1]
        overlaps_cdr = any(r.is_cdr for r in spot_region_residues)
        region_label = "/".join(sorted({r.region for r in spot_region_residues}))
        flags.append(
            Flag(
                check="developability",
                severity="soft",
                region=region_label,
                weight=6.0 if overlaps_cdr else 2.0,
                message=(
                    f"aggregation hot spot '{spot['peptide']}' ({chain_name} {region_label}, "
                    f"mean a4v={spot['mean_a4v']}){' — overlaps CDR' if overlaps_cdr else ''}"
                ),
            )
        )

    pi = isoelectric_point(seq)
    if pi < PI_LOW or pi > PI_HIGH:
        flags.append(
            Flag(
                check="developability",
                severity="soft",
                region="chain",
                weight=2.0,
                message=f"{chain_name} pI={pi} outside typical Fv range [{PI_LOW}, {PI_HIGH}] — may affect solubility/formulation",
            )
        )

    return flags
