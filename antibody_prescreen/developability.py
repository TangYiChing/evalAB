"""Aggregation propensity (AGGRESCAN) and net charge.

Renamed from "developability", which was too vague to be useful: it bundled
three unrelated measurements under one name, and a flag saying
"developability" told the reader nothing about what was actually wrong.
Aggregation and charge are separable, mechanistically different, and each is
something a person can picture.

They are also emitted as separate checks now, because bundling them made both
uninterpretable: the combined check fired on 100% of candidates in both
populations (LR exactly 1.00), and there was no way to tell whether that came
from the aggregation part, the charge part, or the unconditional flag.

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


def aggrescan_na4vss(seq: str) -> float:
    """Aggregate aggregation score (Na4vSS): mean of (a4v - threshold) clipped
    at 0, averaged over the whole sequence. This is the metric AGGRESCAN
    itself uses for whole-protein comparison — almost every real protein has
    a handful of individual hot spots (that's normal, not a defect), so
    per-hot-spot counting is the wrong signal. The aggregate magnitude,
    compared against other candidates in the same batch, is the right one.
    """
    a3v = [A3V.get(aa, 0.0) for aa in seq]
    a4v = _window(a3v, 5)
    return round(sum(max(0.0, v - HOT_SPOT_THRESHOLD) for v in a4v) / len(a4v), 4)


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


def check_aggregation_and_charge(chain: NumberedChain, chain_name: str) -> list[Flag]:
    seq = chain.full_sequence()
    residues = [r for r in chain.residues if r.aa != "-"]

    flags: list[Flag] = []

    hot_spots = aggrescan_hot_spots(seq)
    na4vss = aggrescan_na4vss(seq)
    cdr_overlapping_spots = []
    for spot in hot_spots:
        spot_region_residues = residues[spot["start_idx"] : spot["end_idx"] + 1]
        if any(r.is_cdr for r in spot_region_residues):
            region_label = "/".join(sorted({r.region for r in spot_region_residues}))
            cdr_overlapping_spots.append((spot, region_label))

    # Aggregate signal, not per-hot-spot: almost every real antibody has a
    # handful of individual AGGRESCAN hot spots (normal biology), so the
    # whole-sequence Na4vSS magnitude is the meaningful comparison, weighted
    # modestly since this is one signal among several.
    flags.append(
        Flag(
            check="aggregation",
            severity="soft",
            region="chain",
            weight=na4vss * 15.0,
            message=(
                f"{chain_name} aggregate aggregation score (Na4vSS)={na4vss}, "
                f"{len(hot_spots)} hot spot(s) total, {len(cdr_overlapping_spots)} overlapping a CDR"
            ),
        )
    )

    # Separately and more lightly: flag when a hot spot specifically sits in
    # a CDR loop, since that's the localized-risk case worth a human's
    # attention even when the aggregate score is unremarkable.
    for spot, region_label in cdr_overlapping_spots:
        flags.append(
            Flag(
                check="aggregation",
                severity="soft",
                region=region_label,
                weight=2.0,
                message=(
                    f"aggregation hot spot '{spot['peptide']}' overlaps CDR "
                    f"({chain_name} {region_label}, mean a4v={spot['mean_a4v']})"
                ),
            )
        )

    pi = isoelectric_point(seq)
    if pi < PI_LOW or pi > PI_HIGH:
        flags.append(
            Flag(
                check="charge",
                severity="soft",
                region="chain",
                weight=2.0,
                message=f"{chain_name} pI={pi} outside typical Fv range [{PI_LOW}, {PI_HIGH}] — may affect solubility/formulation",
            )
        )

    return flags
