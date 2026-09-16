"""Sequence-level aggregation propensity and charge.

Three separable measurements, each reported on its own: AGGRESCAN aggregation
propensity, isoelectric point, and net charge. They are kept apart on purpose.
Bundling them under one "developability" heading tells a reader that something
is off without telling them which physical property it was, and the three have
different repair costs — a charge problem in a CDR is a different kind of
problem from a hydrophobic stretch in a framework.

Each function here returns a raw number. Nothing in this module decides whether
that number is acceptable: the value is located against a stated reference band
in `bands.py`, and the deviation is paired with a repair cost in `levels.py`.
That separation is what keeps the pipeline a triage rather than a score.

The AGGRESCAN core (Conchillo-Sole et al. 2007 a3v values, windowed a4v, and
the Na4vSS aggregate) is a port of the pure-Python implementation in the
ToolUniverse `tooluniverse-antibody-engineering` skill's
`scripts/developability.py`, rather than a reimplementation. The CDR-region
awareness is added here; the original is region-blind.

`net_charge` and `isoelectric_point` deliberately share one Henderson-
Hasselbalch model and one pKa table, so the two can never disagree with each
other for the sole reason that they were computed from different constants.
"""

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


def net_charge(seq: str, ph: float = 7.4) -> float:
    """Net charge at a given pH, on the same Henderson-Hasselbalch model as
    `isoelectric_point` — which is the point: pI and net charge must not
    disagree with each other because they were computed from different pKa
    tables.

    Reported for the whole Fv and, separately, for the CDR residues alone.
    The second one is the sequence-level stand-in for TAP's PPC/PNC patch
    metrics: the patches themselves need a structure, but the charge that
    forms them is already visible in the loops.
    """
    counts = {aa: seq.count(aa) for aa in set(seq)}
    pos = 10 ** (N_TERM - ph) / (1 + 10 ** (N_TERM - ph))
    for aa in ("K", "R", "H"):
        pos += counts.get(aa, 0) * 10 ** (PKA[aa] - ph) / (1 + 10 ** (PKA[aa] - ph))
    neg = 1 / (1 + 10 ** (C_TERM - ph))
    for aa in ("D", "E", "C", "Y"):
        neg += counts.get(aa, 0) / (1 + 10 ** (PKA[aa] - ph))
    return round(pos - neg, 3)
