"""Step 7d/7e: read a calibration CSV and report what it says about weights.

Deliberately does NOT auto-write new constants. It prints the evidence and a
recommendation; a human decides what goes into fusion.py, because the label
here is provenance ("reached the clinic" vs "was patented"), not an assay, and
a number fitted to it should be adopted knowingly.

Usage:
    python -m antibody_prescreen.calibration.fit_thresholds calibration_results.csv
"""

import argparse
import csv
import statistics
from pathlib import Path

THERAPEUTIC = "TheraSAbDab"
BACKGROUND = "Patent text"


def _floats(rows, key):
    out = []
    for row in rows:
        try:
            out.append(float(row[key]))
        except (ValueError, KeyError, TypeError):
            pass
    return out


def describe(name, values):
    if not values:
        return f"  {name:24s} (no data)"
    values = sorted(values)
    def pct(p):
        return values[min(len(values) - 1, int(p * len(values)))]
    return (
        f"  {name:24s} n={len(values):4d}  median={statistics.median(values):7.2f}  "
        f"p10={pct(0.10):7.2f}  p90={pct(0.90):7.2f}  max={max(values):7.2f}"
    )


def separation(therapeutic, background):
    """How much of the background sits above the therapeutic median.

    Not an AUC — just the number the threshold decision actually turns on: if
    a cutoff is placed at the therapeutic median, what fraction of each
    population lands above it.
    """
    if not therapeutic or not background:
        return None
    cutoff = statistics.median(therapeutic)
    above_t = sum(1 for v in therapeutic if v > cutoff) / len(therapeutic)
    above_b = sum(1 for v in background if v > cutoff) / len(background)
    return cutoff, above_t, above_b


def auc(therapeutic, background):
    """P(a random background scores higher than a random therapeutic).

    0.5 = the component carries no information about the label. Computed by
    direct pair counting — the populations here are small enough that this is
    cheaper than pulling in scipy.
    """
    if not therapeutic or not background:
        return None
    wins = ties = 0
    for b in background:
        for t in therapeutic:
            if b > t:
                wins += 1
            elif b == t:
                ties += 1
    total = len(background) * len(therapeutic)
    return (wins + 0.5 * ties) / total if total else None


def interpret_auc(value: float, n_t: int, n_b: int) -> str:
    """Plain-language reading of an AUC, with an honest noise floor.

    AUC is a ranking statistic on small samples, so a value near 0.5 says
    nothing and a value slightly above it says almost nothing. The rough
    standard error of AUC under the null is ~sqrt(1/(12*n)) with n the smaller
    group; anything inside ~2 SE of 0.5 is reported as noise rather than
    dressed up as a weak signal.
    """
    n = max(1, min(n_t, n_b))
    noise = 2 * (1.0 / (12 * n)) ** 0.5
    delta = value - 0.5
    if abs(delta) < noise:
        return f"indistinguishable from chance (+/-{noise:.3f} noise floor at n={n})"
    if delta < 0:
        return "INVERTED — background scores LOWER than therapeutics"
    if delta < 0.10:
        return "weak separation"
    if delta < 0.20:
        return "moderate separation"
    return "strong separation"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", type=Path)
    args = parser.parse_args()

    with open(args.csv_path) as handle:
        rows = [r for r in csv.DictReader(handle) if not r["error"]]

    therapeutic = [r for r in rows if r["pairing"] == THERAPEUTIC]
    background = [r for r in rows if r["pairing"] == BACKGROUND]
    print(f"Scored: {len(therapeutic)} therapeutic, {len(background)} background\n")

    components = [
        ("tier1_score", "Tier 1"),
        ("tap_weight", "TAP weight"),
        ("immuno_weight", "Immunogenicity weight"),
        ("total_score", "TOTAL"),
    ]

    print("=== Component distributions ===")
    for key, label in components:
        t, b = _floats(therapeutic, key), _floats(background, key)
        print(f"\n{label}:")
        print(describe("therapeutic", t))
        print(describe("background", b))
        a = auc(t, b)
        if a is not None:
            print(f"  AUC(background > therapeutic) = {a:.3f}  -> {interpret_auc(a, len(t), len(b))}")

    print("\n=== Where would thresholds land? ===")
    total_t, total_b = _floats(therapeutic, "total_score"), _floats(background, "total_score")
    if total_t and total_b:
        sorted_t = sorted(total_t)
        t_low = statistics.median(sorted_t)
        t_high = sorted_t[min(len(sorted_t) - 1, int(0.90 * len(sorted_t)))]
        print(f"  T_LOW  at therapeutic median = {t_low:.2f}")
        print(f"  T_HIGH at therapeutic p90    = {t_high:.2f}")
        for label, values in (("therapeutic", total_t), ("background", total_b)):
            go = sum(1 for v in values if v < t_low)
            cond = sum(1 for v in values if t_low <= v < t_high)
            nogo = sum(1 for v in values if v >= t_high)
            n = len(values)
            print(
                f"    {label:12s} GO {go:3d} ({100*go/n:4.1f}%)  "
                f"CONDITIONAL {cond:3d} ({100*cond/n:4.1f}%)  "
                f"NO-GO {nogo:3d} ({100*nogo/n:4.1f}%)"
            )

    print("\n=== Per-TAP-metric red/amber rates ===")
    for key in ("total", "hydrophobic", "positive", "negative", "sfvcsp"):
        col = "tap_flag_" + key
        line = f"  {key:12s}"
        for label, group in (("thera", therapeutic), ("backg", background)):
            flags = [r[col] for r in group if r.get(col)]
            if not flags:
                continue
            red = sum(1 for f in flags if f == "RED") / len(flags)
            amber = sum(1 for f in flags if f == "AMBER") / len(flags)
            line += f"  {label}: {100*red:4.1f}% red / {100*amber:4.1f}% amber"
        print(line)

    print("\n=== Immunogenicity flag counts ===")
    for key in ("n_immuno_cdr", "n_immuno_framework", "n_immuno_observed",
                "n_predicted_binders", "n_dropped_germline"):
        t, b = _floats(therapeutic, key), _floats(background, key)
        if not t and not b:
            continue
        mt = statistics.mean(t) if t else float("nan")
        mb = statistics.mean(b) if b else float("nan")
        print(f"  {key:22s} therapeutic mean={mt:5.2f}   background mean={mb:5.2f}")


if __name__ == "__main__":
    main()
