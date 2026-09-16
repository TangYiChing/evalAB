"""Likelihood-ratio analysis: decide where each check belongs.

Replaces hand-set weights with an estimated quantity. For a binary check:

    LR+ = P(flag | bad) / P(flag | good)     how much a flag should worry you
    LR- = P(no flag | bad) / P(no flag | good)   how much a pass should reassure you

The thing that decides the triage LEVEL is not the LR. It is
P(flag | good) — the false-rejection rate — because rejecting a candidate
(Level 5) is irreversible and what makes it safe is not firing on good ones. The
LR tells you how much the flag is worth once it fires. These are different
axes and conflating them is a mistake this module exists partly to prevent.

For a CONTINUOUS score (TAP PSH, aggregation, model confidence), one LR for
the whole range hides where the information is. `interval_lr` splits the range
into bins and reports an LR per bin. Measured on TAP's hydrophobic patch
score, that turns a whole-range AUC of 0.519 — chance — into a visible LR of
7.24 in the top bin. The information was always in the tail.

Usage:
    python -m antibody_prescreen.calibration.lr_report results.csv
    python -m antibody_prescreen.calibration.lr_report results.csv --interval tap_hydrophobic
"""

import argparse
import csv
import math
from pathlib import Path

POSITIVE_LABEL = "TheraSAbDab"
NEGATIVE_LABEL = "Patent text"

# Trauma-triage benchmark: under-triage below 5%. A check may only become a
# rejecting gate (Level 5) if the UPPER confidence bound of its false-rejection
# rate on known-good candidates clears this.
FALSE_REJECTION_BUDGET = 0.05


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval — behaves sensibly at k = 0, which normal
    approximation does not, and every interesting gate here has k = 0."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def likelihood_ratios(a: int, n_good: int, b: int, n_bad: int) -> tuple[float, float]:
    """LR+ and LR- with a Haldane-Anscombe 0.5 correction.

    The correction matters: our best gates have zero flags on known-good, and
    an uncorrected ratio would be infinite, which is not a usable number and
    hides that the estimate rests on a small sample.
    """
    pg = (a + 0.5) / (n_good + 1)
    pb = (b + 0.5) / (n_bad + 1)
    return pb / pg, (1 - pb) / (1 - pg)


def recommend_level(flag_rate_good: float, upper_bound: float, lr_plus: float) -> str:
    """Where this check belongs, and why."""
    if upper_bound < FALSE_REJECTION_BUDGET:
        return "REJECTING GATE (routes to Level 5)"
    if flag_rate_good > 0.90:
        return "repair count only — fires on nearly everything"
    if lr_plus > 1.3:
        return "REVIEW (shown, never routes)"
    return "repair count only — carries no discrimination"


def binary_report(rows: list[dict], check_columns: list[str]) -> None:
    good = [r for r in rows if r["pairing"] == POSITIVE_LABEL]
    bad = [r for r in rows if r["pairing"] == NEGATIVE_LABEL]
    print(f"Known-good n={len(good)}   background n={len(bad)}")
    print(f"False-rejection budget: {100 * FALSE_REJECTION_BUDGET:.0f}% "
          f"(upper 95% bound must clear it)\n")

    header = (
        f"{'check':26s} {'P(+|good)':>10s} {'95% upper':>10s} "
        f"{'P(+|bad)':>9s} {'LR+':>7s} {'LR-':>7s}  recommendation"
    )
    print(header)
    print("-" * len(header))

    results = []
    for column in check_columns:
        a = sum(1 for r in good if _truthy(r.get(column)))
        b = sum(1 for r in bad if _truthy(r.get(column)))
        if a == 0 and b == 0 and not any(column in r for r in rows):
            continue
        lo, hi = wilson(a, len(good))
        lr_plus, lr_minus = likelihood_ratios(a, len(good), b, len(bad))
        rate = a / len(good) if good else 0.0
        results.append((column, rate, hi, b / len(bad) if bad else 0.0, lr_plus, lr_minus))

    for column, rate, hi, rate_bad, lr_plus, lr_minus in sorted(
        results, key=lambda x: -x[4]
    ):
        print(
            f"{column:26s} {100 * rate:9.1f}% {100 * hi:9.1f}% {100 * rate_bad:8.1f}% "
            f"{lr_plus:7.2f} {lr_minus:7.3f}  "
            f"{recommend_level(rate, hi, lr_plus)}"
        )

    print(
        "\nNote: every LR here is DERIVATION, measured on the same population the "
        "\nchecks were inspected on. None is externally validated. See "
        "\ndocs/triage_implementation_plan.md Phase 2."
    )


def _truthy(value) -> bool:
    if value is None or value == "":
        return False
    try:
        return float(value) > 0
    except (TypeError, ValueError):
        return str(value).strip().lower() in {"true", "yes", "1", "red", "amber"}


def interval_report(rows: list[dict], column: str, bins: int = 6) -> None:
    """Per-bin LR across a continuous score's range."""
    def value(r):
        try:
            return float(r[column])
        except (KeyError, TypeError, ValueError):
            return None

    good = [v for r in rows if r["pairing"] == POSITIVE_LABEL and (v := value(r)) is not None]
    bad = [v for r in rows if r["pairing"] == NEGATIVE_LABEL and (v := value(r)) is not None]
    if not good or not bad:
        print(f"{column}: not enough data")
        return

    everything = sorted(good + bad)
    edges = [everything[0]] + [
        everything[int(i * len(everything) / bins)] for i in range(1, bins)
    ] + [everything[-1] + 1e-9]

    print(f"\nInterval LR for '{column}'   (good n={len(good)}, bad n={len(bad)})")
    print(f"{'range':>22s} {'good':>6s} {'bad':>6s} {'P(|good)':>9s} {'P(|bad)':>8s} {'LR':>7s}")
    print("-" * 64)
    for lo, hi in zip(edges, edges[1:]):
        a = sum(1 for v in good if lo <= v < hi)
        b = sum(1 for v in bad if lo <= v < hi)
        pg = (a + 0.5) / (len(good) + 1)
        pb = (b + 0.5) / (len(bad) + 1)
        print(
            f"{lo:10.2f}–{hi:<10.2f} {a:6d} {b:6d} "
            f"{100 * a / len(good):8.1f}% {100 * b / len(bad):7.1f}% {pb / pg:7.2f}"
        )
    print(
        "\nA flat column of ~1.00 means this score carries no information as a RANK. "
        "\nAn LR that climbs only in the extreme bin means it works as a BOUNDARY, not "
        "\na ranking — use it as a gate on the tail, not as a score."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--interval", action="append", default=[],
                        help="continuous column to report interval LRs for")
    parser.add_argument("--bins", type=int, default=6)
    args = parser.parse_args()

    with open(args.csv_path) as handle:
        rows = [r for r in csv.DictReader(handle) if not r.get("error")]

    check_columns = [c for c in rows[0] if c.startswith("flag_")]
    if check_columns:
        binary_report(rows, check_columns)
    else:
        print("No flag_* columns found — re-run run_calibration.py to emit them.")

    for column in args.interval:
        interval_report(rows, column, args.bins)


if __name__ == "__main__":
    main()
