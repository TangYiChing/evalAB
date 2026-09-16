"""Screen a batch of designed VH/VL pairs, with an optional sealed outcome label.

    python examples/screen_designs.py data/example/demo_batch.csv
    python examples/screen_designs.py my_batch.csv --shared-scaffold
    python examples/screen_designs.py my_batch.csv --unseal qc_outcome

This is a thin wrapper over the public API — `screen_batch`, `provenance` and
`format_report` — and deliberately adds no decision logic of its own. If you
want the same output without the extras, `python -m antibody_prescreen` does
the whole job.

What the extras are for:

`--shared-scaffold` declares that the batch is a congeneric series: one parent,
one framework, variation in a loop or two. Triage takes the MAX over findings,
so any liability the parent already carried appears in every candidate and
saturates the level for the entire batch. The level stays a true statement
about each candidate and becomes a useless one for choosing between them. The
introduced-only view sets the shared findings aside so the variants can be
compared. It is a declaration and never an inference, because declaring it for
an unrelated batch produces a view with no meaning.

`--unseal` reveals a held-out outcome column AFTER triage. Nothing between the
CSV reader and the triage ever reads that column, which is the point: it is
how you check whether the levels track your own wet-lab endpoint without
letting the endpoint influence the levels. See docs/calibration.md §2 for the
pre-registration rule that makes such a comparison worth anything.
"""

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from antibody_prescreen import format_report, provenance, screen_batch
from antibody_prescreen.bands import BandSet
from antibody_prescreen.screen import load_default_bands


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv_path", type=Path)
    ap.add_argument("--bands", type=Path, default=None,
                    help="band file; default is the repertoire + clinical set shipped with evalAB")
    ap.add_argument("--id-column", default="candidate_id")
    ap.add_argument("--vh-column", default="vh_sequence")
    ap.add_argument("--vl-column", default="vl_sequence")
    ap.add_argument("--shared-scaffold", action="store_true",
                    help="declare a congeneric batch and add the introduced-only view")
    ap.add_argument("--unseal", default=None,
                    help="column holding a held-out outcome label, revealed AFTER triage")
    ap.add_argument("--out", type=Path, default=None, help="write a per-candidate CSV")
    args = ap.parse_args()

    rows = list(csv.DictReader(args.csv_path.open(newline="")))
    if not rows:
        raise SystemExit(f"{args.csv_path} has no rows")

    candidates = [
        {
            "candidate_id": row.get(args.id_column) or f"design-{i + 1:02d}",
            "vh_sequence": row[args.vh_column],
            "vl_sequence": row[args.vl_column],
        }
        for i, row in enumerate(rows)
    ]

    band_set = BandSet.load(args.bands) if args.bands else load_default_bands()
    results = screen_batch(candidates, band_set=band_set)
    run = provenance(results, band_set)

    print(format_report(results, shared_scaffold=args.shared_scaffold, run_provenance=run))

    dist: dict[int, int] = {}
    for result in results:
        dist[result.triage.level] = dist.get(result.triage.level, 0) + 1
    print("\n### Level distribution\n")
    for level in sorted(dist):
        print(f"  L{level}  {dist[level]:3d}  {100 * dist[level] / len(results):5.1f}%")

    if args.out:
        _write_csv(args.out, results)
        print(f"\nwrote {args.out}")

    if args.unseal:
        _unseal(args.unseal, rows, results)


def _write_csv(path: Path, results) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["candidate_id", "level", "level_name", "action",
                         "n_findings", "drivers", "unmeasured", "error"])
        for result in results:
            triage = result.triage
            writer.writerow([
                result.candidate_id, triage.level, triage.name, triage.action,
                len(triage.findings),
                " | ".join(f.message for f in triage.drivers),
                " | ".join(triage.unmeasured),
                triage.error or "",
            ])


def _unseal(column: str, rows: list[dict], results) -> None:
    """Reveal the held-out label and show the level distribution per label.

    This is a description, not a validation. It becomes evidence only if the
    endpoint and its thresholds were written down before the batch ran, and
    only on a holdout that has not been looked at before.
    """
    print(f"\n### Held-out label: {column}\n")
    by_label: dict[str, list[int]] = {}
    for row, result in zip(rows, results):
        by_label.setdefault(row.get(column, "<missing>"), []).append(result.triage.level)
    for label, levels in sorted(by_label.items()):
        levels.sort()
        mean = sum(levels) / len(levels)
        print(f"  {label:20s} n={len(levels):3d}  levels={levels}  mean={mean:.2f}")
    print("\n  Reminder: this is only evidence if the endpoint, its thresholds and\n"
          "  this holdout were fixed before the run — see docs/calibration.md §2.")


if __name__ == "__main__":
    main()
