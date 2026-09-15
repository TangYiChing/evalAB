"""Step 7: score a labelled population and dump per-component results.

Produces one CSV row per candidate with every component of the score kept
separate — Tier 1, each TAP metric, the immunogenicity layers — so thresholds
and weights can be fitted afterwards without re-running the expensive part.

Usage:
    python -m antibody_prescreen.calibration.run_calibration \
        --n-per-group 150 --workers 5 --out calibration_results.csv

Restricted to human-framework antibodies by default. Not because species
confounds the two populations (measured: 74.8% both-human in TheraSAbDab vs
70.4% in patent text — a 4.4 point difference), but because a non-human
framework is correctly scored as almost entirely non-self against human
germline and lands a large, systematic immunogenicity offset. Leaving those
in would let "is it murine" do part of the discrimination, which is a real
signal but a different one from developability. Use --include-non-human to
keep them.
"""

import argparse
import csv
import sys
from concurrent.futures import BrokenExecutor, ProcessPoolExecutor, as_completed
from pathlib import Path

from ..fusion import screen_candidate
from ..numbering import NumberingError, number_antibody
from .plabdab_source import (
    BACKGROUND_PAIRING,
    THERAPEUTIC_PAIRING,
    load_paired_sequences,
)

TAP_METRICS = [
    "Total IMGT CDR Length",
    "Hydrophobic Patch Score",
    "Positive Patch Score",
    "Negative Patch Score",
    "SFvCSP",
]

# Per-check boolean columns so likelihood ratios can be recomputed without
# re-scoring the population. lr_report.py reads every flag_* column.
CHECK_NAMES = [
    "sequence_sanity", "cysteine_pairing", "n_glycosylation", "ptm_liability",
    "v_domain_integrity", "poly_residue_run", "developability", "humanness",
    "tap", "disulfide_pairing", "immunogenicity", "immunogenicity_observed",
]

FIELDNAMES = (
    ["candidate_id", "pairing", "triage_level", "cdr_repairs",
     "framework_repairs", "total_score", "tier1_score", "error"]
    + [f"flag_{c}" for c in CHECK_NAMES]
    + ["flag_tap_any_red", "model_conf_max", "model_conf_mean", "vh_germline_identity"]
    + ["tap_" + m.split()[0].lower() for m in TAP_METRICS]
    + ["tap_flag_" + m.split()[0].lower() for m in TAP_METRICS]
    + [
        "tap_weight",
        "n_immuno_cdr",
        "n_immuno_framework",
        "n_immuno_observed",
        "immuno_weight",
        "n_predicted_binders",
        "n_dropped_germline",
        "vh_v_gene",
        "vl_v_gene",
        "structure_available",
        "immunogenicity_available",
    ]
)


def is_human_framework(vh: str, vl: str) -> bool:
    try:
        chains = number_antibody(vh, vl)
    except NumberingError:
        return False
    return chains["VH"].v_species == "human" and chains["VL"].v_species == "human"


def score_one(candidate: dict, model_cache: str) -> dict:
    """Score one candidate and flatten every component into a flat row."""
    row = {name: "" for name in FIELDNAMES}
    row["candidate_id"] = candidate["candidate_id"]
    row["pairing"] = candidate["pairing"]

    try:
        result = screen_candidate(
            candidate["candidate_id"],
            candidate["vh_sequence"],
            candidate["vl_sequence"],
            run_immunogenicity=True,
            run_structure=True,
            model_cache=model_cache,
        )
    except Exception as e:  # a calibration run must not die on one bad candidate
        row["error"] = f"{type(e).__name__}: {e}"
        return row

    if result.triage_result is not None:
        row["triage_level"] = result.triage_result.level
        row["cdr_repairs"] = result.triage_result.cdr_repairs
        row["framework_repairs"] = result.triage_result.framework_repairs
    fired = {f.check for f in result.all_flags}
    for name in CHECK_NAMES:
        row[f"flag_{name}"] = int(name in fired)
    row["flag_tap_any_red"] = int(
        bool(result.tap_profile) and result.tap_profile.n_red > 0
    )
    if result.model_confidence:
        row["model_conf_max"] = round(result.model_confidence.get("max") or 0, 4)
        row["model_conf_mean"] = round(result.model_confidence.get("mean") or 0, 4)
    row["total_score"] = "" if result.score == float("inf") else round(result.score, 3)
    row["structure_available"] = int(result.structure_available)
    row["immunogenicity_available"] = int(result.immunogenicity_available)
    if result.error:
        row["error"] = result.error

    tap_weight = sum(f.weight for f in result.all_flags if f.check == "tap")
    immuno_flags = [f for f in result.all_flags if f.check == "immunogenicity"]
    observed_flags = [f for f in result.all_flags if f.check == "immunogenicity_observed"]
    immuno_weight = sum(f.weight for f in immuno_flags + observed_flags)

    row["tap_weight"] = round(tap_weight, 3)
    row["immuno_weight"] = round(immuno_weight, 3)
    row["n_immuno_cdr"] = sum(1 for f in immuno_flags if "overlaps CDR" in f.message)
    row["n_immuno_framework"] = len(immuno_flags) - row["n_immuno_cdr"]
    row["n_immuno_observed"] = len(observed_flags)

    soft = sum(f.weight for f in result.all_flags if f.severity == "soft")
    row["tier1_score"] = round(soft - tap_weight - immuno_weight, 3)

    if result.tap_profile:
        for metric in TAP_METRICS:
            key = metric.split()[0].lower()
            row["tap_" + key] = round(result.tap_profile.values.get(metric, float("nan")), 3)
            row["tap_flag_" + key] = result.tap_profile.flags.get(metric, "")

    try:
        chains = number_antibody(candidate["vh_sequence"], candidate["vl_sequence"])
        row["vh_v_gene"] = chains["VH"].v_gene or ""
        row["vl_v_gene"] = chains["VL"].v_gene or ""
        from ..humanness import germline_identity
        h = germline_identity(chains["VH"])
        if h.available and h.identity is not None:
            row["vh_germline_identity"] = round(h.identity, 4)
    except NumberingError:
        pass

    return row


def build_population(
    n_per_group: int, seed: int, include_non_human: bool
) -> list[dict]:
    df = load_paired_sequences()
    df = df[df["pairing"].isin([THERAPEUTIC_PAIRING, BACKGROUND_PAIRING])]
    df = df.dropna(subset=["heavy_sequence", "light_sequence"])
    df = df[
        (df["heavy_sequence"].str.len() >= 90) & (df["light_sequence"].str.len() >= 90)
    ]
    df = df.drop_duplicates(subset=["heavy_sequence", "light_sequence"])

    population = []
    for pairing in (THERAPEUTIC_PAIRING, BACKGROUND_PAIRING):
        group = df[df["pairing"] == pairing]
        # Oversample then filter: the human-framework fraction is ~70-75%, so
        # draw generously and stop once the quota is met.
        pool = group.sample(n=min(len(group), n_per_group * 3), random_state=seed)
        taken = 0
        for record in pool.itertuples():
            if taken >= n_per_group:
                break
            if not include_non_human and not is_human_framework(
                record.heavy_sequence, record.light_sequence
            ):
                continue
            population.append(
                {
                    "candidate_id": str(record.ID),
                    "vh_sequence": record.heavy_sequence,
                    "vl_sequence": record.light_sequence,
                    "pairing": pairing,
                }
            )
            taken += 1
        print(f"  {pairing}: {taken} candidates", file=sys.stderr)
    return population


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-per-group", type=int, default=150)
    parser.add_argument("--seed", type=int, default=0)
    # 3 rather than 5: each worker holds torch + ABodyBuilder2 weights +
    # OpenMM, and five of those exhausted a 16GB machine mid-run.
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--out", type=Path, default=Path("calibration_results.csv"))
    parser.add_argument("--include-non-human", action="store_true")
    parser.add_argument(
        "--resume", action="store_true",
        help="skip candidates already present in --out (default off so a fresh "
             "run does not silently mix results from different code versions)",
    )
    parser.add_argument(
        "--max-attempts", type=int, default=4,
        help="how many times to rebuild the pool after a worker dies",
    )
    parser.add_argument(
        "--model-cache",
        default=str(Path.home() / ".cache" / "evalab" / "model_cache"),
    )
    args = parser.parse_args()

    print("Building population...", file=sys.stderr)
    population = build_population(args.n_per_group, args.seed, args.include_non_human)

    # Resume support, and the reason it exists: a 300-candidate run holds
    # torch + ABodyBuilder2 weights + OpenMM in every worker, and on a 16GB
    # machine five of those is enough to get one killed. The first attempt at
    # this reached 200/300 and then died with BrokenProcessPool, losing
    # everything because results were only written at the end. Rows are now
    # flushed as they complete and finished candidates are skipped on restart.
    done: set[str] = set()
    if args.out.exists() and args.resume:
        with open(args.out, newline="") as handle:
            done = {row["candidate_id"] for row in csv.DictReader(handle)}
        print(f"Resuming: {len(done)} already scored", file=sys.stderr)

    remaining = [c for c in population if c["candidate_id"] not in done]
    print(
        f"Scoring {len(remaining)} candidates on {args.workers} workers...",
        file=sys.stderr,
    )

    write_header = not (args.out.exists() and done)
    with open(args.out, "a" if done else "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        if write_header:
            writer.writeheader()
            handle.flush()

        completed = len(done)
        # A dead pool kills its pending futures, not the work itself — rebuild
        # and carry on with whatever is left rather than discarding the run.
        for attempt in range(1, args.max_attempts + 1):
            if not remaining:
                break
            if attempt > 1:
                print(
                    f"  pool died; attempt {attempt} on {len(remaining)} remaining",
                    file=sys.stderr,
                )
            failed_this_pass = list(remaining)
            try:
                with ProcessPoolExecutor(max_workers=args.workers) as pool:
                    futures = {
                        pool.submit(score_one, candidate, args.model_cache): candidate
                        for candidate in remaining
                    }
                    for future in as_completed(futures):
                        candidate = futures[future]
                        row = future.result()
                        writer.writerow(row)
                        handle.flush()  # survive the next crash
                        failed_this_pass.remove(candidate)
                        completed += 1
                        if completed % 10 == 0:
                            print(f"  {completed}/{len(population)}", file=sys.stderr, flush=True)
            except (BrokenExecutor, OSError) as e:
                print(f"  pool failure: {type(e).__name__}: {e}", file=sys.stderr)
            remaining = failed_this_pass

    if remaining:
        print(
            f"WARNING: {len(remaining)} candidates never completed after "
            f"{args.max_attempts} attempts — rerun with --resume to retry them.",
            file=sys.stderr,
        )
    print(f"Wrote results to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
