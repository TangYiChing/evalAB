"""Example: Tier 1 + Tier 2 (TAP) screening on a small PLAbDab pilot.

Run it:
    conda activate evalab
    python examples/screen_with_tap.py          # 6 candidates, ~30s cold
    python examples/screen_with_tap.py 20       # 20 candidates

What it does, per candidate:
    VH/VL sequence
      -> ABodyBuilder2 -> IMGT-numbered Fv model (.pdb, cached on disk)
      -> TAP -> 5 metrics, each GREEN / AMBER / RED
      -> weighted soft flags, fused with the Tier 1 score into one verdict

The candidates come from PLAbDab's paired_sequences.csv.gz (11 MB, fetched
and cached on first run — NOT the 5 GB archive; see
antibody_prescreen/calibration/plabdab_source.py for why). The pilot is
balanced between clinical-stage therapeutics and unselected patent entries,
which is the population split the threshold calibration will eventually use.

NOTE: the TAP weights in tap_runner.py are PROVISIONAL and uncalibrated, so
the fused score here is not yet meaningful as a go/no-go number. The five TAP
values and their flags, however, are real — those thresholds come from
Raybould et al.'s therapeutic-antibody population, not from us.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from antibody_prescreen import format_report, screen_batch
from antibody_prescreen.calibration import sample_pilot

MODEL_CACHE = Path.home() / ".cache" / "evalab" / "model_cache"


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    per_group = max(1, n // 2)

    candidates = sample_pilot(n_per_group=per_group, seed=0)
    pairing = {c["candidate_id"]: c["pairing"] for c in candidates}
    print(f"Screening {len(candidates)} candidates "
          f"({per_group} therapeutic / {per_group} patent-text)")
    print(f"Model cache: {MODEL_CACHE}\n")

    results = screen_batch(
        candidates,
        run_immunogenicity=False,   # IEDB unreachable in many environments
        run_structure=True,         # <- Tier 2
        model_cache=MODEL_CACHE,
    )

    print("=== Summary ===\n")
    print(format_report(results))

    print("\n=== TAP profile detail ===\n")
    header = f"{'Candidate':24s} {'Source':12s} " + " ".join(
        f"{h:>8s}" for h in ("CDRlen", "PSH", "PPC", "PNC", "SFvCSP")
    )
    print(header)
    print("-" * len(header))
    for r in results:
        if not r.structure_available:
            print(f"{r.candidate_id[:24]:24s} {pairing.get(r.candidate_id,'?')[:12]:12s} "
                  f"(no structure)")
            continue
        v = r.tap_profile.values
        f = r.tap_profile.flags
        cells = []
        for metric in (
            "Total IMGT CDR Length",
            "Hydrophobic Patch Score",
            "Positive Patch Score",
            "Negative Patch Score",
            "SFvCSP",
        ):
            mark = {"GREEN": " ", "AMBER": "*", "RED": "!"}[f[metric]]
            cells.append(f"{v[metric]:7.2f}{mark}")
        print(f"{r.candidate_id[:24]:24s} {pairing.get(r.candidate_id,'?')[:12]:12s} "
              + " ".join(cells))
    print("\n  ' ' = green    '*' = amber    '!' = red")

    flagged = [r for r in results if r.tap_profile and r.tap_profile.worst_flag != "GREEN"]
    if flagged:
        print(f"\n=== {len(flagged)} candidate(s) with a non-green TAP axis ===\n")
        for r in flagged:
            print(f"{r.candidate_id} ({r.verdict}, score={r.score}):")
            for flag in r.all_flags:
                if flag.check == "tap":
                    print(f"  w={flag.weight:<5} {flag.message}")
    else:
        print("\nAll candidates green on all five TAP axes.")


if __name__ == "__main__":
    main()
