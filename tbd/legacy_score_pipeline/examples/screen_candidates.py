"""Example: screen a batch of antibody candidates with the Tier 1 pipeline.

Run it:
    python3 examples/screen_candidates.py

Uses 10 real candidates from the committed sabdab_derived_candidates.json
(see antibody_prescreen/sabdab_source.py for how that was extracted) so this
runs standalone, no network needed. Swap in your own VH/VL sequences by
building the same candidates list shape shown below.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from antibody_prescreen import screen_batch, format_report

DATA_FILE = Path(__file__).parent.parent / "antibody_prescreen" / "data" / "sabdab_derived_candidates.json"


def main():
    with open(DATA_FILE) as f:
        sabdab_candidates = json.load(f)

    # Build the {candidate_id, vh_sequence, vl_sequence} shape screen_batch
    # expects. Your own candidates (from a design tool, a spreadsheet, etc.)
    # just need to match this same shape - antigen_sequence and the other
    # SAbDab-specific fields aren't used by Tier 1, which is sequence-only.
    candidates = [
        {
            "candidate_id": c["pdb_id"],
            "vh_sequence": c["vh_sequence"],
            "vl_sequence": c["vl_sequence"],
        }
        for c in sabdab_candidates[:10]
    ]

    # run_immunogenicity=False: IEDB's API isn't reachable from every
    # environment (see antibody_prescreen/README.md "Data source" and
    # immunogenicity.py) - the check degrades gracefully either way, but
    # explicitly disabling it here skips the network calls for a fast demo.
    results = screen_batch(candidates, run_immunogenicity=False)

    print("=== Summary ===\n")
    print(format_report(results))

    # Full detail for the first candidate that needs work. format_report only
    # shows the top 1-2 reasons by design (that's the whole point - one line
    # per candidate for triage). Drill into result.all_flags for the complete
    # picture before acting on a level.
    print("\n=== Full detail for the first candidate above Level 1 ===\n")
    for r in results:
        if r.triage_result is not None and r.triage_result.level > 1:
            print(f"{r.candidate_id} (Level {r.triage_result.level}: "
                  f"{r.triage_result.action}):")
            for flag in r.all_flags:
                print(f"  [{flag.severity:9s}] {flag.check:16s} {flag.region:10s} "
                      f"w={flag.weight:<5} {flag.message}")
            break


if __name__ == "__main__":
    main()
