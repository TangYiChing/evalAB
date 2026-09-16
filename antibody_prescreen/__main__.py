"""Command-line entry point for evalAB."""

import argparse
import csv
import json
from pathlib import Path

from .bands import BandSet
from .screen import (
    DEFAULT_BANDS, DEFAULT_CLINICAL_BANDS, format_report, load_default_bands,
    provenance, screen_batch,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Range-and-repair triage for paired human VH/VL candidates.")
    parser.add_argument("input", type=Path, help="CSV with candidate_id,vh_sequence,vl_sequence")
    parser.add_argument("--out", type=Path, help="write output here instead of stdout")
    parser.add_argument("--json", action="store_true", help="emit full machine-readable findings")
    parser.add_argument("--shared-scaffold", action="store_true", help="add an introduced-only view for a declared congeneric variant batch")
    parser.add_argument("--run-structure", action="store_true", help="add optional ABodyBuilder2/TAP annotations; never changes Tier 1 level")
    parser.add_argument("--model-cache", type=Path, default=None)
    parser.add_argument("--bands", type=Path, default=DEFAULT_BANDS, help="range band file (default: supplied repertoire bands)")
    parser.add_argument("--clinical-bands", type=Path, default=DEFAULT_CLINICAL_BANDS, help="clinical bands for germline identity")
    args = parser.parse_args()

    with args.input.open(newline="") as handle:
        candidates = list(csv.DictReader(handle))
    bands = load_default_bands() if args.bands == DEFAULT_BANDS and args.clinical_bands == DEFAULT_CLINICAL_BANDS else BandSet.load(args.bands)
    results = screen_batch(candidates, band_set=bands, run_structure=args.run_structure, model_cache=args.model_cache)
    run = provenance(results, bands)
    text = (
        json.dumps({"run": run.to_dict(), "results": [r.to_dict() for r in results]}, indent=2)
        if args.json
        else format_report(results, shared_scaffold=args.shared_scaffold, run_provenance=run)
    )
    if args.out:
        args.out.write_text(text + "\n")
    else:
        print(text)


if __name__ == "__main__":
    main()
