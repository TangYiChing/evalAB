"""Build the reference band file, and measure what it does on unseen antibodies.

Two things happen here and they must not be confused:

  DERIVATION   half the reference set -> the percentile cut points in bands.json
  VALIDATION   the other half -> the level distribution reported at the end

The validation numbers are the only ones that mean anything. A rejection rate
measured on the same antibodies the cut points were fitted to is a restatement
of where the lines were drawn.

Usage:
    python -m antibody_prescreen.calibration.fit_bands --reference therasabdab
    python -m antibody_prescreen.calibration.fit_bands --reference plabdab --limit 8000
"""

import argparse
import json
import sys
from pathlib import Path

from anarci import run_anarci

from ..bands import CLINICAL_METRICS, BandSet, compose, fit_bands
from ..levels import assign
from ..numbering import IMGT_REGIONS, NumberedChain, Residue
from ..profile import build_profile

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

REFERENCES = {
    # name -> (pairing filter or None, human-readable label)
    "therasabdab": (["TheraSAbDab"], "PLAbDab/TheraSAbDab (clinical-stage)"),
    "plabdab": (None, "PLAbDab (all pairings)"),
}


def _region_for(position: int) -> str:
    for name, start, end in IMGT_REGIONS:
        if start <= position <= end:
            return name
    return "FR4"


def number_many(sequences: list[str], expect_heavy: bool) -> list[NumberedChain | None]:
    """One ANARCI call for the whole list.

    `numbering.number_chain` runs ANARCI once per sequence, which is right for
    screening a handful of candidates and wrong for numbering several thousand:
    each call re-launches hmmscan. Same numbering, same scheme, batched.
    """
    out: list[NumberedChain | None] = [None] * len(sequences)
    _, numbered, details, _ = run_anarci(
        [(str(i), s) for i, s in enumerate(sequences)],
        scheme="imgt", assign_germline=True, ncpu=4,
    )
    for i, domains in enumerate(numbered):
        if not domains or len(domains) != 1:
            continue
        d = details[i][0]
        if (d["chain_type"] == "H") != expect_heavy:
            continue
        residues = []
        for (position, insertion), aa in domains[0][0]:
            region = _region_for(position)
            residues.append(Residue(position, insertion.strip(), aa, region, region.startswith("CDR")))
        v_gene = v_species = v_identity = None
        call = (d.get("germlines") or {}).get("v_gene")
        if call:
            (v_species, v_gene), v_identity = call
        out[i] = NumberedChain(
            chain_type=d["chain_type"], species=d["species"], evalue=d["evalue"],
            residues=residues, v_gene=v_gene, v_species=v_species, v_identity=v_identity,
        )
    return out


def load_reference(name: str, limit: int | None, seed: int = 0):
    """Deduplicated (heavy, light) pairs from PLAbDab, shuffled deterministically."""
    from .plabdab_source import MIN_CHAIN_LEN, load_paired_sequences

    pairings, label = REFERENCES[name]
    df = load_paired_sequences()
    if pairings is not None:
        df = df[df["pairing"].isin(pairings)]
    df = df.dropna(subset=["heavy_sequence", "light_sequence"])
    df = df[(df["heavy_sequence"].str.len() >= MIN_CHAIN_LEN)
            & (df["light_sequence"].str.len() >= MIN_CHAIN_LEN)]
    df = df.drop_duplicates(subset=["heavy_sequence", "light_sequence"])
    df = df.sample(frac=1.0, random_state=seed)
    if limit:
        df = df.head(limit)
    return [(str(r.ID), r.heavy_sequence, r.light_sequence) for r in df.itertuples()], label


# Observations that disqualify an entry from being part of a REFERENCE
# population. Not from being a candidate — a candidate with a substituted fold
# anchor is a Level 4 or 5 candidate and the pipeline should say so. But a
# reference population is supposed to answer "what do antibodies look like",
# and a truncated or malformed database row is not an antibody that looks like
# anything. Leaving them in widens every band by exactly the amount the junk
# extends the tails.
#
# Measured before this filter existed, on 20,000 validation-half PLAbDab rows:
# 2.31% tripped `non_standard_aa` and 1.66% `anchor_substituted`. Those are
# database artefacts, not biology.
#
# DELIBERATELY NOT IN THIS SET: `odd_cysteine`. An odd cysteine count is
# unusual biology, not a malformed record, and it is the one observation in
# this repository with a measured association to a wet-lab outcome
# (docs/range-triage-design.md §7). Excluding it from the reference would
# quietly assume the conclusion it is evidence for, and would also bias
# `fv_cys_count`, which is a fitted band.
REFERENCE_DISQUALIFYING = {
    "non_standard_aa",
    "anchor_missing",
    "anchor_substituted",
}


def _is_human(chain) -> bool:
    """ANARCI's germline assignment says human for this chain.

    `v_species` is the species of the closest germline V gene, which is the
    question being asked. `species` (the HMM's own call) is the fallback for
    the rare chain ANARCI numbers without assigning a germline.

    This is a germline call, not a provenance record: a humanised antibody with
    murine back-mutations passes it. That is a limitation of the filter and a
    reason germline identity keeps its own reference population — see
    bands.CLINICAL_METRICS.
    """
    call = (chain.v_species or chain.species or "").lower()
    return call == "human"


def profile_all(rows, human_only: bool = False, clean: bool = False):
    """Number and profile a list of (id, heavy, light).

    human_only / clean are REFERENCE-population filters. They default off, so
    screening a candidate batch never silently drops a candidate — a candidate
    that fails these is exactly the candidate the triage exists to flag.
    """
    heavies = number_many([h for _, h, _ in rows], expect_heavy=True)
    lights = number_many([l for _, _, l in rows], expect_heavy=False)
    profiles = []
    dropped = {"unnumbered": 0, "non_human": 0, "malformed": 0}
    for (cid, _, _), vh, vl in zip(rows, heavies, lights):
        if vh is None or vl is None:
            dropped["unnumbered"] += 1
            continue
        if human_only and not (_is_human(vh) and _is_human(vl)):
            dropped["non_human"] += 1
            continue
        profile = build_profile(cid, {"VH": vh, "VL": vl})
        if clean and any(o.name in REFERENCE_DISQUALIFYING for o in profile.observations):
            dropped["malformed"] += 1
            continue
        profiles.append(profile)
    return profiles, dropped


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", choices=sorted(REFERENCES), default="therasabdab")
    ap.add_argument("--limit", type=int, default=None,
                    help="cap the reference set (use for a fast smoke run)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--clinical-bands", type=Path, default=None,
                    help="compose with a clinical band file before measuring, "
                         "taking bands.CLINICAL_METRICS from it — so the "
                         "validation rates below describe the set that will "
                         "actually be used, not a set nobody runs.")
    ap.add_argument("--human-only", action="store_true",
                    help="reference filter: both chains must call human germline")
    ap.add_argument("--clean", action="store_true",
                    help="reference filter: drop entries with a non-standard "
                         "amino acid or a broken IMGT fold anchor — malformed "
                         "database rows, not antibodies. See "
                         "REFERENCE_DISQUALIFYING.")
    args = ap.parse_args()

    rows, label = load_reference(args.reference, args.limit)
    # The filters change what the population IS, so they belong in the label
    # rather than only in this run's stderr. Every report prints this string.
    applied = [n for n, on in (("human-only", args.human_only), ("cleaned", args.clean)) if on]
    if applied:
        label = f"{label}, {', '.join(applied)}"
    print(f"reference: {label} — {len(rows)} deduplicated pairs", file=sys.stderr)

    half = len(rows) // 2
    derivation, validation = rows[:half], rows[half:]

    print("numbering derivation half...", file=sys.stderr)
    dev_profiles, dev_dropped = profile_all(
        derivation, human_only=args.human_only, clean=args.clean)
    print("numbering validation half...", file=sys.stderr)
    val_profiles, val_dropped = profile_all(
        validation, human_only=args.human_only, clean=args.clean)
    print(f"numbered: {len(dev_profiles)} derivation, {len(val_profiles)} validation",
          file=sys.stderr)
    for name, d in (("derivation", dev_dropped), ("validation", val_dropped)):
        total = sum(d.values())
        if total:
            print(f"  {name} dropped {total}: " +
                  ", ".join(f"{k}={v}" for k, v in d.items() if v), file=sys.stderr)

    observations: dict[str, list[float]] = {}
    tails: dict[str, str] = {}
    for p in dev_profiles:
        for m in p.metrics:
            observations.setdefault(m.name, []).append(m.value)
            tails[m.name] = m.tail

    band_set = fit_bands(
        observations, tails, reference=label,
        n_derivation=len(dev_profiles), n_validation=len(val_profiles),
        notes=("Cut points fitted on the derivation half ONLY. Every rate in "
               "the calibration report is measured on the validation half."),
    )
    out = args.out or DATA_DIR / f"bands_{args.reference}.json"
    out.write_text(band_set.to_json())
    print(f"wrote {out} — {len(band_set.bands)} bands", file=sys.stderr)

    if args.clinical_bands:
        band_set = compose(band_set, BandSet.load(args.clinical_bands),
                           CLINICAL_METRICS)
        print(f"composed: {band_set.reference}", file=sys.stderr)

    # --- what the bands DO, measured on antibodies they have never seen ----
    levels: dict[int, int] = {}
    driver_counts: dict[str, int] = {}
    for p in val_profiles:
        t = assign(p, band_set)
        levels[t.level] = levels.get(t.level, 0) + 1
        if t.level >= 4:
            for f in t.drivers:
                driver_counts[f.source] = driver_counts.get(f.source, 0) + 1

    n = len(val_profiles)
    print("\n=== VALIDATION HALF: level distribution ===")
    for lvl in sorted(levels):
        print(f"  L{lvl}  {levels[lvl]:5d}  {100*levels[lvl]/n:5.1f}%")
    # TWO budgets, because there are two different mistakes.
    #
    # The 5% figure is borrowed from trauma under-triage, where the mistake
    # being bounded is sending a severely injured patient home. Its analogue
    # here is L5 — "do not pursue" — because that is the only level that
    # discards a candidate. L4 hands the candidate to a person; a person
    # looking at an antibody is a cost, not a loss, and bounding it at 5%
    # would be borrowing a number from a situation it does not describe.
    #
    # Both are printed. Which one is the pass condition is a decision for
    # whoever runs this, and it has to be made before the run, not after.
    for label, floor in (("L5 (discard)", 5), ("L4+L5 (discard or escalate)", 4)):
        k = sum(v for lvl, v in levels.items() if lvl >= floor)
        lo, hi = wilson(k, n)
        print(f"\n  {label:30s} {k:4d}/{n} = {100*k/n:5.2f}%  "
              f"(95% CI {100*lo:.2f}-{100*hi:.2f}%)  "
              f"upper bound vs 5%: {'PASS' if hi <= 0.05 else 'FAIL'}")
    print("\n=== What drove L4/L5 ===")
    for k, v in sorted(driver_counts.items(), key=lambda kv: -kv[1]):
        print(f"  {k:26s} {v:5d}  {100*v/n:5.2f}%")


def wilson(k: int, n: int, z: float = 1.96):
    """Wilson score interval. The upper bound is the number the budget is on."""
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    s = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5)
    return max(0.0, (c - s) / d), min(1.0, (c + s) / d)


if __name__ == "__main__":
    main()
