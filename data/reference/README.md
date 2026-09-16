# Reference data

The reference population is the only thing in evalAB that decides what
"unusual" means. Its selection bias is inherited by every band fitted on it,
so it is pinned here rather than downloaded silently at runtime.

## `paired_sequences.csv.gz` — PLAbDab paired sequences

| | |
|---|---|
| Source | Patented and Published Antibody Database (PLAbDab), OPIG, University of Oxford |
| URL | `https://opig.stats.ox.ac.uk/webapps/plabdab/static/downloads/paired_sequences.csv.gz` |
| Size | 11,392,439 bytes (~11 MB) |
| SHA-256 | `10f0fb3b82ddee09209adaff51340380294e3782ba90034af16c9f1d56c458a5` |
| Retrieved | 2026-09-16 |
| Citation | Abanades et al., *PLAbDab: a database of printed and published antibodies*, Nucleic Acids Research (2024) |

This file is **not committed**. It is ~11 MB of third-party data with its own
terms, and redistributing it inside this repository would make evalAB a
mirror of someone else's database. Fetch it instead:

```bash
python -c "from antibody_prescreen.calibration import fetch_paired_sequences as f; print(f())"
```

It lands here, in `data/reference/`, and the checksum above is what you should
get. If it differs, PLAbDab has been updated since these bands were fitted —
which does not make either version wrong, but does mean the shipped bands and
your download no longer describe the same population. Refit before comparing.

Set `EVALAB_REFERENCE_DIR` to keep it somewhere else.

### Why not the 5 GB archive

PLAbDab's documented download is `plabdab_data.tar.gz` (~5.05 GB), almost all
of which is ~65k pre-built ABodyBuilder2 `.pdb` models. Those members come
first in the tar, so streaming and aborting early never reaches the sequence
CSVs. The standalone `paired_sequences.csv.gz` is served from the same
directory and is the only file evalAB needs.

### The `pairing` column is provenance, not an assay

Two of its values are used as weak, opposing populations:

- `TheraSAbDab` (~1.2k) — reached clinical stage
- `Patent text` (~91k) — appears in a patent, no developability filter

"Reached the clinic" is not "manufacturable" and "was patented" is not "bad".
Treating this contrast as an outcome label is exactly the mistake that
produced the AUC 0.587 result in `docs/range-triage-design.md` — the result
that rules out scoring and motivates range triage instead.

## What is deliberately absent

No proprietary or internal design sequences are in this repository, in any
directory, at any point in its history. Everything here is public data or
derived from it. The band files under `antibody_prescreen/data/` are the
derived artefacts, and each records its own reference, split sizes and cut
points.
