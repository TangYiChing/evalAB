# Tier 2 plan: structure-based TAP profile + PLAbDab calibration

Status: **plan only — nothing in this document has been implemented yet.**
Everything marked "verified" below was actually run/probed during planning on
this machine (Apple M5, arm64, 10 cores, macOS 25.6); everything else is
marked as an estimate to be measured.

## Goal

Add a Tier 2 (structure-based) stage to `antibody_prescreen`:

```
VH/VL sequence  ->  ABodyBuilder2 (ImmuneBuilder)  ->  IMGT-numbered .pdb
                ->  TAP 5 metrics  ->  GREEN/AMBER/RED flags
                ->  fused into the existing GO / CONDITIONAL / NO-GO verdict
```

and re-calibrate the fusion thresholds against PLAbDab instead of the current
508-structure SAbDab-derived set.

---

## Findings from recon (these change the shape of the plan)

### 1. The PLAbDab 5 GB download is avoidable

`plabdab_data.tar.gz` is **5,045,787,915 bytes** (verified via HTTP HEAD).
Streaming it and aborting early does **not** work: the tar member order puts
`data-dev/models/*.pdb` (~65k ABodyBuilder2 models, ~280 KB each) **first**,
and the sequence CSVs come after them (verified by streaming the first 300 MB
and listing members).

But the individual files are served standalone from the same directory:

| File | Size | Verified |
|---|---|---|
| `paired_sequences.csv.gz` | **11,392,439 B (11 MB)** | downloaded + parsed OK |
| `unpaired_sequences.csv.gz` | 43,960,421 B (42 MB) | HEAD 200 |
| `config.json` | — | 404 (not needed for CSV-only use) |
| full `plabdab_data.tar.gz` | 5.05 GB | HEAD 200 |

```
https://opig.stats.ox.ac.uk/webapps/plabdab/static/downloads/paired_sequences.csv.gz
```

`paired_sequences.csv.gz` parses to **176,894 rows**, columns:

```
ID, heavy_sequence, light_sequence, heavy_ID, light_ID,
heavy_definition, light_definition, organism,
reference_authors, reference_title, update_date,
cdr_lengths, model, pairing, targets_mentioned
```

**Implication:** we do not need the `PLAbDab` python package, KA-Search, or
the 5 GB archive for calibration. We need `pandas.read_csv` on an 11 MB file.
The PLAbDab package is only required for its *search* features
(`sequence_search` needs the `kasearch_db/` index; `structure_search` needs
`models/`) — neither is on the critical path for "calibrate thresholds".
Keep the package as an optional Stage 4 extra, not a Stage 1 dependency.

### 2. PLAbDab ships its own positive-control label

The `pairing` column partitions the 176,894 rows by provenance:

```
Patent text        91024
Same entry         34983
Unique source      22379
Xtal structure     12778
Ordered entries     9062
Unique chain        4218
Unique word         1252
TheraSAbDab         1198
```

`pairing == "TheraSAbDab"` is **1,198 clinical-stage therapeutics** — i.e. the
current TheraSAbDab test set is a labelled *subset* of PLAbDab. That gives a
two-population calibration for free:

- **enriched-developable**: `pairing == "TheraSAbDab"` (survived far enough to
  get an INN and enter the clinic)
- **unselected background**: `pairing == "Patent text"` (91k sequences that
  were merely patented — no developability filter applied)

A threshold that does not separate these two populations is not measuring
developability. This is a *weak* label (clinical advancement ≠ developability,
and patented ≠ bad), but it is a far better calibration signal than "a crystal
structure exists", which is what `fusion.py` currently uses.

### 3. TAP does not need Rosetta or ChimeraX

`Exscientia/ab-characterisation`'s **top-level** pipeline needs Rosetta +
ChimeraX + mpi4py. But the TAP subtree is self-contained:

```
src/ab_characterisation/developability_tools/tap/
  main.py                    -> run_tap(modelfile, outfile, quiet) -> list[MetricResult]
  definitions.py             -> IMGT CDR ranges, anchors, hydrophobicity, charges
  structure_annotation.py    -> StructureAnnotator (SASA via psa, neighbours, salt bridges)
  metrics/{hydrophobic,positive,negative}_patches.py, sfvcsp.py, total_cdr_length.py
  psa_executables/{psa,psa_mac}
  outputs.py
```

Its only real dependencies are **biopython** and **loguru** plus the bundled
`psa` binary. `run_tap()` takes an already-IMGT-numbered PDB — which is exactly
what `ABodyBuilder2.predict().save()` produces.

Licence: **MIT** (Exscientia, 2024) — vendoring with attribution is permitted.
ImmuneBuilder and PLAbDab are BSD-3-Clause.

### 4. The `psa` binary runs on Apple Silicon

`psa_executables/psa_mac` is `Mach-O 64-bit executable x86_64`. **Verified:
it executes on this arm64 M5 under Rosetta 2** and prints its usage banner
(`psa version 2.0 (Feb 2000)`). `structure_annotation.py` already selects it
automatically via `sys.platform == "darwin"`.

Prerequisite: Rosetta 2 must be installed (`softwareupdate --install-rosetta
--agree-to-license`). If a future macOS drops Rosetta 2, the fallback is to
run TAP remotely (Step R below) — there is no arm64 `psa` build.

### 5. Consequence for the "do I need a GPU / remote box?" question

**No, not for the pilot.** ABodyBuilder2 is a small model; its slow step is
OpenMM refinement, which is CPU. `openmm 8.5.2` and `torch 2.13.0` (MPS
available) are already present in the base python on this machine. A pilot of
50–500 antibodies is a local, single-machine job.

Remote is needed only for the **full** calibration sweep (tens of thousands of
ABodyBuilder2 models). That is Step R, kept deliberately separate so the
pipeline can be proven correct before anyone pays for compute.

Timing per antibody is **not measured yet** — Step 6 exists to measure it
before committing to any batch size.

---

## Target architecture

```
antibody_prescreen/
  numbering.py            (unchanged)
  checks.py               (unchanged)   <- Tier 1 sequence liabilities
  developability.py       (unchanged)
  immunogenicity.py       (unchanged)
  fusion.py               (MODIFIED: accepts optional Tier 2 flags)
  structure/                              <- NEW
    __init__.py
    modelling.py          ABodyBuilder2 wrapper + on-disk model cache
    tap_runner.py         thin adapter: pdb path -> list[Flag]
    vendor/
      tap/                vendored from Exscientia/ab-characterisation (MIT)
      VENDOR.md           commit SHA, licence text, list of local edits
  calibration/                            <- NEW
    plabdab_source.py     download + parse paired_sequences.csv.gz, sample subsets
    run_calibration.py    score a population, dump distributions
  data/
    sabdab_derived_candidates.json        (kept, existing)
    plabdab_pilot.json                    NEW: the small pilot subset
```

**Design rule:** Tier 2 stays optional and lazily imported. `screen_batch(...)`
with no structure flag must keep working with zero new dependencies, so the
existing 15 tests and `examples/screen_candidates.py` are unaffected.

---

## Step-by-step

### Step 0 — decide where the code lives (5 min)

`gh auth status` on this machine is logged in as `TangYiChing` with `repo`
scope, and `gh repo view TangYiChing/evalAB` reports
`viewerPermission: ADMIN`. The clone at `antibody_evaluation/evalAB` has
`origin` pointing at the GitHub repo. So: work in the local clone, commit,
and `git push` / open a PR. No separate copy needed.

Current branch is `claude/pensive-dijkstra-gylwgr`, which is also the repo's
default branch. Work on a new branch:

```bash
cd ~/Projects/antibody_evaluation/evalAB
git checkout -b tier2-tap
```

### Step 1 — create the conda environment (20 min, one-off)

Rosetta 2 first (needed for the `psa` binary):

```bash
softwareupdate --install-rosetta --agree-to-license
```

Then the env. Pin python 3.10 to match what ab-characterisation and ANARCI
are tested against:

```bash
conda create -n evalab python=3.10 -y
conda activate evalab

# Tier 1 (existing pipeline)
conda install -c bioconda -c conda-forge hmmer anarci -y

# Tier 2 (new)
conda install -c conda-forge biopython pandas numpy openmm pdbfixer loguru pytest -y
pip install torch                 # CPU/MPS wheel is fine, no CUDA on a Mac
pip install ImmuneBuilder

# repo itself
pip install -e .                  # if/when a pyproject.toml is added; else just run from the repo root
```

Verify before going further:

```bash
python -c "import anarci, Bio, torch, openmm, loguru; print('deps ok')"
python -c "from ImmuneBuilder import ABodyBuilder2; ABodyBuilder2(); print('ABB2 weights ok')"
```

The second command downloads ABodyBuilder2 weights on first run (~100 MB to
`~/.cache`). Do it now so it is not a surprise mid-batch.

Write the exact resolved versions into `environment.yml` in the repo so this
is reproducible — do not leave it as prose in a README.

### Step 2 — vendor the TAP subtree (30 min)

Copy `src/ab_characterisation/developability_tools/tap/` into
`antibody_prescreen/structure/vendor/tap/`, then:

1. Rewrite the absolute imports (`from ab_characterisation.developability_tools.tap...`)
   to relative ones (`from .definitions import ...`). This is mechanical —
   about 8 files.
2. Keep both `psa_executables/psa` (Linux) and `psa_mac` (macOS) — the
   platform switch in `structure_annotation.py.__post_init__` already handles
   the choice, and keeping both is what makes Step R work unchanged.
3. `chmod +x` both binaries and make sure git records the executable bit
   (`git update-index --chmod=+x`). A vendored binary that loses its `+x` bit
   fails at runtime with the unhelpful `PSAError("psa executable was not found.")`.
4. Write `VENDOR.md` recording: upstream URL, the exact commit SHA pinned, the
   MIT licence text, and the list of edits made (import rewrites only).

Why vendor rather than `pip install` upstream: `ab-characterisation` has no
PyPI release, and installing it from git drags in the Rosetta/ChimeraX/mpi4py
side of the package that we do not want. Vendoring one MIT subtree with a
recorded SHA is the smaller commitment.

Checkpoint:

```bash
python -c "from antibody_prescreen.structure.vendor.tap.main import run_tap, list_metrics; list_metrics()"
```

Expected output is the 5 metrics with their green/amber regions:

| Metric | GREEN | AMBER |
|---|---|---|
| Total IMGT CDR Length | 43 – 55 | 37–43, 55–63 |
| Hydrophobic Patch (PSH) | 137.61 – 200.71 | 106.44–137.61, 200.71–225.85 |
| Positive Patch (PPC) | 0 – 1.19 | 1.19 – 3.58 |
| Negative Patch (PNC) | 0 – 1.67 | 1.67 – 3.50 |
| SFvCSP | ≥ −4.20 | −20.50 – −4.20 |

Anything outside green or amber is RED. Note these regions were set by
Raybould et al. against a therapeutic-antibody population — they are already
calibrated, which is precisely what the Tier 1 thresholds are not.

### Step 3 — `structure/modelling.py` (1 h)

```python
def build_model(candidate_id, vh, vl, cache_dir) -> Path
```

- returns the cached `.pdb` if `cache_dir/{candidate_id}.pdb` exists — model
  building is the expensive step and calibration will re-run many times
- otherwise `ABodyBuilder2().predict({"H": vh, "L": vl})`, then `.save(path)`
- catches ImmuneBuilder failures and raises a typed `ModellingError` so the
  caller can degrade to Tier-1-only for that candidate instead of crashing
  the batch — same graceful-degradation contract `immunogenicity.py` already
  follows with `available=False`

The cache directory belongs in `.gitignore`; PDBs are build artefacts.

### Step 4 — `structure/tap_runner.py` (1 h)

Adapter from TAP's `MetricResult` to the existing `Flag` dataclass, so
`fusion.py` does not need to learn a second result type:

- `flag == "RED"` → `Flag(severity="soft", weight=RED_WEIGHT, ...)`
- `flag == "AMBER"` → `Flag(severity="soft", weight=AMBER_WEIGHT, ...)`
- `flag == "GREEN"` → no flag

**Do not make any TAP RED a hard gate.** The Tier 1 cysteine story in the
README is the precedent: a blunt hard gate was already downgraded once
because real approved antibodies trip it. TAP RED means "outside the range of
known therapeutics on one axis" — several approved drugs sit in amber/red on
at least one metric. Weighted soft flags keep that recoverable.

Put `RED_WEIGHT` / `AMBER_WEIGHT` as named module constants with a comment
saying they are **provisional until Step 7 sets them from data** — not
hand-tuned numbers hidden in a function body.

### Step 5 — wire into `fusion.py` (1 h)

Add an opt-in parameter, mirroring the existing `run_immunogenicity` flag:

```python
def screen_candidate(candidate_id, vh_sequence, vl_sequence,
                     run_immunogenicity=True,
                     run_structure=False,        # NEW
                     model_cache=None):          # NEW
```

- import `structure.*` **inside** the `if run_structure:` branch, so a user
  without ImmuneBuilder installed is unaffected
- add `structure_available: bool` to `CandidateResult`, exactly like
  `immunogenicity_available`
- extend `format_report` with a TAP column only when at least one result has
  `structure_available` — do not widen the default table for everyone

Checkpoint: the existing test suite must still pass untouched.

```bash
python -m pytest tests/ -v
```

Then add new tests: a known-clean therapeutic pair should come back all-green
or near-all-green; a deliberately mutated variant (inject a hydrophobic patch
into CDR-H3) should move PSH toward red. Same mutation-injection style the
existing 15 tests already use.

### Step 6 — measure the cost before scaling (30 min)

Before any batch, time it honestly:

```bash
python -m antibody_prescreen.calibration.run_calibration --n 20 --time-only
```

Record wall-clock per antibody for (a) ABB2 predict, (b) OpenMM refine,
(c) TAP. Multiply by the calibration set size. This number decides Step R.
Do not guess it.

### Step 7 — PLAbDab pilot and calibration (half a day)

**7a. Fetch the sequences (seconds, not a 5 GB download):**

```python
# antibody_prescreen/calibration/plabdab_source.py
PLABDAB_PAIRED_URL = (
    "https://opig.stats.ox.ac.uk/webapps/plabdab/static/downloads/"
    "paired_sequences.csv.gz"
)
```

Cache it to `~/.cache/evalab/plabdab/paired_sequences.csv.gz` (11 MB), never
into the repo.

**7b. Build the pilot set (~100 antibodies):**

- 50 from `pairing == "TheraSAbDab"` (enriched-developable)
- 50 from `pairing == "Patent text"` (unselected background)
- drop rows with a null/short `heavy_sequence` or `light_sequence`
- drop rows ANARCI cannot number
- deduplicate on `(heavy_sequence, light_sequence)` — `Same entry` and
  `Unique source` pairings mean the same antibody can appear more than once
- fix the random seed and commit the resulting **ID list** (not the sequences)
  to `data/plabdab_pilot.json` so the pilot is reproducible

100 antibodies is deliberately small: the point of the pilot is to prove the
pipeline end-to-end (sequence → model → TAP → fused verdict) without waiting
on a long run. Correctness first, statistics second.

**7c. Run both tiers over the pilot** and dump per-candidate: the 5 TAP raw
values, the 5 flags, the Tier 1 score, and the fused score.

**7d. Look at the two distributions.** The questions to answer, in order:

1. Does the pipeline run to completion on 100 real sequences without a crash,
   and what is the ANARCI/ABB2 failure rate?
2. Do the TheraSAbDab entries sit greener than the patent-text entries on the
   5 TAP metrics? If *not*, stop — either the runner is wrong, or the label is
   weaker than assumed. Do not proceed to retune thresholds on a signal that
   is not there.
3. Only if (2) holds: where do `T_LOW` / `T_HIGH` and the TAP weights have to
   sit to keep most TheraSAbDab entries at GO while pushing the worst of the
   background to NO-GO?

**7e. Scale up** to a few thousand per group once the pilot is clean, and
re-fit. Record the population, the date, and the seed in `fusion.py`'s
threshold comment — the existing comment block is a good template for this
and should be updated, not appended to.

### Step 8 — documentation and honesty pass (1 h)

Update `antibody_prescreen/README.md`:

- move "Tier 2 (structure-based checks)" out of "What's not done"
- add the Tier 2 setup block and the Rosetta-2 requirement
- state the new calibration population and, explicitly, what its label does
  **not** mean: "reached the clinic" is not "manufacturable", and "appears in
  a patent" is not "bad". This is a *weak, provenance-based* label and the
  README should say so as plainly as it currently says the SAbDab set is
  selected for crystallisability.
- keep the existing note that Tier 2 does not resolve the cysteine-pairing
  question unless the disulfide check is actually implemented against the
  model — the README currently promises "deferred to a future Tier 2 structure
  check", so either implement it or restate the deferral

### Step 9 — ship

```bash
git add -A && git commit -m "Add Tier 2 structure screening: ABodyBuilder2 -> TAP profile"
git push -u origin tier2-tap
gh pr create --fill
```

---

## Step R — remote execution (only if Step 6 says local is too slow)

Trigger condition: measured per-antibody wall time × target set size exceeds
what is acceptable locally. For a 100-antibody pilot this will almost
certainly not trigger. For a 20,000-antibody calibration it will.

This job is **CPU-parallel, not GPU-bound**: ABodyBuilder2's network is small
and its bottleneck is OpenMM refinement plus the single-threaded `psa`
subprocess call per structure. So the right remote shape is *many CPU cores*,
not an expensive GPU.

**Machine shape:** a 32–64 vCPU Linux VM (or a SLURM array job on a cluster).
GPU optional; it would only speed the predict step, not refinement or TAP.

**Environment — identical recipe, one substitution:**

```bash
conda create -n evalab python=3.10 -y && conda activate evalab
conda install -c bioconda -c conda-forge hmmer anarci -y
conda install -c conda-forge biopython pandas numpy openmm pdbfixer loguru pytest -y
pip install torch ImmuneBuilder
```

On Linux the vendored `psa` (not `psa_mac`) is selected automatically by
`sys.platform`, so no code change is needed — this is why Step 2 keeps both
binaries. Verify the Linux binary is executable after checkout:

```bash
chmod +x antibody_prescreen/structure/vendor/tap/psa_executables/psa
./antibody_prescreen/structure/vendor/tap/psa_executables/psa   # should print the usage banner
```

**Where the data goes:**

| What | Where | Size |
|---|---|---|
| `paired_sequences.csv.gz` | `$SCRATCH/plabdab/` | 11 MB |
| ABodyBuilder2 weights | `~/.cache/` (first run) | ~100 MB |
| generated `.pdb` models | `$SCRATCH/model_cache/` | ~280 KB each |
| results CSV | rsync'd back to the Mac | small |

Model cache sizing: at ~280 KB per model (measured from the PLAbDab archive's
own ABB2 models), 20,000 candidates ≈ **5.6 GB**. Provision scratch for it, or
add a `--no-cache` mode that deletes each PDB after TAP reads it. Do **not**
download the 5 GB PLAbDab tarball to the remote box just to get models —
building them is cheaper than transferring them, and PLAbDab's models may not
be from the same ABB2 version.

**Sharding:** split the input CSV into N chunks, run one process per chunk
(`--shard i --of N`), concatenate the result CSVs. Each candidate is
independent, so this is embarrassingly parallel — no MPI needed despite
upstream's `mpi4py` dependency, which exists for the Rosetta stage we are not
using.

**Bring back only the results CSV.** Models stay remote unless a specific one
needs inspecting.

---

## What this plan does not do

- **Does not validate TAP against wet-lab outcomes.** TAP's own thresholds
  come from a therapeutic-antibody population; PLAbDab's `pairing` column is
  a provenance label, not a measured developability assay. Neither is an
  experimental ground truth. The pipeline gets *better* calibrated, not
  *validated*.
- **Does not use PLAbDab's KA-Search or structure search.** Both need the
  5 GB archive. If germline-precedent or nearest-therapeutic lookup becomes
  wanted (it is listed as an outstanding bonus signal in the README), that is
  a separate decision with a separate download.
- **Does not solve the immunogenicity gap.** IEDB is still unreachable per the
  existing README; Tier 2 does not change that.
- **Does not add a disulfide-pairing check** even though a structure is now
  available. Worth doing, but it is its own piece of work with its own tests.
