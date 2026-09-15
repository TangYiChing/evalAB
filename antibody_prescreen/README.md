# Antibody Pre-Screening (Tier 1 + Tier 2) — Prototype

**Status: prototype. Routing is by triage LEVEL, not by a score.** The summed
score turned out not to support a go/no-go decision (no gap between the
buckets, AUC 0.587, a false-rejection rate that was defined rather than
measured), so it is now a diagnostic field that nothing branches on. See
`antibody_prescreen/triage.py` and `docs/triage_implementation_plan.md`.

- **Five triage levels.** Level 1 is the candidate you act on first, Level 5
  the one you do not pursue. Note this is the INVERSE of CTAS/ESI severity,
  where Level 1 is the resuscitation patient: an emergency department spends
  its scarce resource on the sickest, drug discovery spends bench time on the
  cleanest, and both mappings cannot be preserved at once.

  | Level | Name | Condition |
  |---|---|---|
  | 1 | **Ready** | No findings, or only framework-located repairs. One round of point mutagenesis at most; binding untouched, no re-measurement. |
  | 2 | **Rework** | A repairable liability sits in a CDR or CDR vicinity. Fixable, but the fix may take the binding with it — re-measure affinity. |
  | 3 | **Immune risk** | A liability with no repair locus bearing on immunogenicity: substantially non-human framework, or Fv-wide charge asymmetry. A re-engineering decision, not a repair. |
  | 4 | **Hold** | A structural finding that is serious but unresolved — an unpaired cysteine in the predicted fold. A person must decide. |
  | 5 | **Out of range** | Outside the range spanned by every known therapeutic. The flag names which boundary and by how much. |

  **The level is the worst finding, not the sum of findings** — an emergency
  department triages a patient with a broken finger and chest pain on the chest
  pain. That is also why there are no weights in `triage.py`: a max needs none.

  A check may only route to Level 5 if the **upper** 95% bound of its
  false-rejection rate on known-good antibodies clears 5% — the trauma-triage
  standard. Measured on 150 clinical-stage therapeutics vs 150 patent-text:

  | | L1 Ready | L2 Rework | L3 Immune risk | L4 Hold | L5 Out of range |
  |---|---|---|---|---|---|
  | TheraSAbDab | 16.0% | 60.0% | 16.0% | 5.3% | **2.7%** |
  | Patent text | 12.7% | 56.7% | 13.3% | 8.0% | **9.3%** |

- **There is no GO / CONDITIONAL / NO-GO verdict any more.** Collapsing five
  levels onto three words destroyed the distinction the levels exist to carry:
  Level 2 and Level 3 both read "CONDITIONAL", but whether the fix lands in a
  CDR — and therefore whether affinity has to be re-measured — is the whole
  reason they are separate. Framework-repair candidates also read as "GO",
  which invites sending them to the bench unrepaired.
- **Immunogenicity is report-only (weight 0).** It runs, and the germline
  correction makes its flags meaningful, but calibration measured AUC 0.500 —
  exactly chance — and including it made the combined score *worse*
  (0.587 -> 0.557). Flags are reported in full; they do not move the score.
- **TAP flag regions are calibrated** (Raybould et al.) and contribute a
  small, statistically real improvement (AUC 0.578 -> 0.587, permutation
  p = 0.005). The cleanest single result in the whole calibration:
  **0 of 150 clinical-stage therapeutics carry any TAP RED flag, versus 5 of
  145 patent-text antibodies** (Fisher exact p = 0.028).

Read the flags, especially any TAP red, rather than the number. See
"What's not done" below before using this for a real go/no-go decision.

## What this is

Pre-screening for a batch of antibody VH/VL candidates,
before spending wet-lab budget on binding/neutralization assays, and collapses
the result into one verdict per candidate — **GO / CONDITIONAL / NO-GO** —
instead of a wall of per-check flags.

Two tiers:

- **Tier 1** (default): sequence-only. Cheap and fast — no GPU, no structure
  prediction, milliseconds per candidate.
- **Tier 2** (opt-in, `run_structure=True`): predicts an Fv structure with
  ABodyBuilder2 and computes the TAP developability profile. ~8s per
  candidate on a laptop CPU, cached thereafter. No GPU required.

Full detail is retained and available, just not shown by default. See
`docs/antibody-prescreen-agent-plan.md` for the Tier 1 design rationale and
`docs/tier2-tap-plabdab-plan.md` for Tier 2.

## Quick start

Runnable end-to-end example on real data, no setup beyond the [Setup](#setup)
step below — screens 10 real candidates from the committed dataset:

```bash
python3 examples/screen_candidates.py
```

For your own candidates, build the same input shape:

```python
from antibody_prescreen import screen_batch, format_report

candidates = [
    {"candidate_id": "Ab_001", "vh_sequence": "EVQL...", "vl_sequence": "DIQM..."},
    {"candidate_id": "Ab_002", "vh_sequence": "...", "vl_sequence": "..."},
]

results = screen_batch(candidates, run_immunogenicity=False)  # see note below
print(format_report(results))
```

### With Tier 2 (structure / TAP)

Tier 2 is **opt-in**. Pass `run_structure=True` and every candidate also gets
an ABodyBuilder2 model and a TAP profile:

```python
results = screen_batch(
    candidates,
    run_immunogenicity=False,
    run_structure=True,
    model_cache="~/.cache/evalab/model_cache",
)
```

Runnable example on real PLAbDab candidates:

```bash
python examples/screen_with_tap.py
```

The report grows a TAP column (`1R/2A` = one red, two amber; `all green`),
and `result.tap_profile` carries the five raw values and their flags. With
`run_structure=False` (the default) the table and the dependencies are exactly
as they were before Tier 2 existed.

`format_report` prints a markdown table, ranked GO first:

```
| Candidate | Verdict | Score | Why |
|---|---|---|---|
| Ab_001 | **GO** | 8.2 | — |
| Ab_002 | **CONDITIONAL** | 32.1 | deamidation motif 'NG' at VH CDR2 position 62 |
```

For per-candidate detail beyond the summary line, inspect `result.all_flags`
(every flag from every check, with region, severity, and weight) on the
`CandidateResult` objects `screen_batch` returns.

## What each module does

| Module | Does |
|---|---|
| `numbering.py` | IMGT numbering + CDR/framework region assignment, via ANARCI |
| `checks.py` | Sequence sanity, cysteine pairing, N-glycosylation, PTM liability (deamidation/isomerization/oxidation), V-domain integrity |
| `developability.py` | `aggregation` (AGGRESCAN, ported from ToolUniverse's antibody-engineering skill) and `charge` (pI), emitted as separate checks — "developability" was too vague a name to act on |
| `triage.py` | Assigns the triage level: the worst finding, never a sum |
| `immunogenicity.py` | Non-germline MHC-II binder scan (IEDB prediction + IEDB observed-epitope evidence) |
| `germline.py` | Per-residue germline / mutated / junctional classification, from ANARCI's own germline tables |
| `fusion.py` | Combines everything into one verdict + score + top reasons |
| `structure/modelling.py` | **Tier 2**: ABodyBuilder2 (ImmuneBuilder) VH/VL → IMGT-numbered Fv model, disk-cached |
| `structure/tap_runner.py` | **Tier 2**: TAP 5-metric profile → weighted soft flags |
| `structure/vendor/tap/` | Vendored TAP implementation (MIT, Exscientia) — see `structure/vendor/VENDOR.md` |
| `calibration/plabdab_source.py` | PLAbDab paired-sequence table as a calibration population (not used during screening) |

## Setup

One environment covers both tiers:

```bash
conda env create -f environment.yml
conda activate evalab
```

**On Apple Silicon, also install Rosetta 2 once** — TAP shells out to a
bundled x86_64 `psa` binary and there is no arm64 build:

```bash
softwareupdate --install-rosetta --agree-to-license
```

Verified working on osx-arm64 (Apple M5): `python` 3.10.20, `hmmer` 3.4,
`anarci` 2026.2.13.2, `biopython` 1.88, `openmm` 8.5.2, `pytorch` 2.10.0
(conda-forge), `ImmuneBuilder` 1.2.

Two install gotchas, both already handled in `environment.yml`:

- **Install `pytorch` from conda-forge, not pip.** The pip wheel bundles its
  own `libomp`, which collides with the one conda's `openmm` links against;
  the symptom is an abort at import time with `OMP: Error #15 ... libomp.dylib
  already initialized`.
- **`hmmscan` must be on `PATH`.** ANARCI shells out to it. Running the env's
  `python` by absolute path without activating the env is enough to break
  this, with a confusing `FileNotFoundError: 'hmmscan'` surfacing from inside
  ABodyBuilder2.

Tier 1 alone only needs `python`, `hmmer` and `anarci` — the heavy Tier 2
dependencies are imported lazily, so a Tier-1-only install still works.

## Running the tests

```bash
python -m pytest tests/ -v                          # all 39
python -m pytest tests/ -m "not slow and not network"   # 31 offline, ~2s
python -m pytest tests/ -m slow                     # 4 Tier 2 model builds, ~18s
python -m pytest tests/ -m network                  # 4 live IEDB, ~12s
```

**39 tests, all passing as of this build.**

Tier 1 (15): numbering correctness, each individual check (verified via
targeted mutation — inject a liability, assert it's caught), and end-to-end
fusion routing (clean/broken/garbage candidates, plus a regression fixture
asserting score ordering).

Immunogenicity (8, in `tests/test_immunogenicity.py`): germline assignment
and the three residue states, the FR3 core recognised as self, the grafted
epitope core recognised as non-self, V-domain anchor breakage and truncation,
plus four `network` tests running the real IEDB endpoints — the negative
control (clean antibody, 8 binders predicted, 8 dropped, 0 flags), the
positive control (grafted flu epitope survives), evidence reported
separately, and graceful degradation when the API is unreachable.

Tier 2 (10, in `tests/test_structure.py`): flag mapping and weighting, cache
keying, graceful degradation when ImmuneBuilder is missing or modelling
fails, plus four `slow` integration tests that build real models and run the
real `psa` binary. Those four follow the same mutation-injection style: the
reference pair comes back all-green, replacing CDR-H3 with poly-arginine
drives PPC from 0.59 to ~9.5 (red), and poly-aspartate drives PNC from 0.09
to ~10.9 (red) while dragging SFvCSP into amber. The `slow` tests skip
automatically if ImmuneBuilder is not installed.

## Data source

`TheraSAbDab`, `SAbDab`, `RCSB`, and `PDBe` are all unreachable from this
build environment (network egress allowlist blocks them — confirmed via
direct `curl` and the ToolUniverse `TheraSAbDab_get_therapeutic_sequences`
tool's underlying HTTP client, which hits the same block). GitHub
(`raw.githubusercontent.com`, `github.com`) is reachable.

`kevinmicha/ANTIPASTI` hosts, per PDB entry, a `lists_of_residues/{pdb}.npy`
array: every residue as an `<amino_acid><chain_letter><position>` token
(e.g. `'DH  1 '` = Asp at chain H position 1), bounded by `START-Ab`/`END-Ab`
markers separating the antibody (H+L) block from the antigen chain that
follows. `sabdab_source.py` parses these against `sabdab_summary_all.tsv`
(matching a PDB's actual npy chain letters to the correct TSV row, since a
structure with multiple Fab copies has one npy per copy, not per row) into
real VH/VL/antigen sequences.

`antibody_prescreen/data/sabdab_derived_candidates.json` is the extracted
result: 509 candidates with complete heavy, light, and antigen chains.
Regenerate it with:

```python
from antibody_prescreen.sabdab_source import load_all_candidates
candidates = load_all_candidates(
    "path/to/sabdab_summary_all.tsv",
    "path/to/ANTIPASTI/data/lists_of_residues",
)
```

**Known limitation**: sequences are read off residues actually resolved in
the crystal structure, so disordered/missing loops are silently absent from
the extracted sequence — not the same as the true, fully expressed sequence.
This can produce a spurious flag (e.g. an apparent odd cysteine count if one
member of a real pair wasn't resolved). One investigated case (`1sy6`) turned
out to be real biology (a genuine extra CDR-H3 cysteine, not a missing-residue
artifact), but this hasn't been checked for all 42 hard-gated candidates.

**Selection bias**: this population is "antibodies real enough to solve a
structure for," not "antibodies with a known real-world developability or
clinical outcome." The original Stage D ask — cross-reference threshold
placement against literature-documented liabilities — is still outstanding.

## Tier 2: structure-based TAP profiling

`run_structure=True` adds a second stage:

```
VH/VL sequence -> ABodyBuilder2 -> IMGT-numbered Fv model -> TAP -> 5 flags
```

The five TAP metrics and their flag regions (Raybould et al. 2019, PNAS
116(10):4025-4030), taken verbatim from the vendored implementation:

| Metric | GREEN | AMBER | else |
|---|---|---|---|
| Total IMGT CDR Length | 43 – 55 | 37–43, 55–63 | RED |
| Hydrophobic Patch Score (PSH) | 137.61 – 200.71 | 106.44–137.61, 200.71–225.85 | RED |
| Positive Patch Score (PPC) | 0 – 1.19 | 1.19 – 3.58 | RED |
| Negative Patch Score (PNC) | 0 – 1.67 | 1.67 – 3.50 | RED |
| SFvCSP | ≥ −4.20 | −20.50 – −4.20 | RED |

Unlike the Tier 1 thresholds, **these regions are already calibrated** — they
were set against a population of known therapeutics, which is the thing the
Tier 1 thresholds still lack.

### Why no TAP flag is a hard gate

Every TAP flag enters fusion as a *weighted soft flag*, never a hard gate.
A RED means "outside the range spanned by known therapeutics on this one
axis", and real approved antibodies do sit outside it on at least one axis.
Making it an unconditional reject would repeat the odd-cysteine hard-gate
mistake documented below. The weights live in `structure/tap_runner.py` as
named constants.

### The TAP weights are provisional

`AMBER_WEIGHT = 2.0` and `RED_WEIGHT = 6.0` are **placeholders, not fitted
values.** They were chosen only so the behaviour is predictable: against
`T_LOW = 31.0`, a single RED (6.0) is smaller than the odd-cysteine flag
(8.0) and cannot on its own push an otherwise-clean candidate out of GO,
while three REDs (18.0) moves it materially. Until Step 7 of
`docs/tier2-tap-plabdab-plan.md` runs, **the fused score with Tier 2 enabled
is not a calibrated go/no-go number.** The five TAP values and their
green/amber/red flags are meaningful on their own; the number they roll up
into is not, yet.

### Models are cached, and not bit-reproducible

Models are cached on disk keyed by a hash of the VH/VL sequence pair (not by
candidate id), so re-running a batch or re-scoring the same population does
not rebuild them. Measured on an Apple M5: **~6s to build a model, ~1.6s to
run TAP** — roughly 8s per candidate cold, ~1.6s warm.

ABodyBuilder2's refinement is **not deterministic**. Building the same
sequence pair twice produces different .pdb files and slightly different TAP
values — measured spreads of ~0.4 and ~6.3 on PSH for two different
antibodies across rebuilds. Within a project the disk cache makes results
stable, but a candidate sitting near a flag boundary can change colour if its
model is rebuilt. Do not treat a single TAP value as exact to two decimals.

### Graceful degradation

If ImmuneBuilder is missing, or modelling/TAP fails for a candidate, that
candidate still gets its Tier 1 verdict with `structure_available=False`, and
a zero-weight flag recording *why* there is no structure signal. The batch
does not die, and a silently-Tier-1 candidate is visible in `all_flags` —
the same contract `immunogenicity.py` already follows.

## Calibration data: PLAbDab

`calibration/plabdab_source.py` exists to eventually replace the SAbDab-derived
threshold calibration. **It is not used during screening**, and no calibration
has been run yet.

Two things worth recording about it:

**The documented 5 GB download is avoidable.** PLAbDab's
`plabdab_data.tar.gz` is 5,045,787,915 bytes, almost all of which is ~65k
pre-built ABodyBuilder2 models. Those tar members come *first*, so streaming
and aborting early never reaches the sequence tables either. But the files are
also served standalone from the same directory, and `paired_sequences.csv.gz`
is **11 MB** — 176,894 paired antibodies, which is all the calibration needs.
The `PLAbDab` python package, KA-Search, and the model archive are only
required for its *search* features and are not dependencies here.

**PLAbDab ships its own weak label.** The `pairing` column partitions the
table by provenance, including `TheraSAbDab` (~1,198 clinical-stage
therapeutics) and `Patent text` (~91,024 merely-patented entries, with no
developability filter applied). That gives two opposing populations for free.

Be clear about what that label is not: "reached the clinic" is not
"manufacturable", and "appears in a patent" is not "bad". It is a provenance
label, not an assay. It is a better signal than "a crystal structure exists"
— which is what `T_LOW`/`T_HIGH` are currently fitted to — and it is still
not experimental ground truth.

## Immunogenicity: the germline correction

`run_immunogenicity=True` asks whether the antibody will be seen as foreign
and destroyed by anti-drug antibodies (ADA). The chain of events:

```
APC engulfs the drug -> chops it into ~15-mers -> presents them on MHC-II
-> a CD4+ T cell recognises one -> B cells make ADA -> the drug is cleared
```

This sits next to the developability checks because it is another way the
molecule dies before it helps anyone — and the two are mechanistically
linked, not merely adjacent: aggregates are taken up by APCs far more
efficiently than monomer, so a high Tier 2 hydrophobic patch score feeds
straight into step 1.

### Why a raw MHC-II scan does not work

MHC-II binding is necessary, not sufficient. Between presentation and a
T-cell response sits **central tolerance** — T cells reactive to self peptides
are deleted in the thymus, so a germline framework peptide can be presented
beautifully and provoke nothing.

IEDB cannot supply that. Measured on the reference VH:

| | |
|---|---|
| Strong predicted binders (rank < 2.0) | 8 |
| ...overlapping a CDR | **0** |
| ...in framework | **8** |

The strongest carries `YLQMNSLRAEDT`, the IGHV3 FR3 germline motif — present
in **24.6%** of PLAbDab's 176,894 heavy chains and **34.1%** of its 1,198
clinical-stage therapeutics. Flagging it does not distinguish a risky
candidate from an approved drug.

### The filter works on the 9-mer core, not the 15-mer window

MHC-II binding is decided by a 9-residue core sitting in the groove; the
flanks hang outside it. The API reports that core in `core_peptide`.

| Filter | False positives surviving |
|---|---|
| window-level (any mutation in the 15-mer) | **8 / 8** — useless |
| **core-level** (any mutation in the 9-mer core) | **0 / 8** — correct |

All 8 shared the identical, fully germline core `YLQMNSLRA`; each window
merely happened to contain one mutated flanking residue.

**Positive control** (a filter that drops everything would pass the negative
test): grafting influenza HA306-318 into CDR-H3 gives 15 strong binders — the
8 germline framework ones are dropped, and all 7 windows of the grafted
epitope (core `YVKQNTLKL`) are kept.

Germline reference comes from ANARCI's own tables (249 human IGHV alleles,
IMGT-position-aligned), so this needs no new data source, no nucleotide
round-trip, and no web service. Three residue states, not two:

| State | Meaning | Tolerance |
|---|---|---|
| germline | matches the assigned V gene | covered |
| mutated | differs from it | **not covered** |
| beyond V | V-D-J junction, no germline exists | **not covered** |

Note that CDR3 is only *partly* junctional: IMGT 105-106 (the C-A-R stem) is
V-gene encoded; 107 onward is junction.

### Tolerance belongs to the patient, not to the molecule

The germline reference is **human by default**, regardless of what species
ANARCI assigns the chain. This distinction inverts the check if got wrong.

ANARCI gives a murine chain a mouse V gene. Comparing a murine framework
against *mouse* germline marks it "self" and drops exactly the peptides that
make murine antibodies immunogenic in people — the classic HAMA failure.
Measured on PLAbDab murine antibodies, strong binders dropped as self:

| Antibody | VH V gene | Strong binders | Dropped vs mouse ref | Dropped vs human ref |
|---|---|---|---|---|
| UCY74699_UCY74691 | IGHV5-9-3*01 | 10 | 6 | **0** |
| AVY28417_AVY28469 | IGHV2-2*03 | 11 | 9 | **0** |
| QRJ72240_QRJ72247 | IGHV1-71-16*01 | 7 | 5 | **0** |

A non-human framework now correctly reads as almost entirely non-self. Pass
`reference_species=` to model a different patient population.

### Two layers, reported separately

| Flag | Source | Claim |
|---|---|---|
| `immunogenicity` | IEDB prediction API | this peptide *might* bind MHC-II |
| `immunogenicity_observed` | IEDB IQ-API | this peptide *was observed* to provoke a T-cell response |

They are never summed into one number. The germline filter runs before
**both**: `YLQMNSLRA` has 16 positive human T-cell assay records in IEDB, so
the evidence layer false-positives on germline exactly as the predictor does.
An observed positive record is evidence, not a verdict — the study context
may be unrelated to therapeutic use.

### Network and performance

Both endpoints are live and used over **https** (the http URLs 308-redirect
and `urllib` will not re-POST — that, not a firewall, is why this check used
to report itself unavailable). Responses are cached to
`~/.cache/evalab/iedb`.

The whole chain goes in **one** POST with a comma-separated allele list: a
120-residue VH against 5 alleles returns 530 rows in ~3s. Evidence lookups
are ~0.8s each and run only on peptides that survived the filter.

Overlapping cores are clustered into one flag per liable region. Without
that, one stretch produces several shifted cores (`LLISAASSL`, `LISAASSLQ`,
`ISAASSLQS` — measured on a real candidate) and gets weighted three times.

### This check is report-only: weight 0

All three weights are **0.0**. Flags are produced with full detail and appear
in `all_flags`; they contribute nothing to the fused score.

That is a calibration result, not an oversight. Across 300 human-framework
PLAbDab antibodies (150 clinical-stage therapeutics vs 150 patent-text):

| | therapeutic | background | AUC |
|---|---|---|---|
| immunogenicity weight | median 10.00 | median 10.00 | **0.500** |
| non-germline CDR binders (mean) | 1.93 | 1.94 | 0.505 |
| observed T-cell matches (mean) | 0.03 | 0.05 | 0.510 |

And including it made the combined score worse: tier 1 + TAP scores AUC
0.587, tier 1 + TAP + immunogenicity scores 0.557.

The reading: after the germline correction what remains is predicted MHC-II
binding in CDRs, and that is about equally common in antibodies that reached
the clinic and ones that did not. Clinical antibodies are not selected
against *predicted* epitopes, and real ADA rates depend on dose, route,
duration and patient HLA — none of which a sequence scan sees.

The flags still earn their place: a specific CDR epitope with published
positive T-cell assay records is worth a human's attention. The arithmetic
does not. **Do not raise these weights without a population that shows them
predicting something.** See `docs/calibration-step7-results.md`.

## V-domain integrity (and why "productivity" is not checkable here)

`check_v_domain_integrity` verifies the conserved IMGT anchors that hold the
immunoglobulin fold together:

| IMGT position | Expected | Role |
|---|---|---|
| 23 | C | 1st-CYS, intradomain disulfide |
| 41 | W | CONSERVED-TRP, packs the hydrophobic core |
| 104 | C | 2nd-CYS, other half of the disulfide |
| 118 | F or W | J-PHE/J-TRP, marks a complete V domain |

Measured on 400 random PLAbDab pairs: **16 (4%) fail** — mostly truncation
before position 118, plus genuine substitutions (C23→S, W41→R, C104 missing).
ANARCI numbers all 400 without complaint, so it does not catch these alone.
All are soft flags, per the odd-cysteine precedent.

**This is not IMGT/V-QUEST productivity, and that check cannot be performed
on this input.** V-QUEST's definition — in-frame, no premature stop codon —
is a property of a *nucleotide* sequence. This pipeline takes amino acids,
and reverse-translating them first (as the V-QUEST tutorial workflow does)
makes the question tautological: no amino acid maps to a stop codon, and
every indel is a whole codon, so the result is in-frame and stop-free by
construction. Confirmed: **zero** of PLAbDab's 176,894 heavy chains contain
a `*`. If real sequencing-derived nucleotide data ever enters the pipeline,
V-QUEST productivity becomes meaningful again — it is not meaningful for
amino acid input.

## Scaling: what a 1,000-candidate batch costs

Measured per-unit costs on an Apple M5, and what they project to for 1,000
candidates (2,000 chains):

| Stage | per unit | 1k serial | 1k, 3 workers |
|---|---|---|---|
| ABodyBuilder2 model | 6.0 s | 100 min | **33 min** |
| TAP (psa + metrics) | 1.6 s | 27 min | 9 min |
| ANARCI numbering (×2) | 0.07 s | 1 min | <1 min |
| IEDB prediction, per chain | 2.82 s | 94 min | **94 min** (not parallelisable) |
| IEDB prediction, batched @200 | 0.31 s | 10 min | **10 min** |

**Total: ~137 min per-chain, ~53 min batched.**

### IEDB was the bottleneck. It is not any more.

The API accepts multi-FASTA and tags each result row with `seq_num`, so a
whole batch of chains goes in one request. Measured per-chain cost against
batch size (5 alleles, ~120-residue chains):

| batch | s/chain | | batch | s/chain |
|---|---|---|---|---|
| 1 | 2.82 | | 100 | 0.353 |
| 5 | 0.839 | | 200 | 0.313 |
| 10 | 0.630 | | 400 | 0.297 |
| 25 | 0.472 | | 600 | 0.304 |

~9× faster per chain, flattening around 100–200. `BATCH_SIZE = 200` is the
default: the response grows with the batch (6 MB at 200, 19 MB at 600) and a
failed request costs the whole batch — at 600 that is three minutes to redo.

`screen_batch(..., run_immunogenicity=True)` does this automatically — it
prefetches the whole batch before scoring. The prefetch writes exactly the
cache entries the per-chain path reads, so it is a pure speedup: results are
identical either way, and there is a test asserting that.

This also removes any reason to run IEDB calls concurrently, which is what
produced the order-correlated throttling that nearly invalidated the Step 7
calibration (`docs/calibration-step7-results.md`).

### So the bottleneck is now ABodyBuilder2

Options, roughly in order of effort:

1. **Turn Tier 2 off for the first pass.** Tier 1 alone screens 1,000
   candidates in about a minute. Run it first, then run TAP only on what
   survives. This is usually the right answer: the point of a cheap tier is
   to avoid paying for the expensive one.
2. **Skip immunogenicity in bulk.** It is report-only (weight 0), so it
   changes no verdict. Run it on the shortlist, where you will actually read
   the flags.
3. **Reuse the model cache.** Models are keyed by a hash of the VH/VL pair, so
   re-runs and duplicate sequences are free. In an affinity-maturation
   campaign where the light chain is constant, that is ~1,001 distinct chains
   for 1,000 candidates rather than 2,000.
4. **More cores.** Modelling is embarrassingly parallel. Three workers is the
   ceiling on a 16 GB machine — each holds torch, ABodyBuilder2 weights and
   OpenMM, and five of them got one killed mid-run. A 32–64 vCPU Linux box
   runs the same code unchanged (see Step R in
   `docs/tier2-tap-plabdab-plan.md`).
5. **Local IEDB predictor**, if you want zero network dependency at scale:
   `IEDB_MHC_II-3.0.1` is Linux-x86_64-only, so it belongs on that same
   remote box, not on a Mac. See `docs/iedb-immunogenicity-evaluation.md`.

## What's not done

- **Literature cross-referencing for calibration.** See "Data source" above —
  thresholds are now set from a real structural population's score
  distribution, which is a real step forward from an illustrative 3-candidate
  set, but still not cross-referenced against documented real-world
  developability/immunogenicity failures the way the original Stage D plan
  called for.
- ~~Cysteine hard-gate may be too blunt.~~ **Done**: downgraded from a hard
  gate to a weighted soft flag (`CYSTEINE_ODD_COUNT_WEIGHT = 8.0` in
  `checks.py`) — an odd count can be a genuine non-canonical disulfide (real
  biology, confirmed via `1sy6`) as well as a real defect, and sequence alone
  can't tell these apart. Verified against the full real population: 8 of 21
  odd-cysteine candidates now land GO/CONDITIONAL instead of an automatic
  reject. Real pairing confirmation (which Cys bonds to which) is still
  deferred to a future Tier 2 structure check — this only stops the sequence
  check from acting as an unconditional rule in the meantime.
- ~~Immunogenicity is not scored.~~ **Done**: IEDB is reachable over https,
  the check runs, and the germline correction makes it discriminate. But see
  "The weights are provisional and currently decisive" above — the signal is
  real and the weighting is not calibrated.
- **Immunogenicity uses a 5-allele HLA-DR panel.** That is a starting point,
  not population coverage, and it is a real parameter of the check rather
  than an incidental constant. No DQ/DP alleles are scanned at all.
- **The local standalone predictor is not wired in.** `IEDB_MHC_II-3.0.1`
  (515 MB) is Linux-x86_64-only, so it cannot run on Apple Silicon even under
  Rosetta 2, and it bundles NetMHCIIpan 3.2 against the API's 4.x. It is the
  right choice on the remote Linux box (Step R), not here. See
  `docs/iedb-immunogenicity-evaluation.md`.
- **Germline precedent check** (IMGT/TheraSAbDab lookup, bonus signal in the
  original design) — not built today.
- ~~Tier 2 (structure-based checks).~~ **Done**: ABodyBuilder2 → TAP profile,
  opt-in via `run_structure=True`. See "Tier 2" above — but note the weights
  feeding it into the fused score are still provisional.
- ~~Tier 2 thresholds are not calibrated.~~ **Done**: 300 human-framework
  PLAbDab antibodies scored, thresholds re-anchored, immunogenicity demoted to
  report-only. Results in `docs/calibration-step7-results.md`.
- **The calibration label is provenance, not an assay.** "Reached the clinic"
  is not "manufacturable" and "was patented" is not "bad" — most patented
  antibodies are real programs, not failures. That puts a ceiling on the
  measurable separation, and AUC 0.587 is consistent with hitting that ceiling
  rather than with the checks being uninformative. A stronger calibration
  needs measured outcomes (developability assays, manufacturing results,
  clinical ADA rates), which PLAbDab does not contain.
- **Disulfide pairing is still not checked**, even though a structure is now
  available. The cysteine note below promises this as "a future Tier 2
  structure check" — Tier 2 now exists, and this specific check still does
  not. It is its own piece of work.
- **CDR/framework boundaries use ANARCI directly**, not AbNumber's nicer
  wrapper API — AbNumber's only documented install path is conda, not
  present in this build environment. If conda becomes available, revisiting
  this is a contained change (only `numbering.py` would need to change).

## Reference sequences used during development

The VH/VL pair used throughout the test suite and in ad hoc testing during
this build (`VH_REF`/`VL_REF` in `tests/test_pipeline.py`) is **not verified
against a canonical database in this session** — RCSB and IMGT were
unreachable from this environment. Treat it as a structurally valid
human-framework test pair, not a confirmed match to any specific approved
drug, until independently verified.
