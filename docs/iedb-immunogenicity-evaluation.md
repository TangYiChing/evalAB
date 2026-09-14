# Evaluating IEDB options for the immunogenicity check

Decision record. Everything marked "measured" was actually run from this
machine (Apple M5, osx-arm64) on the date of this commit.

## TL;DR

1. **IEDB is reachable from this machine.** The README's "IEDB unreachable"
   note is stale — it described the old build sandbox. The current code's
   real blocker is that it requests `http://`, gets a 308 redirect, and
   `urllib` does not re-POST to the redirect target. Switching to `https://`
   returns 200 immediately.
2. **The current code makes ~1,060 HTTP calls per candidate to get data that
   one call returns.** The API accepts a whole sequence and a comma-separated
   allele list and does the windowing itself. Measured: full 120-residue VH ×
   5 alleles = 530 result rows in **3 seconds, one POST**.
3. **But fixing both of those still leaves the check nearly useless as
   written**, for a reason that has nothing to do with the API — see
   "The real problem" below. That finding should drive the design, not the
   API choice.

**Recommendation:** legacy `tools_api` for prediction (option A) + IQ-API as a
second, independent evidence layer (option C), with the germline correction
described at the end. Standalone (option D) is for the remote Linux box only.

---

## The real problem: the predictor mostly finds germline framework

Running the full reference VH against 5 common HLA-DR alleles gives 530
peptide/allele predictions. Of those, 8 are strong binders (percentile rank
< 2.0, the field's conventional cutoff). Measured breakdown:

| Strong binders (rank < 2.0) | Count |
|---|---|
| Overlapping a CDR | **0** |
| Framework only | **8** |

The five strongest:

```
rank=0.01  NTAYLQMNSLRAEDT  [FR3]  HLA-DRB1*04:01
rank=0.02  KNTAYLQMNSLRAED  [FR3]  HLA-DRB1*04:01
rank=0.03  NTAYLQMNSLRAEDT  [FR3]  HLA-DRB1*01:01
rank=0.05  KNTAYLQMNSLRAED  [FR3]  HLA-DRB1*01:01
rank=0.20  SKNTAYLQMNSLRAE  [FR3]  HLA-DRB1*04:01
```

That `YLQMNSLRAEDT` is the canonical IGHV3 framework-3 motif. Measured against
the PLAbDab paired table (176,894 antibodies):

| Motif | Present in |
|---|---|
| `YLQMNSLRAEDT` | 43,507 / 176,894 heavy chains (24.6%) |
| `QMNSLRAEDT` | 45,095 / 176,894 (25.5%) |
| `QMNSLRAEDT`, among **TheraSAbDab clinical-stage therapeutics** | **408 / 1,198 (34.1%)** |

So the single strongest predicted MHC-II binder in our candidate is a germline
peptide carried by a third of the antibodies that already reached the clinic.
Flagging it does not distinguish a risky candidate from an approved drug.

**Consequence for design:** a raw count of predicted MHC-II binders is not a
triage signal for antibodies. What matters is binders that are **not**
germline-encoded — i.e. in CDRs, at junctions, or in engineered/mutated
framework positions. The existing weighting (CDR 5.0 / framework 1.0) points
in the right direction but is not enough: on this candidate it would emit 8
framework flags worth 8.0 total and 0 CDR flags, on an antibody whose
framework is ordinary human germline.

This is a Goodhart problem, not an API problem. Picking a different IEDB
endpoint does not fix it.

---

## The four options

### A. Legacy tools API — `tools-cluster-interface.iedb.org/tools_api/mhcii/`

Prediction. What the current code targets.

```bash
curl --data "method=netmhciipan_el&sequence_text=$VH&\
allele=HLA-DRB1*01:01,HLA-DRB1*03:01,HLA-DRB1*04:01,HLA-DRB1*07:01,HLA-DRB1*15:01&\
length=15,15,15,15,15" \
  https://tools-cluster-interface.iedb.org/tools_api/mhcii/
```

Returns TSV: `allele seq_num start end length core_peptide peptide score rank`,
**sorted by score descending, not by position**.

- Measured: 200 OK, 530 rows, 3.0s for a full VH × 5 alleles.
- Synchronous, one request, no job polling, no auth, no SDK.

Two bugs this exposes in the current `immunogenicity.py`:

1. `IEDB_MHCII_URL` uses `http://` → 308 Permanent Redirect → `urllib` does
   not re-POST → every call returns `None` → `available=False`. This is the
   entire reason the check reports "unreachable" here.
2. `_query_iedb` posts one 15-mer for one allele at a time, and parses only
   `lines[0]`. With the whole sequence sent at once that parser would read a
   single row of a 530-row, score-sorted response — so the batching fix and
   the parser fix have to land together.

**Verdict: adopt.** Simplest thing that works, and it works today.

### B. Next-generation tools API — `nextgen-tools.iedb.org/api/v1/pipeline`

Prediction, newer method versions.

- Measured: endpoint live (`GET` → 405, so it exists and wants `POST`).
- Asynchronous pipeline: `POST /pipeline` → poll `GET /pipeline` → `GET /stage`
  → `GET /results`.

**Verdict: reject for now.** Three-to-four round trips, job state to manage,
and polling to get an answer option A returns synchronously in one call. That
complexity buys newer method versions we have no specific need for. Revisit if
the legacy endpoint is retired.

### C. IQ-API — `query-api.iedb.org`

**Experimental assay data, not prediction.** This is the API behind the
`comp-immunology/iedb-assay-gap` skill referenced in the discussion.

**First, about that skill**: it answers a different question from ours. It
finds epitopes that have T-cell assay data but lack MHC-binding assay data,
to quantify a *curation gap* in IEDB and prioritise which epitopes are worth
validating experimentally. That is a question for someone planning IEDB
assays, not for someone triaging antibody candidates. Using it as-is for
candidate screening would be a category error.

The underlying API, however, is directly useful to us for a different
question: **has any peptide in this antibody actually been observed as a
T-cell epitope in a published assay?**

```bash
curl -G --data-urlencode "linear_sequence=like.*YLQMNSLRAEDT*" \
     --data "select=structure_id,linear_sequence,qualitative_measure,source_organism_name" \
     https://query-api.iedb.org/tcell_search
```

Measured:

- 577,739 T-cell assay records total (from the `Content-Range` header with
  `Prefer: count=exact`).
- Substring queries take **0.6–0.9s** each.
- Exact-match (`eq.`) on a 15-mer returns nothing useful; substring (`like.*X*`)
  is the query that works.
- Real hit for our reference antibody's FR3 region: `SLYLQMNSLRAEDTA`,
  `qualitative_measure: Positive`, `Homo sapiens` — a published positive human
  T-cell response to an antibody framework peptide.

**Verdict: adopt, as a second layer — but only on peptides already flagged by
A.** At ~0.8s per query, scanning all 106 windows × 2 chains would cost ~170s
per candidate. Querying only the handful A flags costs ~5-15s.

The two layers have genuinely different epistemic status and should be
reported separately, never summed into one number:

| | Says |
|---|---|
| A (prediction) | "this peptide *might* bind MHC-II" |
| C (assay records) | "this peptide *was observed* to provoke a T-cell response" |

Caveat to document: a `Positive` T-cell assay in IEDB does not mean
therapeutically immunogenic. The record may come from a vaccine study, a
particular donor pool, or a non-therapeutic context. It is evidence, not a
verdict — and, per the germline finding above, a positive record on a
framework peptide that a third of marketed drugs carry is not a reason to
reject a candidate.

### D. Standalone `IEDB_MHC_II-3.0.1.tar.gz`

Local prediction, no network.

- Measured: 514,600,263 bytes (515 MB), downloadable without a licence form.
- **README: "Linux 64-bit environment", "linux 64-bit specific binaries".**
  These are Linux ELF binaries. Rosetta 2 translates x86-64 *macOS* binaries;
  it does not run Linux ELF. So this **cannot run on the MacBook** without
  Docker or a VM.
- Bundles NetMHCIIpan **3.2**; the web API serves `netmhciipan_el` 4.x.
  The local option is the *older* predictor.

**Verdict: reject for the laptop, adopt for Step R.** It is the right choice
on the remote Linux box from the Tier 2 plan: no per-candidate network calls,
fully reproducible, and version-pinnable. Not on this machine, and not at the
cost of dropping to an older method for the interactive path.

---

## Recommended shape

```
VH/VL
  |
  +-- A: one POST, full sequence x allele panel   (~3s, 1 call)
  |     -> all 15-mer windows scored
  |
  +-- germline correction                          (local, free)
  |     -> drop / downweight windows that are ordinary human germline
  |     -> keep CDR, junction, and engineered-framework binders
  |
  +-- C: IQ-API lookup on the SURVIVING peptides   (~0.8s each, ~5-15 calls)
        -> "predicted" vs "actually observed" reported separately
```

Design points to carry into implementation:

- Keep the `available=False` graceful-degradation contract. Network checks
  must never fail a batch — same as Tier 2.
- Cache IEDB responses on disk keyed by (sequence, allele panel, method).
  The batch is re-run often during calibration and the input is immutable.
- Report predicted and observed separately in `all_flags`; do not sum them
  into one immunogenicity number.
- The allele panel is a real parameter, not an incidental constant. Five
  HLA-DR alleles is a starting point, not population coverage.

## RESOLVED: how to decide "is this window germline?"

The earlier open question is now answered and validated in both directions.

### The answer: filter on the 9-mer binding CORE, not the 15-mer window

MHC-II binding is determined by a 9-residue core that sits in the groove;
the flanking residues hang outside it. The API already reports that core in
the `core_peptide` column, so no extra work is needed to find it.

**Measured on the reference VH — this distinction is the whole ballgame:**

| Filter | Strong binders surviving |
|---|---|
| window-level ("any mutation in the 15-mer → keep") | **8 / 8** — useless |
| core-level ("any mutation in the 9-mer core → keep") | **0 / 8** — correct |

All 8 false positives share the identical binding core `YLQMNSLRA`, which is
100% germline. Each 15-mer window happened to contain one mutated flanking
residue, which is why the window-level filter kept everything.

### Germline reference: ANARCI already has it

No new data source, no nucleotide round-trip, no Docker. ANARCI ships 249
human IGHV alleles (plus IGKV/IGLV) as **IMGT-position-aligned 128-character
strings** — the same coordinate system `numbering.py` already produces, so
comparison is a direct index lookup:

```python
from anarci import run_anarci, germlines
det = run_anarci([("vh", vh)], scheme="imgt", assign_germline=True)[2][0][0]
species, v_gene = det["germlines"]["v_gene"][0]      # ('human', 'IGHV3-66*01')
gl = germlines.all_germlines["V"]["H"][species][v_gene]
gl[residue.position - 1]                              # germline aa at this IMGT position
```

Three residue states fall out of this, and they are not the same thing:

| State | Meaning | Treatment |
|---|---|---|
| germline | matches the assigned V gene | self — tolerised |
| mutated | differs from the assigned V gene | engineered/somatic — **not** covered by tolerance |
| beyond V | CDR3 / J region, no V-gene germline | junctional — **not** covered by tolerance |

### Positive control: the filter keeps real signal

Grafting the textbook universal MHC-II epitope (influenza HA306-318,
`PKYVKQNTLKLAT`) into CDR-H3 and re-running:

```
15 strong binders
  8 x core YLQMNSLRA  [FR3]   germline core     -> DROPPED
  7 x core YVKQNTLKL  [CDR3]  0 germline residues -> KEPT
```

Noise out, signal in. `YVKQNTLKL` is exactly the known HA306-318 core.

### IQ-API needs the SAME filter — it is not a tolerance-aware oracle

Cross-checking both cores against IEDB's experimental T-cell assay records:

| Core | Records | Outcomes | Source organism |
|---|---|---|---|
| `YVKQNTLKL` (grafted flu epitope) | 200 | 189 positive, 11 negative | Influenza A |
| `YLQMNSLRA` (germline FR3) | 17 | **16 positive**, 1 negative | **Homo sapiens** |

The germline framework core has 16 positive human T-cell assay records. So
the experimental-evidence layer false-positives on germline exactly the way
the prediction layer does.

**Consequence:** the germline filter sits *before both layers*, not between
them. Neither IEDB endpoint knows anything about immune tolerance — that is
the part we have to supply.

## On IMGT/V-QUEST and the existing igpipeline

`IMGT_V-QUEST_Tutorial.docx` documents an existing route to the same
information: amino acid -> nucleotide (`aa2nt.py`) -> IMGT/V-QUEST via
`immcantation/suite:4.8.0` in Docker -> `output/germpass.tsv`, with
`diffcompare.py` comparing a design against germline and parent.

That route is more rigorous than the ANARCI approach — it works at nucleotide
level, resolves D and J segments, and reports V-domain productivity (in-frame,
no stop codons), which ANARCI does not address at all.

It is not usable for this check on this machine: the pipeline lives at
`/isilon/ytang4/antibody_evaluation/igpipeline` (not present locally) and
Docker is not installed here. It also answers a question adjacent to ours —
"is this sequence translatable and which genes built it" rather than "is this
residue self".

**Recommendation: use ANARCI in-process for the immunogenicity filter, and
treat V-QUEST/immcantation as the authority when the two disagree.** ANARCI
runs already, costs nothing extra, and returns the per-position germline
identity this filter needs. If `germpass.tsv` is available for a given
candidate set, prefer it — but do not make the check depend on infrastructure
that is not on this machine.

Worth noting separately: V-QUEST's productivity check (in-frame, no premature
stop) is a genuinely useful Tier 1 gate that `checks.py` does not currently
perform. That is its own piece of work, not part of this one.
