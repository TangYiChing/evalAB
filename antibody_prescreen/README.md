# Antibody Tier 1 Pre-Screening — Prototype

**Status: prototype, calibrated against real structural data but NOT against
known clinical/manufacturing outcomes.** Thresholds (`T_LOW=31`, `T_HIGH=40`
in `fusion.py`) are set from the score distribution of 508 real,
structurally-solved antibodies (see "Data source" below) — a real
improvement over an earlier 3-candidate placeholder pass, but this
population is selected for "a structure was solved," not "known good or bad
in the clinic/manufacturing." See "What's not done" below before using this
on real candidates for a real go/no-go decision.

## What this is

Sequence-only (Tier 1) pre-screening for a batch of antibody VH/VL candidates,
before spending wet-lab budget on binding/neutralization assays. Runs cheap,
fast checks (no GPU, no structure prediction) and collapses the result into
one verdict per candidate — **GO / CONDITIONAL / NO-GO** — instead of a wall
of per-check flags. Full detail is retained and available, just not shown by
default. See `docs/antibody-prescreen-agent-plan.md` for the design rationale.

## Quick start

```python
from antibody_prescreen import screen_batch, format_report

candidates = [
    {"candidate_id": "Ab_001", "vh_sequence": "EVQL...", "vl_sequence": "DIQM..."},
    {"candidate_id": "Ab_002", "vh_sequence": "...", "vl_sequence": "..."},
]

results = screen_batch(candidates, run_immunogenicity=False)  # see note below
print(format_report(results))
```

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
| `checks.py` | Sequence sanity, cysteine pairing, N-glycosylation, PTM liability (deamidation/isomerization/oxidation) |
| `developability.py` | Aggregation propensity (AGGRESCAN, ported from ToolUniverse's antibody-engineering skill) + pI |
| `immunogenicity.py` | Own-sequence MHC-II binding scan (IEDB API) — **not scored in this environment**, see below |
| `fusion.py` | Combines everything into one verdict + score + top reasons |

## Setup

```bash
apt-get install -y hmmer
pip install anarci pytest
```

Verified working versions: `hmmer` 3.4+dfsg-2, `anarci` 2026.2.13.2.

## Running the tests

```bash
python3 -m pytest tests/ -v
```

15 tests, all passing as of this build: numbering correctness, each
individual check (verified via targeted mutation — inject a liability, assert
it's caught), and end-to-end fusion routing (clean/broken/garbage candidates,
plus a regression fixture asserting score ordering).

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

## What's not done (explicitly out of scope for today)

- **Literature cross-referencing for calibration.** See "Data source" above —
  thresholds are now set from a real structural population's score
  distribution, which is a real step forward from an illustrative 3-candidate
  set, but still not cross-referenced against documented real-world
  developability/immunogenicity failures the way the original Stage D plan
  called for.
- **Cysteine hard-gate may be too blunt.** An odd cysteine count can be a
  genuine non-canonical disulfide (real biology, confirmed in one case) as
  well as a real defect or a missing-residue artifact — sequence alone can't
  tell these apart. Worth downgrading to a strong soft flag requiring human
  review rather than an automatic hard gate, once there's a way to
  distinguish the cases (e.g. structure-based confirmation in Tier 2).
- **Immunogenicity is not scored.** IEDB's MHC-II API host
  (`tools-cluster-interface.iedb.org`) is not on this environment's network
  allowlist. The check degrades gracefully (returns `available=False`, no
  crash) rather than silently skipping — but no candidate in this build has
  actually been scored for immunogenicity. Allowlist that host, or wire in a
  local NetMHCIIpan install, before relying on this signal.
- **Germline precedent check** (IMGT/TheraSAbDab lookup, bonus signal in the
  original design) — not built today.
- **Tier 2 (structure-based checks)** — explicitly out of scope per the
  original plan, tracked separately.
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
