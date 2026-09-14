# Antibody Tier 1 Pre-Screening — Prototype

**Status: prototype, NOT calibrated for production use.** Thresholds
(`T_LOW=15`, `T_HIGH=40` in `fusion.py`) are provisional, set from a small
illustrative test set during today's build, not a real Stage D calibration
against approved therapeutics. See "What's not done" below before using this
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

## What's not done (explicitly out of scope for today)

- **Real Stage D calibration.** Today's thresholds come from a 3-candidate
  illustrative set built during development, not a curated set of approved
  therapeutic antibodies cross-referenced against literature. Before trusting
  this on real candidates: pull 15-30 verified VH/VL sequences from
  TheraSAbDab/IMGT, run them through the pipeline, and set `T_LOW`/`T_HIGH`
  so that set lands mostly GO/CONDITIONAL. `RCSB`, `IMGT`, and `TheraSAbDab`
  were all unreachable from this build environment's network egress
  allowlist — this has to happen from an environment that can reach them, or
  with sequences supplied directly rather than fetched.
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
