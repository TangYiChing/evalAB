# Antibody Pre-Screening Agent — Implementation Plan

**Status**: Draft, ready for hand-off
**Scope of this phase**: Sequence-only checks (Tier 1). Structure-based checks (Tier 2) are explicitly out of scope and tracked separately.
**Owner**: _unassigned_
**Last updated**: 2026-09-11

---

## 1. Motivation

Every antibody candidate that reaches a binding/neutralization assay costs real time and real money — reagents, instrument time, and a scientist's or physician's attention to interpret the result. Not all candidates deserve that spend. A meaningful fraction of any AI-generated or panel-derived candidate list contains sequences with fatal or near-fatal defects that are detectable from the sequence alone, before any wet-lab step: unpaired cysteines, CDR-embedded glycosylation sites, aggregation-prone regions, chemical degradation motifs (deamidation, oxidation, isomerization) sitting in the binding loops, and immunogenic sequence liabilities.

Today, when these checks are run at all, they're run as independent tools that each report their own flags. Handed to a scientist or physician as a stack of per-tool reports, this creates **decision fatigue**: a reviewer facing 15 yellow flags per candidate across 50 candidates either stops reading closely (and misses real problems) or spends hours triaging flags that don't actually matter (e.g., a PTM motif buried in the framework, far from the binding site, that's cheap to fix later and irrelevant to a go/no-go call right now).

The goal of this project is not to add more checks. It's to fuse the checks that already exist (or are cheap to add) into a **single, defensible verdict per candidate** — GO / CONDITIONAL / NO-GO — with the two or three reasons that actually drove the decision, so a human reviewer can act on 50 candidates in minutes instead of hours, while still being able to drill into full detail for any candidate that needs a closer look.

## 2. How this reduces decision burden

The design principle is: **run everything, show little, keep everything available.**

- All sequence-based checks run automatically and in parallel across a full candidate batch — no manual tool-by-tool invocation.
- Hard-gate defects (things that are essentially always disqualifying, e.g. an unpaired cysteine) short-circuit straight to NO-GO without needing a human to interpret a dozen soft signals.
- Everything else is fused into one risk score, weighted so that **CDR-region liabilities count for more than framework-region liabilities** — because CDR liabilities are expensive to fix (risk losing binding) and framework liabilities are cheap to fix (swap the residue, no affinity risk). This mirrors how an experienced reviewer already triages by hand, just done consistently and at scale.
- Output is one row per candidate: verdict, score, and the 1–2 reasons that mattered. Full per-check detail is available one click away, not shown by default.
- The end deliverable is a **ranked shortlist** of GO candidates plus a short "why" per CONDITIONAL candidate (so a scientist can decide whether the suggested fix is worth a redesign cycle) — not a wall of raw tool output.

## 3. Scope for this phase: Tier 1 only

This phase implements **sequence-only checks** exclusively. No structure prediction, no GPU dependency, no target sequence requirement. This is a deliberate scoping decision: Tier 1 is cheap, fast, and catches the majority of clearly-bad candidates before any structure-based (Tier 2) compute is spent on them. Tier 2 (structure-informed checks, complex prediction against a target, interface-aware liability re-scoring) is designed conceptually but **not implemented here** — see `docs/antibody-prescreen-agent-tier2-plan.md` (to be written separately) when this phase is complete.

---

## 4. Design Stages

### Stage A — Input & data model

**Goal**: define exactly what a candidate looks like going into the pipeline, so every downstream check reads the same shape.

- [ ] Define the candidate schema (suggested fields below; adjust to match what your intake actually provides):
  ```json
  {
    "candidate_id": "string, required, unique",
    "vh_sequence": "string, required, amino acid, IMGT or Kabat numbering not required at input",
    "vl_sequence": "string, required, amino acid",
    "target_sequence": "string, optional — not used in Tier 1, carried through for future Tier 2 use",
    "source": "string, optional — e.g. which design tool/campaign generated this candidate",
    "notes": "string, optional, free text"
  }
  ```
- [ ] Decide and document input format: CSV/TSV with one row per candidate, or FASTA pairs with a naming convention (`<id>_VH`, `<id>_VL`), or JSON lines. Pick one as primary; a converter to/from the others is a nice-to-have, not a blocker.
- [ ] Write input validation: reject (with a clear per-row error, not a silent skip) any row with a non-standard amino acid character, empty sequence, or missing required field, before it reaches the checks below.
- [ ] Decide how CDR boundaries are determined for every downstream check that needs CDR-vs-framework distinction. Recommended: IMGT numbering scheme (unambiguous CDR/FR boundary definitions) via a standard numbering library (e.g. ANARCI) rather than hand-rolled regex guessing — get this right once, since every weighting decision downstream depends on it.

### Stage B — Individual sequence-based checks

Each check should be its own small, independently testable function/script that takes a sequence (or VH/VL pair) and returns a structured result — not free text. Implement in this order (each is independently useful, so partial progress still ships value):

- [ ] **Cysteine pairing check** (hard gate). Count cysteines in VH and VL; flag odd counts or any cysteine outside the expected conserved positions as a hard gate.
- [ ] **Sequence sanity check** (hard gate). Reject non-standard amino acids, embedded stop codons/frameshift artifacts, and CDR-H3 length outside a configurable normal range (default suggestion: 4–30 residues; confirm against your own design tool's typical output range before finalizing).
- [ ] **N-glycosylation motif scan** (hard gate if in CDR, soft flag if in framework). Regex for N-x-[S/T] where x ≠ P. Requires CDR/framework boundary from Stage A.
- [ ] **PTM liability scan** (soft flag, CDR-weighted). Regex-based motif scan: deamidation (N-G, N-S, N-T, N-H), isomerization (D-G, D-S), oxidation-prone residues in exposed positions (M, W, unpaired C — note overlap with the cysteine check, don't double-penalize the same residue). This does not exist yet in the ToolUniverse antibody-engineering skill (confirmed gap) — build from scratch here; it's straightforward regex, no external dependency.
- [ ] **Aggregation propensity / developability**. Reuse the existing AGGRESCAN-style implementation (`scripts/developability.py` from the ToolUniverse `tooluniverse-antibody-engineering` skill) rather than rewriting — it already computes aggregation propensity, pI, and hydrophobic patches in pure Python with no external dependency. Wrap its output into this pipeline's structured result format.
- [ ] **Immunogenicity scan (own-sequence MHC-II binding)** — corrected approach, not the skill's target-based IEDB search. Slide a 15-mer window across VH and VL, score each window against a panel of common human HLA-DR alleles using IEDB's Analysis Resource MHC-II binding prediction tool (NetMHCIIpan or equivalent), independent of target sequence. Flag high-affinity predicted binders, weighted higher if inside a CDR.
- [ ] **Germline precedent check** (bonus signal, not a penalty). Query IMGT/TheraSAbDab for whether the VH/VL framework matches a well-precedented clinical-stage germline. Treat as a positive signal in the score, not a required gate — novel frameworks aren't automatically bad, just less de-risked by precedent.

For each check, define and document: (a) what it returns (numeric score, boolean, or categorical), (b) whether it's a hard gate or soft flag, (c) its default weight if soft, (d) which region (CDR vs framework) it applies to when relevant.

### Stage C — Score fusion & verdict logic

- [ ] Implement the fusion function: hard gates short-circuit to NO-GO with the triggering reason; otherwise compute a weighted risk score from all soft flags, with CDR-region flags weighted 3–5x framework-region flags (exact multiplier to be set in Stage D calibration, not hardcoded here).
- [ ] Implement the three-bucket classification (GO / CONDITIONAL / NO-GO) against two thresholds (`T_low`, `T_high`) — thresholds are placeholders until Stage D calibration sets real values.
- [ ] For CONDITIONAL and NO-GO verdicts, implement "top-2-reasons" extraction: sort triggered flags by weighted contribution, surface only the top 1–2 in the summary output, keep the rest accessible in a detail view/expanded record.
- [ ] For CONDITIONAL verdicts specifically, generate a one-line suggested fix where mechanically derivable (e.g., "framework PTM at position 48 (N-G) — consider N48Q substitution"); leave as "manual review" where the fix isn't a simple substitution (e.g., CDR-embedded liabilities).

### Stage D — Calibration on approved mAbs

**Goal**: set real, defensible thresholds instead of arbitrary numbers, by running the pipeline against antibodies whose real-world outcome is already known.

- [ ] Assemble a calibration set of 15–30 approved therapeutic antibody sequences (VH/VL), pulled from TheraSAbDab or IMGT. Include a spread: some early first-generation mAbs (chimeric/humanized, known to have had manufacturability or immunogenicity issues historically) alongside recent, highly optimized clinical-stage antibodies. The spread matters — you want to see the score distribution across both "acceptable but imperfect" and "excellent" cases, not just the best-in-class ones.
- [ ] Run the full Tier 1 pipeline (Stages B + C, with placeholder thresholds) against every sequence in the calibration set.
- [ ] Record each calibration antibody's raw score and which flags contributed, in a shared spreadsheet/table — this becomes the evidence base for threshold decisions, not just a one-time computation.
- [ ] Cross-reference against any known real-world manufacturability/immunogenicity issues for antibodies in the set (literature search — e.g., published ADA rates, known reformulation efforts) where available, to sanity-check that the score direction makes sense (antibodies with known real issues should score worse than clean ones).
- [ ] Set `T_low` and `T_high` such that the bulk of the approved-mAb calibration set lands in GO/CONDITIONAL, not NO-GO. If more than a small handful of approved, marketed antibodies land in NO-GO, the thresholds (or a check's weight) are miscalibrated — revisit before finalizing.
- [ ] Document the final thresholds and the reasoning behind them (not just the numbers) in this file's Appendix, so a future reviewer understands why `T_low`/`T_high` are what they are, not just what they are.
- [ ] Re-run the full calibration set after any weight/threshold change and confirm the distribution still makes sense — calibration is not a one-shot step, treat it as a regression check.

### Stage E — Testing

- [ ] Unit tests per individual check (Stage B): known-bad sequence in → expected flag out, for each check independently. Include edge cases: sequence with zero cysteines, sequence with a CDR-H3 at each length boundary, sequence with multiple overlapping liability motifs in the same window.
- [ ] Unit tests for CDR/framework boundary assignment (Stage A) against a few manually-verified reference sequences — this is the check most likely to silently produce wrong weighting if it's wrong, so verify it directly rather than only through downstream score tests.
- [ ] Integration test: full pipeline run on a small fixture set (5–10 sequences) with known expected verdicts, checked into the repo as a regression fixture.
- [ ] Regression test against the Stage D calibration set: after any pipeline change, re-run the calibration set and confirm no unexpected verdict flips (a script that diffs current run against the last-recorded calibration run, not just eyeballing).
- [ ] Manual review pass: have a domain scientist review a sample of CONDITIONAL verdicts and confirm the suggested fixes are sensible, before treating the pipeline's output as trustworthy for real candidate triage.

### Stage F — Wrap-up

- [ ] Write a short usage README: how to run the pipeline on a candidate batch, what the output format looks like, how to interpret GO/CONDITIONAL/NO-GO.
- [ ] Document known limitations explicitly, so users don't mistake this for a full developability/immunogenicity assessment: this is a sequence-only pre-screen; it does not predict binding affinity, does not assess the target-binding interface, and does not replace wet-lab confirmation. It exists to cut obviously bad candidates before spending assay budget, not to pick the final winner.
- [ ] Hand off open items for Tier 2 (structure-based checks) as a separate, explicitly out-of-scope follow-up — do not let Tier 2 scope creep into this phase.
- [ ] Confirm with the receiving team member: data input format, where the calibration set and its results live, and who owns threshold updates going forward.

---

## Appendix: Open decisions for the incoming owner

- [ ] Final CDR numbering scheme choice (IMGT recommended, confirm no conflicting internal convention already in use).
- [ ] Final `T_low` / `T_high` threshold values and the calibration evidence behind them (fill in after Stage D).
- [ ] Final CDR-vs-framework weight multiplier (fill in after Stage D).
- [ ] Primary input file format decision (CSV/FASTA/JSONL).
- [ ] Whether the germline precedent check should ever act as a hard requirement rather than a bonus signal (recommendation in this plan: no, keep it a bonus — novel frameworks shouldn't be auto-rejected).
