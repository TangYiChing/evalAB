# Triage implementation plan: from weighted sum to likelihood-ratio gating

**Status: design document. Nothing here is implemented yet.**
Every number quoted was measured on this repository's own data; the provenance
of each is stated inline. Where a check does not exist, it says so — there are
no placeholders standing in for unimplemented work anywhere in this plan, and
none may be added. A check that is not implemented must return "not available",
never a fabricated value.

**If you are picking this up: read §1 to understand why, §5 to see what exists,
then start at Phase 0 in §6.** Everything else is supporting material.

---

## 1. Rationale: what is wrong with the current approach

### 1.1 What we built, and what it actually does

`fusion.py` sums weighted soft flags into one score and cuts it at two
thresholds (`T_LOW = 29.41`, `T_HIGH = 38.31`) to produce GO / CONDITIONAL /
NO-GO. Those thresholds are percentile placements — median and p90 — on 150
clinical-stage therapeutics (`docs/calibration-step7-results.md`).

Three measurements say this cannot support the decision it is being asked to
support.

**(a) There is no gap between GO and NO-GO.** The score distribution on known-
good antibodies is a single smooth unimodal hill:

```
 15-19 ##
 20-24 ####################
 25-29 ####################################################
 30-34 ##########################################
 35-39 #####################
 40-44 ##########
 45-49 ###
```

The thresholds are two lines drawn across the middle of one population. 17% of
known-good antibodies sit within ±3 points of `T_HIGH` — close enough that one
extra deamidation motif moves a candidate across it.

**(b) The false-rejection rate is defined, not measured.** `T_HIGH` was set at
the p90 of the therapeutic population, so "10% of known-good antibodies are
NO-GO" is a restatement of where the line was drawn, not a finding. Recomputed
with current scoring (immunogenicity weight 0), the actual figure is **12.7%**
of clinical-stage therapeutics rejected.

**(c) Two of the checks carry no information at all.** Likelihood ratios
measured across 150 therapeutics vs 150 patent-text antibodies:

| check | P(flag \| good) | P(flag \| bad) | LR⁺ |
|---|---|---|---|
| `sequence_sanity` | 0.0% | 3.3% | 11.00 |
| TAP any RED | 0.0% | 3.3% | 11.00 |
| `n_glycosylation` | 3.3% | 9.3% | 2.64 |
| `cysteine_pairing` | 5.3% | 11.3% | 2.06 |
| `ptm_liability` | 100.0% | 100.0% | **1.00** |
| `developability` (AGGRESCAN) | 100.0% | 100.0% | **1.00** |
| `v_domain_integrity` | 0.7% | 0.0% | 0.33 |

(Haldane–Anscombe 0.5 correction applied so a zero cell does not produce an
infinite ratio.)

An LR of exactly 1.00 means observing the flag should not change your belief at
all. `developability` and `ptm_liability` fire on **100% of clinical-stage
therapeutics and 100% of patent-text antibodies**. They are contributing weight
to a score while contributing zero information to the decision.

### 1.2 The distinction that matters: ranking vs triage

A continuous aggregation score is not useless — it is being *read* the wrong
way. Measured distributions of the AGGRESCAN-derived score:

| population | min | median | max |
|---|---|---|---|
| therapeutic (n=150) | 8.58 | 17.58 | 24.56 |
| patent text (n=150) | 10.65 | 18.05 | 24.53 |

Nearly complete overlap. There is no cut point anywhere in that range that a
bad antibody crosses and a good one does not — a gate at 24.6 would reject 0
known-good and catch 0 background.

But the same underlying physics read as a **discrete boundary crossing** does
work: TAP's hydrophobic-patch RED flag is part of the `TAP any RED` gate above,
which has LR 11.

So the finding is not "aggregation is uninformative". It is:

> **The same physical quantity carries no information when used to rank a
> population, and real information when used to ask whether a candidate has
> left the region occupied by every known therapeutic.**

That sentence is the whole rationale for this redesign.

### 1.3 Why weighted sums are the wrong shape for this decision

A weighted sum assumes the weights are commensurable — that 2 points of
aggregation trades off against 2 points of PTM liability in some meaningful
way. Nothing in our data supports that assumption, and the weights currently in
the codebase were chosen by hand.

Worse, a sum lets many weak signals overwhelm one decisive one. A candidate
whose immunoglobulin fold is broken (2nd-CYS absent) and a candidate with six
framework deamidation motifs can land on the same score. They are not the same
kind of problem and should not produce the same output.

---

## 2. Literature: the two ideas we are borrowing

### 2.1 Likelihood ratio — the principled replacement for hand-set weights

The diagnostic-testing literature has a standard answer to "how much should
this finding move my belief":

$$\text{LR}^+ = \frac{P(\text{flag} \mid \text{bad})}{P(\text{flag} \mid \text{good})}
\qquad
\text{LR}^- = \frac{P(\text{no flag} \mid \text{bad})}{P(\text{no flag} \mid \text{good})}$$

Two properties make this the right tool here:

1. **It is estimated from data, not chosen.** Each check's weight becomes an
   empirical quantity with a confidence interval, not a constant someone typed.
2. **It composes correctly.** Log-LRs add: posterior log-odds = prior log-odds
   + Σ log LR. So if you must have a score, `log₂ LR` *is* the principled
   weight — and applying that rule automatically sets `developability` and
   `ptm_liability` to weight 0.00, because that is what their data says.

The caveat, which must be stated wherever this is used: **log-LRs only add
cleanly when the findings are conditionally independent.** Our checks are not —
CDR deamidation motifs and CDR length are obviously correlated. Section 6
Phase 3 addresses this.

**Clinical decision rules** are the applied form of this. The Ottawa Ankle
Rules report ~98% pooled sensitivity and >99% NPV, and — this is the part we
have never done — they were **derived on one population and then externally
validated on independent ones**. A rule that has only been derived is a
hypothesis. Ours is currently a hypothesis.

### 2.2 Emergency Severity Index — triage is not scoring

ESI, the five-level triage algorithm used in emergency departments, has **no
weights anywhere in it**. It is a decision tree asked in a fixed order:

1. Does this patient need an immediate life-saving intervention? → Level 1
2. High-risk situation, confused/lethargic, or severe pain? → Level 2
3. Otherwise: **how many resources will this patient consume?**
   → many = Level 3, one = Level 4, none = Level 5

Two ideas transfer, and the second one is the interesting one.

**Idea 1: a level is defined by the action it triggers, not by a number.**
Level 1 does not mean "score above X". It means "goes to resuscitation now".
The antibody analogue: a triage level should name what happens to the candidate
— goes to the bench, goes back for one round of engineering, is not pursued —
not report a number the user then has to interpret.

**Idea 2: the lower levels triage by predicted resource consumption, not by
severity.** ESI counts concrete resources: labs, imaging, ECG, IV fluids, IV/IM
medications, specialist consults.

Translated to antibody triage, the question stops being "is this antibody
good?" and becomes:

> **"How many engineering cycles will this candidate consume before it is worth
> putting on the bench?"**

- One framework N-glyc motif → one round of site-directed mutagenesis → **1 resource**
- Three CDR deamidation motifs plus charge asymmetry → CDR redesign → **many resources**
- 2nd-CYS absent, fold not supported → **unbounded** — not "bad", but *not on this path*

This reframing dissolves the weighting problem, because counting resources
counts commensurable things. It also produces output a bench scientist can act
on directly.

### 2.3 Over-triage vs under-triage — the design target we lacked

The American College of Surgeons Committee on Trauma sets an explicit,
quantified, asymmetric standard: **under-triage below 5%, with over-triage of
25–35% accepted as the price**. Missing a severely injured patient is worse
than mobilising a trauma team unnecessarily, and the numbers say by how much.

Applying that standard to what we have:

| mechanism | false rejection of known-good | meets <5%? |
|---|---|---|
| score-based NO-GO (`T_HIGH = 38.31`) | **12.7%** | no |
| TAP any RED gate | **0.0%** (0/150) | yes |
| `sequence_sanity` gate | **0.0%** (0/150) | yes |

This gives us the acceptance criterion the current design never had: **a gate
may only enter Level 1 if its measured false-rejection rate on known-good
antibodies is below 5%, with a confidence interval, on a population it was not
fitted to.**

### 2.4 Alert fatigue — "click fatigue" has a literature and a fix

Clinicians override 49–96% of drug-interaction alerts. The evidence-based
remedies are not better-looking alerts; they are **tiering by risk**, **moving
low-value alerts from interruptive pop-ups to passive in-workflow display** (one
study reduced alert volume 80% this way), and **role-tailoring**.

The mechanism matters for us: an alert with LR = 1 that fires on every
candidate does not merely waste attention — it *trains the user to ignore the
LR = 11 alert sitting next to it*. Our `developability` and `ptm_liability`
checks fire on 100% of candidates. They are the click-fatigue source, and
demoting them is a correctness fix, not a cosmetic one.

### 2.5 Prior art, and where this sits

Multi-parameter optimisation, desirability functions (Derringer–Suich), and
Pareto / non-dominated sorting are all established in antibody and small-
molecule candidate selection. Every one of them is a framework for **ranking**.
Pareto approaches are notable for needing no user-specified weights at all, and
are a reasonable fallback if the LR route stalls.

What we did not find in the literature: treating **triage itself** as the
object of study for antibody design — deriving gates, validating them
externally, and holding them to an explicit under-triage budget the way trauma
systems do. That is the contribution this plan is aiming at.

### 2.6 A published design pipeline already uses this shape

AbSci's Origin-1 supplementary material (Supplementary Table 7,
`260114_Origin1_Final_Submission.pdf`) lists the sequence liabilities used
while designing their libraries. The structure is worth reading closely,
because it is not a weighted score:

| Liability | Definition (verbatim) | Category |
|---|---|---|
| N-Glycosylation | the motif NXS or NXT with X any residue other than proline | **Critical** |
| 3-W | Five consecutive tryptophan residues in a single CDR | **Critical** |
| 4-G | Four consecutive glycine residues in a single CDR | **Critical** |
| 5-Y | Five consecutive tyrosine residues in a single CDR | **Critical** |
| 5-S | Five consecutive serine residues in a single CDR | **Critical** |
| 3-G | Three consecutive glycine residues in a single CDR | Other (4 max) |
| 3-S | Three consecutive serine residues in a single CDR | Other (13 max) |
| 4-S | Four consecutive serine residues in a single CDR | Other (3 max) |

Two things to notice.

**It is two-tier, not summed.** "Critical" liabilities are filtered out
outright during Sequence Search; "Other" liabilities are permitted under a
budget. That is the Level-1/Level-3 split in §3, arrived at independently by a
group actually shipping designs.

**"Other" carries explicit quotas.** The paper describes selecting the top 95
structure/sequence pairs "restricting the number of 'Other' sequence
liabilities permissible". Read alongside "(13 max)" against a 95-design
library, this reads as a **per-library quota** rather than a per-sequence
limit — i.e. resource allocation across a batch, which is exactly the ESI
Level-3/4/5 idea. **This reading is an inference, not something the text states
outright; confirm it before implementing quota semantics.**

**Discrepancy to resolve before implementing**: the row labelled `3-W` is
defined as *five* consecutive tryptophans. The label and the definition
disagree. Implement to the written definition (five), record the ambiguity in
the code comment, and do not silently pick one.

---

## 3. Target architecture

Three levels, each defined by the action it triggers.

**IMPLEMENTED. Numbering follows ESI: 1 is the one you act on first, 5 the one
you do not pursue.** An earlier draft of this document had it backwards, with
Level 1 meaning "reject" — the opposite of what Level 1 means to anyone who has
seen ESI, where Level 1 is the resuscitation patient.

| Level | Meaning | Action |
|---|---|---|
| 1 | no repairs needed | straight to the bench |
| 2 | framework repairs only | one round, no re-measurement |
| 3 | repairs land in a CDR | re-measure affinity afterwards |
| 4 | a human must adjudicate | serious, not disqualifying |
| 5 | crossed a boundary | do not pursue |

**Level eligibility is decided by P(flag | known-good), not by LR magnitude.**
An earlier version of this section said Level 2 was "checks with 2 ≲ LR ≲ 11",
which was wrong twice over: the upper bound was just the largest value we
happened to observe, and more importantly LR magnitude is the wrong axis. What
makes an automatic rejection safe is not firing on good candidates. A check
with LR 3 and P(flag|good) = 20% cannot gate; one with LR 3 and
P(flag|good) = 0.5% can. LR tells you what a flag is worth once it fires.

**Level 1 — hard exclude.** Entry requirement: measured
P(flag | known-good) < 5% with a confidence interval, on a population the gate
was not fitted to. Output names the specific boundary crossed. Interruptive.

**Level 2 — flagged for human review.** Checks with real but non-decisive LR.
Displayed passively next to the candidate. Never blocks. This is where the
alert-fatigue remedy lives.

**Level 3/4/5 — resource count.** Replaces the score. Counts *fixable*
liabilities and where they sit:

| Level | Meaning | Typical content |
|---|---|---|
| 5 | 0 fixable liabilities | straight to bench |
| 4 | 1 framework-region liability | one round of point mutagenesis |
| 3 | multiple, or any in CDR | CDR redesign; queue behind 4 and 5 |

`developability` and `ptm_liability` move here. "This sequence has 3 CDR
deamidation motifs" is a useful workload estimate. "Its aggregation score is
18.4" is not a useful ranking.

**The fused score is removed from the decision path entirely.** It may remain
as a diagnostic value in `CandidateResult`; nothing may route on it.

---

## 4. Data inventory

Every dataset, when it is used, and what for. **A dataset used for deriving a
gate may never also be used to validate it.**

| # | Dataset | Size | Where it lives | Used at | Used for |
|---|---|---|---|---|---|
| D1 | PLAbDab `paired_sequences.csv.gz` | 11 MB, 176,894 paired | `~/.cache/evalab/plabdab/`, fetched by `calibration/plabdab_source.py` | Phases 1–3 | source of everything below |
| D2 | **Derivation set** — 150 TheraSAbDab + 150 patent-text, human framework, seed 0 | 300 | `docs/calibration-step7-results.csv` | Phase 1 | estimating LR per check. **Already used. Burned for validation.** |
| D3 | **Validation set** — TheraSAbDab entries NOT in D2 | ~1,048 | to be built, seed ≠ 0 | Phase 2 | measuring false-rejection rate of Level 1 gates. **Must not be looked at until gates are frozen.** |
| D4 | **Background validation set** — patent-text not in D2 | ~sample 1,000 | to be built | Phase 2 | measuring what fraction of background the gates catch (over-triage side) |
| D5 | Non-human-framework subset | ~25% of D1 | derivable via `--include-non-human` | Phase 4 | separate derivation for murine/chimeric input; see §7.2 |
| D6 | **de novo design set** | none yet | — | Phase 5 | **does not exist.** See §7.1 — gates derived on D2/D3 are not validated for de novo input |
| D7 | ABodyBuilder2 model cache | ~280 KB/model | `~/.cache/evalab/model_cache/` | Phases 1–3 | TAP input; content-addressed, safe to share across runs |
| D8 | IEDB response cache | varies | `~/.cache/evalab/iedb/` | Phases 1–3 | immunogenicity; report-only, not on the decision path |

**The D2/D3 split is the single most important line in this table.** Every LR
in §1.1 was estimated on D2. Re-measuring those same gates on D2 and calling it
validation would be circular. D3 exists to be the independent test, and its
value is destroyed the moment someone looks at it before the gates are frozen.

---

## 5. Liability check inventory

Status vocabulary, used strictly:

- **IMPLEMENTED** — code exists, tested, LR measured on D2
- **IMPLEMENTED, LR UNMEASURED** — code exists, no LR yet
- **NOT IMPLEMENTED** — no code. Must return "not available" if queried.
  **Must not be faked with a constant, a random draw, or a stand-in metric.**

### 5.1 Implemented

| Check | Module | LR⁺ on D2 | P(flag\|good) | Proposed level |
|---|---|---|---|---|
| `sequence_sanity` | `checks.py` | 11.00 | 0.0% | **Level 1 candidate** |
| TAP any RED | `structure/tap_runner.py` | 11.00 | 0.0% | **Level 1 candidate** |
| `n_glycosylation` | `checks.py` | 2.64 | 3.3% | Level 2 |
| `cysteine_pairing` | `checks.py` | 2.06 | 5.3% | Level 2 |
| `v_domain_integrity` | `checks.py` | 0.33 (n too small) | 0.7% | Level 1 candidate — LR needs re-estimation on D3 |
| `ptm_liability` | `checks.py` | 1.00 | 100.0% | **Level 3/4/5 resource count** |
| `developability` (AGGRESCAN) | `developability.py` | 1.00 | 100.0% | **Level 3/4/5 resource count** |
| `immunogenicity` (MHC-II + germline) | `immunogenicity.py` | 0.500 AUC | — | report-only, weight already 0 |

Note on `v_domain_integrity`: its D2 LR of 0.33 comes from 1 therapeutic vs 0
background tripping it. That is noise, not an inversion. Its 0.7% flag rate on
known-good is what qualifies it for Level 1; the LR must be re-estimated on the
larger D3.

### 5.2 AbSci Origin-1 liabilities — NOT IMPLEMENTED

None of the poly-residue-run checks below exist in the codebase. Implementing
them is Phase 1 work.

| Check | Definition | AbSci category | Status |
|---|---|---|---|
| poly-W run in a CDR | 5 consecutive W (see label discrepancy, §2.6) | Critical | **NOT IMPLEMENTED** |
| poly-G run in a CDR | 4 consecutive G | Critical | **NOT IMPLEMENTED** |
| poly-Y run in a CDR | 5 consecutive Y | Critical | **NOT IMPLEMENTED** |
| poly-S run in a CDR (5) | 5 consecutive S | Critical | **NOT IMPLEMENTED** |
| poly-G run in a CDR (3) | 3 consecutive G | Other (4 max) | **NOT IMPLEMENTED** |
| poly-S run in a CDR (3) | 3 consecutive S | Other (13 max) | **NOT IMPLEMENTED** |
| poly-S run in a CDR (4) | 4 consecutive S | Other (3 max) | **NOT IMPLEMENTED** |

N-glycosylation is the one Origin-1 Critical liability we already have
(`check_n_glycosylation`), and our definition matches theirs (NXS/NXT, X ≠ P).

**Measured base rates of these motifs on D2, so whoever implements them knows
what to expect** (computed by scanning CDR region sequences directly; no check
code was written for this):

| Motif | P(\|good) | P(\|bad) | LR⁺ |
|---|---|---|---|
| 5×W | 0.0% | 0.0% | 1.00 |
| 4×G | 1.3% | 0.0% | 0.20 |
| 5×Y | 0.0% | 1.3% | 5.00 |
| 5×S | 0.7% | 1.3% | 1.67 |
| 3×G | 4.0% | 4.0% | 1.00 |
| 3×S | 14.0% | 19.3% | 1.37 |
| 4×S | 0.7% | 4.0% | 4.33 |

**Read this table carefully before doing anything with it.** The Critical
motifs essentially never occur in natural antibodies — 5×W occurs zero times in
all 300. Zero-versus-zero gives no usable LR at any sample size we can reach
from PLAbDab.

That is not a defect in the checks. It tells you what they are *for*: these are
**generative-model failure modes**, not natural-antibody failure modes. A
protein language model or diffusion model emitting `SSSSS` in CDR-H3 is a
sampling artefact that evolution does not produce. They cannot be calibrated on
natural antibodies, and their LR must be derived on de novo output (D6, §7.1).

So they enter the system as **Level 1 gates justified by mechanism and by
precedent, with their LR explicitly marked as unmeasured** — not as calibrated
gates. That distinction must survive into the code comments and the report
output.

### 5.3 Known gaps — NOT IMPLEMENTED, nothing pretending to stand in

| Gap | Why it matters | Status |
|---|---|---|
| pTM / per-residue confidence from ABodyBuilder2 | the model's own uncertainty is a strong "this does not look like an antibody" signal, especially for de novo | **NOT IMPLEMENTED** — ABodyBuilder2 produces per-residue error estimates; we currently save only the PDB and discard them |
| Humanness score (e.g. OASis-style) | AbSci used BioPhi in lead optimisation; we have no humanness measure | **NOT IMPLEMENTED** |
| Disulfide pairing from structure | promised in the README as "a future Tier 2 structure check"; Tier 2 now exists, this does not | **NOT IMPLEMENTED** |
| V-QUEST nucleotide productivity | not assessable from amino-acid input — see `README.md`. `check_v_domain_integrity` is the AA-level analogue and is **not** a substitute | **N/A by construction, documented** |
| Viscosity / self-association | not covered by any current check | **NOT IMPLEMENTED** |
| External validation of any gate | every LR above is derivation-only | **NOT DONE** — this is Phase 2 |

---

## 6. Implementation phases

### Phase 0 — freeze and instrument (start here)

No behaviour change. Make the current state measurable.

1. Add `calibration/lr_report.py`: given a scored population CSV with a
   `pairing` label, emit per-check P(flag|good), P(flag|bad), LR⁺, LR⁻, and
   Wilson confidence intervals. The numbers in §1.1 were computed ad hoc; this
   makes them reproducible and re-runnable.
2. Emit per-check boolean columns from `run_calibration.py` (currently only
   aggregate weights are written), so LR can be recomputed without re-scoring.
3. **Do not touch D3.**

Done when: `lr_report.py` reproduces the §1.1 table from
`docs/calibration-step7-results.csv`.

### Phase 1 — implement the missing checks

1. `check_poly_residue_runs` in `checks.py`, covering the Origin-1 motifs in
   §5.2. CDR-only, per-CDR (a run split across CDR1 and CDR2 does not count).
   Resolve the 3-W label discrepancy per §2.6 and comment it.
2. Extract ABodyBuilder2's per-residue error into `structure/modelling.py` and
   surface it. This is the single highest-value missing signal for de novo
   input.
3. Every new check returns flags only; no weights, no scoring.

Done when: new checks are tested the way the existing suite tests — inject the
liability by mutation, assert it is caught — and each has a measured base rate
on D2.

### Phase 2 — derive gates, then validate externally

1. Freeze the Level 1 candidate set from §5.1 plus the mechanism-justified
   Origin-1 Critical motifs.
2. Build D3 and D4. **Only now.**
3. Measure false-rejection rate of the frozen gate set on D3, with CIs.
4. Accept a gate into Level 1 only if the **upper** confidence bound of its
   false-rejection rate is below 5%.

Done when: there is a table of gates with FRR and CI on a population they were
not fitted to, and an explicit accept/reject decision per gate.

This is the step that converts the pipeline from a hypothesis into a rule. If
`TAP any RED` holds at 0% across ~1,048 unseen therapeutics, that is the
confidence the whole design rests on.

### Phase 3 — restructure fusion into three levels

1. Replace `screen_candidate`'s score-and-threshold with the §3 tree.
2. `CandidateResult` grows `triage_level` and `resource_count`; `verdict` is
   derived from `triage_level`, not from `score`.
3. Keep `score` as a diagnostic field. Route nothing on it.
4. `format_report` shows the level and the specific gate crossed, not a number.
5. Check conditional independence among Level 2 checks before combining any
   log-LRs; if correlated, report them separately rather than summing.

Done when: the existing test suite passes, and a candidate rejected at Level 1
reports which boundary it crossed.

### Phase 4 — non-human frameworks

Derive a separate gate set on D5. The germline reference is already the
patient's species (human by default), so murine input carries a large
systematic immunogenicity offset. Do not reuse human-framework gates unchanged.

### Phase 5 — de novo (blocked on data)

See §7.1. Cannot start without D6.

---

## 7. Limitations

### 7.1 De novo input — the important one

**Every number in this document was measured on natural or engineered human
antibodies. A sequence from RFDiffusion, an inverse-folding model, or a
protein language model is out of that distribution, and none of the measured
rates transfer.**

What that means concretely, if a user feeds in RFDiffusion-derived designs
tomorrow:

**The gates will still run, and they will still be meaningful — but for a
different reason than the numbers suggest.** The LR-based justification comes
from natural antibodies. The mechanism-based justification — "IMGT position 104
is the second cysteine of the intradomain disulfide; without it there is no
immunoglobulin fold" — does not depend on the calibration population at all.
Level 1 gates are mechanism-justified first and LR-justified second, which is
deliberate, and it is why they are the ones that survive the distribution
shift.

**The measured false-rejection rate does not transfer.** "0/150 known-good
antibodies rejected" tells you the gate does not reject *natural* antibodies.
It says nothing about how often it rejects a *good de novo design*. A de novo
binder might legitimately sit outside the TAP green region and still work. **We
have no evidence either way, and the report must not imply we do.**

**The rates will move in opposite directions, and this is the interesting
part.** For natural antibodies, `developability` has LR = 1.00 because every
natural antibody's hydrophobicity sits inside the same narrow band —
evolution put it there. De novo designs are under no such constraint. The same
check that carries zero information on naturals may carry real information on
de novo output. Symmetrically, the Origin-1 poly-residue motifs occur ~0% in
naturals but are a known generative-sampling artefact.

So the honest summary is:

> **De novo triage needs its own derivation study. The gates transfer as
> hypotheses; the numbers do not transfer at all. A check that is useless on
> natural antibodies may be the most useful one on de novo output, and vice
> versa.**

Until D6 exists, any de novo run must be labelled in its output as
**uncalibrated for this input distribution**. That label is not optional
politeness; without it the user will read a 0% false-rejection rate that was
never measured on their data.

There is one further wrinkle specific to co-folding pipelines. A design that
comes with a predicted co-structure has already been filtered by the generator
on interface quality (ipTM or similar). Conditioning on that changes the prior,
so the LRs — which are ratios conditional on population membership — are being
applied to a population that was selected by a correlated criterion. Phase 5
must account for this rather than assuming independence.

### 7.2 Limitations of the label itself

The D2/D3 label is provenance, not an assay. "Reached the clinic" is not
"manufacturable"; "was patented" is certainly not "bad" — most patented
antibodies are real programmes. The background population is therefore not
truly negative, which caps any separation measurable here and means LR⁺ values
are likely **underestimates** of true discrimination against a genuinely bad
population.

This cuts both ways and should not be used as an excuse: it makes a low LR
ambiguous (weak check, or weak label?), but it makes a *high* LR trustworthy —
discriminating against a population that is mostly fine is harder, not easier.

### 7.3 Statistical limitations

- Zero-cell LRs (0/150) are bounded only from one side. The Haldane–Anscombe
  correction gives a point estimate; the confidence interval is what matters
  and Phase 0 must report it.
- Log-LR addition assumes conditional independence, which our checks violate.
- 150 per group supports estimating a ~3% event rate poorly. D3 at ~1,048 is
  the fix for Level 1 gates specifically.

### 7.4 What this plan does not address

Expression titer, viscosity, self-association, thermal stability, polyreactivity,
and actual binding. Triage decides what not to test. It cannot tell you what
will work.

---

## 8. Where to start

Phase 0, §6. Two files, no behaviour change, and it makes every claim in §1
reproducible rather than something a previous session asserted.

The one thing to protect above all: **do not look at D3 until the Phase 2 gates
are frozen.** It is the only independent population we have, and it can only be
spent once.
