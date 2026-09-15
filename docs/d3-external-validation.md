# D3 external validation — pre-registration

**Written BEFORE the run. Nothing below may be edited after the result is
known; corrections go in an appended section.**

Gate set frozen at commit `a19808e`. The point of this document is that a
prediction recorded before a measurement is a test, and a prediction recorded
after one is a description.

## What is being validated

Only one number: **the false-rejection rate of the Level 5 gate set on
clinical-stage therapeutic antibodies that no part of this pipeline has seen.**

Every likelihood ratio in the repository so far was measured on D2 — the same
300 antibodies the checks were inspected and adjusted on. That makes the
pipeline a hypothesis. This run is what converts it into a rule, or refutes it.

## Frozen gate set

A candidate routes to Level 5 if any of these fires:

| Check | Tier | Rejecting condition | Justified by |
|---|---|---|---|
| `sequence_sanity` | 1 | hard gate: non-standard amino acid, or CDR-H3 length outside [4, 30] | measured on D2: 0/150, LR+ 11.00 |
| `v_domain_integrity` | 1 | any conserved IMGT anchor (23/41/104/118) absent or substituted | mechanism: no intradomain disulfide, no fold |
| `poly_residue_run` | 1 | AbSci Origin-1 Critical tier: 5xW, 4xG, 5xY, 5xS in one CDR | mechanism and precedent; LR unmeasurable on natural antibodies |
| `tap_psh` | 2 | RED | measured on D2: 0/150, LR+ 7.00 |
| `tap_pnc` | 2 | RED | measured on D2: 0/150, LR+ 5.00 |
| `tap_cdr_length` | 2 | RED | zero RED events on D2 in either group |
| `tap_ppc` | 2 | RED | zero RED events on D2 in either group |
| `tap_sfvcsp` | 2 | RED | zero RED events on D2 in either group |

## Population

**D3** = every `pairing == "TheraSAbDab"` entry in PLAbDab, **excluding the 150
used in D2**.

Sizing, corrected before the run: PLAbDab holds 1,198 TheraSAbDab rows, 1,069
after deduplicating on (heavy, light) and dropping fragments. Removing D2's 150
leaves **919**. The plan document said "~1,048", which was the raw count minus
150 without deduplication.

Both arms are scored, because the compute is the same either way:

- **primary**: the human-framework subset, which is what D2 was filtered to and
  therefore the only arm comparable to the derivation numbers
- **secondary**: the non-human-framework remainder, reported separately. The
  gates are supposed to work on these too, but nothing was derived on them.

Deliberately not included: a background arm. This run measures under-triage
only. Over-triage needs D4 and is a separate question.

## The standard being tested

American College of Surgeons Committee on Trauma: under-triage below 5%. Our
budget is the upper 95% confidence bound, not the point estimate — which is
the rule that demoted the CDR N-glycosylation gate on D2.

## Pre-registered predictions

Recorded before the run, with the D2 point estimate each is extrapolated from:

| Quantity | D2 observed | D3 prediction |
|---|---|---|
| Level 5 rate, human-framework arm | 2.7% (4/150) | **2–5%**, upper 95% bound below 5% |
| `tap_psh` RED | 0/150 | 0–1% |
| `tap_pnc` RED | 0/150 | 0–1% |
| `poly_residue_run` Critical | 2.0% (3/150) | 1–4% — the largest single contributor |
| `v_domain_integrity` | 0.7% (1/150) | 0–2% |
| `sequence_sanity` hard gate | 0/150 | 0–1% |

**Pass condition**: upper 95% bound of the overall Level 5 rate < 5%.

**What failure would mean, stated in advance so it cannot be rationalised
afterwards**: if the rate lands materially above 5%, the gate set does not
generalise from the 300 antibodies it was derived on, and the honest response
is to demote whichever gate carries the excess — not to widen the budget, and
not to re-sample D3.

## Deviation from a fully blind test, declared

`run_immunogenicity=False` for this run. The immunogenicity checks are
report-only (weight 0, and they route to no level), so excluding them cannot
change a single Level 5 assignment. They are skipped because 1,048 candidates
against a live IEDB endpoint is the exact setup that produced order-correlated
throttling during Step 7, and re-introducing that risk buys nothing here.

D2's numbers were produced with immunogenicity enabled. Since those checks
route nothing, the Level 5 comparison is unaffected; the Level 1-4
distribution is comparable but not identical, and is reported as secondary.

---

# RESULT — appended after the run

Nothing above this line was edited. Gate set as frozen at `a19808e`;
per-candidate data in `d3-external-validation-results.csv`.

## The pre-registered pass condition was not met

| Arm | n | Level 5 | 95% CI | Budget (upper bound < 5%) |
|---|---|---|---|---|
| **human framework (primary)** | 669 | **5.53%** (37) | 4.04 – 7.53% | **FAIL** |
| non-human framework (secondary) | 250 | 2.00% (5) | 0.86 – 4.60% | pass |

D2 predicted 2–5% on the primary arm. The observed 5.53% exceeds the budget on
the point estimate, not merely on the upper bound — the gate set rejects
roughly twice as many unseen clinical-stage therapeutics as the derivation
sample suggested it would.

## Which gate carried the excess

Per-gate rates on the 669 human-framework therapeutics:

| Gate | D3 rate | D2 rate |
|---|---|---|
| **poly_residue_run** | **19.88%** (any tier) | 2.0% |
| tap_psh RED | 0.90% | 0.0% |
| v_domain_integrity | 0.75% | 0.7% |
| tap_pnc RED | 0.60% | 0.0% |
| tap_cdr_length RED | 0.30% | 0.0% |
| sequence_sanity | 0.15% | 0.0% |
| tap_ppc RED | 0.15% | 0.0% |
| tap_sfvcsp RED | 0.00% | 0.0% |

## The real error: the budget was applied per gate, not to the union

Broken down per motif, **every AbSci Critical liability clears 5% on its own**:

| Motif | D3 rate | 95% upper |
|---|---|---|
| 5×W | 0.00% | 0.57% |
| 4×G | 0.45% | 1.31% |
| 5×Y | 1.35% | 2.54% |
| 5×S | 1.35% | 2.54% |

And yet the set failed. A candidate is rejected if **any** gate fires, so the
quantity the 5% standard bounds is the union — and eight gates at roughly 1%
each is 8%. Both the plan document and the pre-registration stated the rule
per check ("a check may only enter Level 5 if its false-rejection rate clears
5%"), which is not the constraint that governs what actually happens to a
candidate.

This is the most useful thing the run produced, and it would not have shown up
in any amount of review of the derivation numbers: on D2 the union happened to
stay under budget, so the distinction never bit.

## Action taken — the one recorded in advance

`poly_residue_run` demoted out of the rejecting set to review-only. It
contributed 26 of the 37 rejections. Recomputing the union with the remaining
gates:

| Arm | n | Level 5 | 95% CI | Budget |
|---|---|---|---|---|
| human framework | 669 | **2.54%** (17) | 1.59 – 4.03% | **PASS** |
| non-human framework | 250 | 1.60% (4) | 0.62 – 4.04% | pass |

No threshold was widened and D3 was not re-sampled, which the pre-registration
ruled out in advance.

**This 2.54% is not an independent validation.** It is the post-demotion rate
on the same data that prompted the demotion. The surviving gate set has been
validated to the extent that its failure mode was found and removed; confirming
the corrected set needs a population it has not touched.

## What the demotion does and does not mean

It is not a verdict on the motifs. It is a statement about which population
they describe. AbSci filters its Critical tier outright when selecting designs
out of a generative model, and that remains right: a run of five serines in
CDR-H3 is a sampling artefact a language model produces and evolution does
not. But **3.9% of antibodies that reached the clinic carry one**, so against
natural and humanised sequences the motifs are a description rather than a
defect.

For de novo input this check should be a gate again. Its false-rejection rate
on de novo designs is unmeasured, and re-enabling it without measuring that
would repeat the mistake this run caught.

## Secondary observations

The non-human arm passed at 2.00% and was never derived on, which makes it the
closest thing here to a clean test of the non-poly gates. Its level
distribution is very different — 70.0% at Level 3 (Immune risk) against 15.8%
for the human arm — which is the germline-identity check doing exactly what it
is for.

`tap_psh` RED fired on 0.90% of unseen therapeutics after 0/150 on D2. Still
comfortably inside budget, but a reminder that a zero on 150 does not mean a
zero.

## What is now spent, and what is left

D3 is spent. The derivation set D2 was already spent. **The corrected gate set
has no unused population left to validate it against.** D4 (background arm) can
measure over-triage but not under-triage, since it contains no known-good
antibodies.

Options, in order of preference: PLAbDab's `Xtal structure` pairing class as a
third known-good-ish population; a fresh PLAbDab release; or wet-lab outcomes.
