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
