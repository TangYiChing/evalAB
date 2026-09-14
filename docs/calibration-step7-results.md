# Step 7 calibration: results

Run date: this commit. Raw per-candidate data: `calibration-step7-results.csv`.
Reproduce with:

```bash
python -m antibody_prescreen.calibration.run_calibration \
    --n-per-group 150 --workers 2 --seed 0 --out calibration_results.csv
python -m antibody_prescreen.calibration.fit_thresholds calibration_results.csv
```

## Headline

**The immunogenicity check carries no information about the label and has been
set to report-only (weight 0). The combined score separates the two
populations only weakly (AUC 0.587). Thresholds were re-anchored to a
clinical-grade population, which improves what the scale *means* without
making it a classifier.**

## Population

300 human-framework antibodies from PLAbDab, 150 per group:

| Group | `pairing` | Meaning |
|---|---|---|
| therapeutic | `TheraSAbDab` | reached clinical stage |
| background | `Patent text` | patented, no developability filter applied |

### Why human-framework only

Not because species confounds the populations — it barely does:

| | both-chain human framework |
|---|---|
| TheraSAbDab | 74.8% |
| Patent text | 70.4% |

But because the germline reference is the patient's species (human), so a
murine framework correctly reads as almost entirely non-self and carries a
large systematic offset. Leaving those in would let "is it murine" do part of
the discrimination — a real signal, but a different one from developability.
`--include-non-human` keeps them.

PLAbDab's own `organism` column is unusable here: it is `Unknown.` for all
92,222 rows across both populations, so species must be derived from ANARCI's
germline assignment.

## Results

### Per-component separation

AUC is the probability that a random background antibody scores worse than a
random therapeutic. 0.5 = the component says nothing about the label.

| Component | therapeutic median | background median | AUC |
|---|---|---|---|
| Tier 1 | 29.41 | 30.87 | 0.578 |
| TAP weight | 0.00 | 0.00 | 0.552 |
| **Immunogenicity weight** | 10.00 | 10.00 | **0.500** |
| Total (all three) | 41.95 | 43.23 | 0.552 |

Combinations:

| Combination | AUC |
|---|---|
| Tier 1 only | 0.578 |
| **Tier 1 + TAP** | **0.587** (permutation p = 0.005) |
| Tier 1 + TAP + immunogenicity | 0.557 |
| TAP only | 0.552 |
| Immunogenicity only | 0.500 |

Adding immunogenicity makes the combined score **worse**.

### Immunogenicity carries nothing

| | therapeutic | background | AUC |
|---|---|---|---|
| non-germline CDR binders (mean) | 1.93 | 1.94 | 0.505 |
| non-germline framework binders (mean) | 0.53 | 0.47 | 0.488 |
| observed T-cell assay matches (mean) | 0.03 | 0.05 | 0.510 |

The reading: after the germline correction, what remains is predicted MHC-II
binding in CDRs — and that is about equally common in antibodies that reached
the clinic and antibodies that did not. Which is not surprising. Clinical
antibodies are not selected against *predicted* epitopes, and real ADA rates
depend on dose, route, duration and patient HLA, none of which a sequence
scan sees.

The check was not wasted: the germline correction is what makes its output
meaningful at all (uncorrected it fires on every antibody), and a specific
CDR epitope with published positive T-cell assay records is worth a human's
attention. But it should not move an automated score.

### The one clean categorical signal

| | any TAP RED flag |
|---|---|
| therapeutic | **0 / 150** |
| background | **5 / 145** |

Fisher exact, two-sided: **p = 0.028**.

No clinical-stage antibody in this sample trips a single TAP red flag.
Sensitivity is low (it catches 3.4% of background), so it cannot carry a
screen by itself, but as a categorical "no approved-stage antibody looks like
this" marker it is the cleanest result in the run.

Per-metric rates:

| Metric | therapeutic red/amber | background red/amber |
|---|---|---|
| Total CDR length | 0.0% / 12.0% | 0.0% / 18.6% |
| PSH (hydrophobic) | 0.0% / 10.0% | 2.1% / 11.7% |
| PPC (positive) | 0.0% / 6.7% | 0.0% / 6.9% |
| PNC (negative) | 0.0% / 2.0% | 1.4% / 3.4% |
| SFvCSP | 0.0% / 2.0% | 0.0% / 5.5% |

Raw TAP values barely separate at all (AUC of absolute deviation from the
therapeutic median: 0.506–0.540). The discrete flags carry what little signal
there is; the continuous values do not.

## Changes adopted

1. **Immunogenicity weights -> 0.** `CDR_BINDER_WEIGHT`,
   `FRAMEWORK_BINDER_WEIGHT`, `OBSERVED_EPITOPE_WEIGHT` are all 0.0. Flags are
   still produced with full detail and appear in `all_flags`; they contribute
   nothing to the score. Raising them again requires a population that shows
   them predicting something.
2. **`T_LOW` 31.0 -> 29.41, `T_HIGH` 40.0 -> 38.31.** Percentile placements on
   the therapeutic population (median and p90). The numbers barely moved; what
   changed is what they are anchored to — "antibodies that reached the clinic"
   rather than "antibodies someone crystallised".
3. **TAP weights unchanged.** They contribute a small but statistically real
   improvement (0.578 -> 0.587, p = 0.005) and are left alone.

Routing under the adopted thresholds:

| | GO | CONDITIONAL | NO-GO |
|---|---|---|---|
| therapeutic | 50.0% | 40.0% | 10.0% |
| background | 40.7% | 40.0% | 19.3% |

### One behaviour change worth knowing

The repository's own reference test pair (`VH_REF`/`VL_REF`) scores 39.92 and
is now **NO-GO**, where it was CONDITIONAL under `T_HIGH = 40.0`. It carries
several CDR deamidation motifs and sits in the worst decile of the
therapeutic population. That is a legitimate result, not a threshold bug —
and the README already records that this pair is not a verified match to any
approved drug.

`test_clean_pair_is_not_no_go_on_hard_gates` asserted a verdict, which coupled
it to threshold constants; it now asserts the hard-gate property its name
describes, and is renamed `test_clean_pair_is_not_hard_gated`.

## What this does not establish

**The label is provenance, not an assay.** "Reached the clinic" is not
"manufacturable", and "was patented" is certainly not "bad" — most patented
antibodies are real programs from real companies, not developability
failures. A background population that is not actually negative puts a
ceiling on any separation measurable here, and AUC 0.587 is consistent with
that ceiling rather than with the checks being uninformative.

So: the score is a **triage ordering**, not a classifier. Use it to rank a
batch and to decide what to read first. Do not use a threshold crossing as a
decision by itself, and read the flags — especially any TAP red — rather than
the number.

A stronger calibration needs a population with measured outcomes:
developability assays, manufacturing results, or clinical ADA rates. That
data is not in PLAbDab.

## Data-quality note: a bias that nearly invalidated this run

The first completed run had immunogenicity for only 205/300 candidates, split
**100% (150/150) for therapeutics and 36.7% (55/150) for background**.

That is not a property of the data. IEDB throttles under concurrent load, and
the two populations were scored in order, so the group that ran second
absorbed all the failures. Retrying those same sequences afterwards succeeded
immediately (515, 530, 520 rows in 3.7s, 2.8s, 3.0s).

Had it been analysed as-is, "background shows less immunogenicity" would have
been an artefact of not measuring it. A transient failure correlated with
scoring order is far more dangerous than a random one.

Fixed by adding retry with jittered exponential backoff, and by validating
responses **before** caching them — the original code cached the response text
first, which would have frozen a transient error page into the cache
permanently. IEDB's plain-text rejection of a sequence containing an
ambiguous residue (`X`, `B`, `Z`) *is* cached, since that is a reproducible
answer about that sequence rather than a transient failure.

Final coverage: 150/150 therapeutic, 149/150 background (the one miss is a
genuinely rejected ambiguous sequence).
