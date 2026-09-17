# Range-based triage: design, calibration, and the first design batch

This is the design record for evalAB's range-and-repair triage: why the
system does not score candidates, what was measured to reach that conclusion,
and how the bands and the level matrix were derived.

It supersedes an earlier weighted-sum implementation, which is archived under
`tbd/legacy_score_pipeline/` for provenance. That code is outside the import
path, the test suite and the public workflow, and nothing in this document
imports it. Where the text below names `triage.py` or `fusion.py`, it is
referring to those archived files as the historical measurement context — not
to anything evalAB runs today.

## 1. What changed, in one sentence

The system stopped trying to say which antibody is bad, and started saying how
far outside the common range a candidate sits and what it would cost to bring
it back.

## 2. Why the old question was unanswerable

Nothing in this repository has an outcome label. Every attempt to separate
"good" from "bad" therefore had to borrow one, and the only labels available
were provenance — reached the clinic versus appeared in a patent — which
produced an AUC of 0.587. Statistically real, practically nothing.

The TAP episode is the clearest case, and it is what forced the redesign. TAP
reports GREEN/AMBER/RED. Those bands were drawn around the distribution of
known therapeutic antibodies. Measuring their likelihood ratio on a population
of known therapeutic antibodies therefore measured how much that distribution
resembles itself, and returned 1.00-2.43. That is not a weak signal; it is the
same ruler applied twice. No better statistic fixes it, because the problem is
not the statistic.

## 3. The question that replaced it

> Where does this candidate sit relative to antibodies we have seen, and what
> would it cost to move it back inside?

Descriptive, and answerable with a reference population and nothing else. Two
axes, and one table.

**Axis 1 — deviation (`bands.py`).** For each measured quantity, percentile
cut points fitted on a stated reference population:

| tier | where | meaning |
|---|---|---|
| 0 | inside p5-p95 | typical |
| 1 | p95-p99 | unusual |
| 2 | p99-p99.9 | strongly unusual |
| 3 | beyond p99.9 | outside the range the reference spans |

Two-sided by default. A CDR-H3 of 4 residues is as unusual as one of 28.

**Axis 2 — repair cost (`profile.py`).** What undoing the finding costs, and
what it costs you:

| cost | what a person does |
|---|---|
| 0 | nothing |
| 1 | one framework substitution; binding untouched |
| 2 | one CDR substitution; affinity must be re-measured |
| 3 | redesign a loop or a surface — not a substitution |
| 4 | no repair path by substitution |

**The table (`levels.py`).** `level = max(MATRIX[deviation][cost])` over every
finding. Still a max, still no weights.

```
                         cost ->
                   0     1       2        3         4
    deviation 0    1      1      2        3         4
              1    1      2      2        3         4
              2    2      2      3        4         5
              3    3      3      4        5         5
```

No check name appears anywhere in the routing logic. Adding a check means
answering the two questions every other check answered, not arguing for a new
routing rule.

## 4. Two bugs the first runs found, and what they taught

**(a) A finding that is typical AND expensive is not a finding.** The first
calibration put 98% of clinical-stage therapeutics in Level 3. The cause was
an oxidation-prone tryptophan in FR3 being emitted as deviation 0, cost
REDESIGN — on every antibody ever made. Fixed in `profile.py` by not emitting
it (tryptophans still count toward the banded `cdr_oxidation_sites` metric),
not in the table, so the table stays monotone in both arguments.

**(b) A tier cannot be finer than the reference population's resolution.** The
first design run put 23 of 23 candidates at Level 5, all on `cdrh1_len=7`. The
TheraSAbDab reference has p0.1 = p1 = p5 = 8 for that metric: the entire lower
tail is one number. The p99.9 cut was the minimum of the sample wearing a
percentile's name, because 534 x 0.001 = 0.53 antibodies.

`bands.MIN_TAIL_OBS` now requires a tail to hold at least 5 reference
antibodies before a cut placed in it means anything, and each band records the
finest tier its `n` can assert. **Consequence, stated plainly: with a reference
of ~1,000 antibodies, no range metric can reach deviation tier 3.** Tier 3
becomes reachable only from mechanism — a non-standard amino acid, a broken
fold anchor. Wanting tier 3 from a distribution means wanting a reference of
~5,000 or more; it is not a reason to pretend 534 is enough.

## 5. Calibration

Reference: `PLAbDab/TheraSAbDab`, 1,069 deduplicated pairs, split
deterministically into 534 derivation / 535 validation. **Cut points are fitted
on the derivation half only. Every rate below is measured on the validation
half**, which the cut points have never seen — the fix for critique (b) in
`triage_implementation_plan.md`, where the old `T_HIGH` rejection rate was a
restatement of where the line was drawn.

Reproduce:

```
python -m antibody_prescreen.calibration.fit_bands --reference therasabdab
```

16 bands fitted. Level distribution on the validation half:

| level | n | % |
|---|---|---|
| L1 Ready | 37 | 6.9% |
| L2 Point repair | 207 | 38.7% |
| L3 Redesign | 223 | 41.7% |
| L4 Human decision | 67 | 12.5% |
| L5 Out of range | 1 | 0.2% |

**Two budgets, because there are two different mistakes.** The 5% figure is
borrowed from trauma under-triage, where the bounded mistake is sending a
severely injured patient home. Its analogue here is L5 — "do not pursue" — the
only level that discards a candidate. L4 hands the candidate to a person, which
is a cost, not a loss.

| bound | rate | 95% CI | upper bound vs 5% |
|---|---|---|---|
| L5 (discard) | 0.19% | 0.03-1.05% | **PASS** |
| L4+L5 (discard or escalate) | 12.71% | 10.15-15.80% | **FAIL** |

Which of these is the pass condition is a decision that has to be made before a
run, not after. Recorded here unresolved rather than settled silently.

The single L5 was a substituted conserved IMGT anchor — mechanism, not
distribution, exactly as §4(b) predicts.

## 6. The design batch

23 designed VH variants sharing one light chain and one framework, varying in
CDR-H3.

**These sequences are internal and are not in this repository.** The numbers
below are reported as a record of what the design produced on a real batch;
they are not reproducible from a clone, and nothing in the shipped code or data
depends on them. The reproducible worked example is `data/example/`, which uses
public data only.

The run was the equivalent of:

```
python examples/screen_designs.py <internal_batch>.csv \
    --shared-scaffold --unseal Status
```

**Result: 23 of 23 at Level 4.** All on findings inherited from the shared
scaffold — `cdrh1_len=7` (never seen in 534 clinical antibodies) and
`cdrh3_len=25` (above p99). The parent molecule is a Level 4 antibody, so every
child is.

This is a real property of a max-based triage on a congeneric series, not a
defect: the level is a true statement about each candidate and a useless one
for choosing between them. `screen_designs.py` therefore reports the batch
twice — once as-is, and once with the 11 findings common to all 23 set aside.
Nothing is re-scored; the same findings are partitioned by whether the design
introduced them.

With the scaffold set aside: 1 candidate at L1, 19 at L3, 3 at L4 — driven by
Fv net charge, Fv cysteine count, isoelectric point, and odd cysteine parity.

## 7. Held-out outcome labels

The `Status` column (7 "Failed purity QC", 16 "Ready") was not read by any code
path between profiling and level assignment. Unsealed after the run.

**The level does not separate them.**

| | n | full levels | introduced-only |
|---|---|---|---|
| Failed purity QC | 7 | all L4, mean 4.00 | all L3, mean 3.00 |
| Ready | 16 | all L4, mean 4.00 | mean 3.06 |

**One individual finding does.** Post-hoc, exploratory, on 23 congeneric
sequences from one parent — this is an observation, not a validated check, and
no threshold in this repository was moved because of it:

```
VH cysteine parity      failed   ready
  odd                        4       0
  even                       3      16
  Fisher one-sided p = 0.0040
```

Every other finding source was null or inverted (`fv_pi` p=0.98,
`fv_net_charge` p=0.99, `fv_cys_count` p=1.00).

The relevance is that `cysteine_pairing` was demoted to report-only in the
archived score pipeline (`tbd/legacy_score_pipeline/.../triage.py`) because it
fires on 5.3% of clinical-stage therapeutics — over the rejection budget. That demotion was correct for the population it was measured
on and appears to be wrong for this one. Natural selection does not produce
unpaired cysteines, so the check has no variance in a natural-antibody
reference and its measured likelihood ratio there says nothing about what it
does on generative output.

**This is the same error as the TAP circularity, one level up: not the wrong
statistic, the wrong population.** A check aimed at de novo designs cannot be
calibrated on antibodies that were never designed.

## 8. What this run does not establish

- 23 sequences from one parent is not a population. The cysteine result needs
  an independent batch.
- "Failed purity QC" is one assay at one step. It is not developability.
- Tier 2 (TAP, structure) was not run. Raybould's published AMBER/RED bands
  remain inherited, with the circularity of §2 unresolved for those four
  metrics — by decision, and it is recorded as a decision.
- No range metric can currently reach tier 3. Whether that is acceptable
  depends on a reference-population choice that has not been remade.

---

# 9. Appended: widening the ranges, and what it cost

Raised after the run above: CDR loop length is something antibody engineering
changes routinely, and framing a length change as an excursion against 534
clinical antibodies is the wrong frame.

That is correct, and the fix is not a wider multiplier. It is that **CDR length
is a property of the repertoire, not of what survived a development pipeline.**
Measuring it against TheraSAbDab asks "what do clinical-stage drugs look like"
and then reads the answer as "what do antibodies look like".

## 9.1 The two populations, side by side

Refitted on 40,000 deduplicated PLAbDab pairs (20,000 derivation / 20,000
validation), against the same 534/535 TheraSAbDab split:

| metric | TheraSAbDab p5-p95 | PLAbDab p5-p95 | TheraSAbDab p99.9 | PLAbDab p99.9 |
|---|---|---|---|---|
| cdrh1_len | 8 - 9 | 8 - 10 | 10 | 12 |
| cdrh2_len | 7 - 10 | 7 - 10 | 12.3 | 12 |
| **cdrh3_len** | 8 - 18.4 | 8 - **21** | 24.4 | **63** |
| cdrl1_len | 6 - 12 | 6 - 12 | 12 | 12 |
| cdrl2_len | 3 - 3 | 3 - 3 | 3.47 | 7 |
| cdrl3_len | 8 - 11 | 8 - 11 | 13.5 | 14 |
| fv_cys_count | 4 - 5 | 4 - **6** | 6 | 10 |
| fv_pi | 4.95 - 9.26 | 4.83 - 9.23 | 9.64 | 9.81 |
| fv_net_charge | -3.27 - 5.89 | -4.26 - 5.84 | 8.78 | 9.84 |
| vh_aggregation | 0.154 - 0.238 | 0.154 - 0.243 | 0.288 | 0.285 |
| vl_aggregation | 0.140 - 0.211 | 0.139 - 0.216 | 0.252 | 0.255 |
| vh_germline_identity | 0.764 - 0.99 | **0.724** - 0.99 | — | — |
| vl_germline_identity | 0.804 - 1.0 | **0.739** - 1.0 | — | — |

Three readings, and the third is the one worth keeping.

**Loop lengths and cysteine count move a lot.** CDR-H3's p99.9 goes from 24.4
to 63. Six cysteines in an Fv is above p95 clinically and inside it in the
repertoire.

**Surface properties barely move at all.** `vh_aggregation` has a p5 of 0.154
in both populations, to three significant figures. pI, net charge and
aggregation propensity are nearly identical in 534 clinical antibodies and
20,000 mostly-patent ones. That is a finding in its own right: on these
sequence-level measures, the development pipeline does not measurably narrow
the distribution relative to the patent corpus. It also means the population
choice is nearly free for these metrics — the same lines, with enough n to
support tier 3.

**Germline identity moves, and that is exactly why it must not be widened.**
PLAbDab's patent corpus contains murine and chimeric sequences. Calibrating a
humanness band on a population that includes them widens the band precisely far
enough to stop catching them. It stays on the clinical population.

## 9.2 What §9.1 led to — and what it did NOT settle

> **Superseded by §10.2.** The composition described here was an intermediate
> arrangement: TheraSAbDab as the base, with loop lengths and cysteine count
> taken from PLAbDab. The shipped configuration is the *inverse* — PLAbDab
> (human-only, cleaned) as the base, with only germline identity taken from
> the clinical population. Read §10.2 for what evalAB actually runs. This
> section is kept because the measurement below is what motivated the change.

The reading that survived from §9.1 is the third one: loop lengths and cysteine
count move a lot between the two populations, surface properties barely move at
all, and germline identity must not be widened. Once surface properties were
shown to be nearly population-independent, there was no longer a reason to pay
TheraSAbDab's resolution cost for them — which is what turned the composition
around. Note also that §9.1's table was fitted on 40,000 pairs **without** the
`--human-only --clean` filters; the frozen bands in §10.3 were not.

Validation half, intermediate composition vs TheraSAbDab-only:

| | TheraSAbDab only | intermediate composition |
|---|---|---|
| L5 (discard) | 0.19% (PASS) | 0.19% (PASS) |
| L4+L5 | 12.71% | **10.28%** |
| `cdrh3_len` as an L4/L5 driver | 2.62% | **0.56%** |

## 10.1 Declared scope

**Human antibodies, supplied as a paired VH/VL** — from an IgG, an scFv, or a
paired-chain design. Nothing else. Non-human germlines are explicitly out of
scope rather than silently mis-measured: a candidate whose closest germline is
not human raises `out_of_scope_species` and routes to Level 4, with a message
saying every band was fitted on human antibodies and does not describe it.

That is a deliberate choice over letting `germline_identity` handle it. The
germline identity band is fitted on 393 antibodies and therefore tops out at
deviation tier 1, so a murine framework would come back "unusual" — which is
not what is true about it. What is true is that nobody calibrated anything for
it, and Level 4 says that.

## 10.2 The reference populations

Both drawn from PLAbDab's `paired_sequences.csv.gz` (78,267 pairs after
deduplicating on (heavy, light) and dropping fragments).

**Base — cleaned human repertoire.** 58,177 profiled, split 29,147 derivation /
29,030 validation. Two filters, both `--human-only --clean`:

| filter | dropped (derivation) | why |
|---|---|---|
| non-human germline call | 8,928 | declared scope |
| malformed | 1,058 | `non_standard_aa`, `anchor_missing`, `anchor_substituted` — truncated or broken database rows, not antibodies |

`odd_cysteine` is **deliberately not** a disqualifier. An odd cysteine count is
unusual biology, not a malformed record; it is the one observation here with a
measured association to a wet-lab outcome (§7); and excluding it would bias
`fv_cys_count`, which is a fitted band. Removing it would quietly assume the
conclusion it is evidence for.

**Override — human clinical.** TheraSAbDab, same filters, 393 derivation / 420
validation. Used for `bands.CLINICAL_METRICS` only: `vh_germline_identity` and
`vl_germline_identity`. A humanness band fitted on a corpus containing murine
sequences widens exactly far enough to stop catching them, and "human" here is
ANARCI's germline call rather than a provenance record, so a humanised antibody
with murine back-mutations passes the filter. Anchoring the band to clinical
antibodies is what makes it encode "human enough to dose a person".

Everything else comes from the base, on the evidence in §9.1: loop lengths and
cysteine count are repertoire properties and move a lot between the two
populations, while surface properties barely move at all (`vh_aggregation` p5 =
0.154 in both, to three significant figures) — so the larger population wins on
resolution for free. 29,147 antibodies support deviation tier 3; 393 support
only tier 1.

## 10.3 The frozen bands

```
metric                 source        n     maxtier   p5      p95     p99.9
cdrh1_len              repertoire    29147    3       8      10      13
cdrh2_len              repertoire    29147    3       7      10      13
cdrh3_len              repertoire    29147    3       9      22      35.7
cdrl1_len              repertoire    29147    3       6      11      12
cdrl2_len              repertoire    29147    3       3       3       7
cdrl3_len              repertoire    29147    3       8      12      15
fv_cys_count           repertoire    29147    3       4       6       8
fv_pi                  repertoire    29147    3       4.87    9.25    9.82
fv_net_charge          repertoire    29147    3      -4.05    5.85    9.85
cdr_net_charge         repertoire    29147    3      -5.07    2.16    5.93
vh_aggregation         repertoire    29147    3       0.156   0.244   0.286
vl_aggregation         repertoire    29147    3       0.140   0.213   0.255
cdr_ptm_motifs         repertoire    29147    3       0       4       6
cdr_oxidation_sites    repertoire    29147    3       0       4       6
vh_germline_identity   clinical        393    1       0.768   0.99    —
vl_germline_identity   clinical        393    1       0.815   1.0     —
```

Files: `antibody_prescreen/data/bands_human_repertoire.json` and
`bands_clinical.json`. The uncleaned §9 intermediates are in `data/superseded/`
and must not be used.

## 10.4 Reproduce

```bash
python -m antibody_prescreen.calibration.fit_bands --reference therasabdab \
    --human-only --clean --out antibody_prescreen/data/bands_clinical.json

python -m antibody_prescreen.calibration.fit_bands --reference plabdab \
    --human-only --clean \
    --out antibody_prescreen/data/bands_human_repertoire.json \
    --clinical-bands antibody_prescreen/data/bands_clinical.json

python examples/screen_designs.py <batch.csv> \
    --bands antibody_prescreen/data/bands_human_repertoire.json \
    --clinical-bands antibody_prescreen/data/bands_clinical.json
```

## 10.5 What the frozen system does, on 29,030 unseen human antibodies

| level | n | % |
|---|---|---|
| L1 Ready | 2,112 | 7.3% |
| L2 Point repair | 10,301 | 35.5% |
| L3 Redesign | 13,086 | 45.1% |
| L4 Human decision | 3,116 | 10.7% |
| L5 Out of range | 415 | 1.4% |

| bound | rate | 95% CI | vs 5% upper-bound rule |
|---|---|---|---|
| **L5 (discard)** | **1.43%** | 1.30-1.57% | **PASS** |
| L4+L5 (discard or escalate) | 12.16% | 11.79-12.54% | FAIL |

**The pass condition is L5, and it passes.** The 5% figure comes from trauma
under-triage, where the bounded mistake is sending a severely injured patient
home. Its analogue is the only level that discards a candidate. L4 hands the
candidate to a person, which is a cost, not a loss — bounding it at 5% would be
borrowing a number from a situation it does not describe. Recorded here as the
standing decision rather than re-argued per run.

After cleaning, **no mechanism observation appears among the L4/L5 drivers at
all** — the reference no longer contains broken records, so every excursion is
now distributional. Top drivers: `vh_aggregation` 2.22%, `cdr_net_charge`
1.84%, `fv_net_charge` 1.84%, `vl_aggregation` 1.84%, `fv_pi` 1.71%.

## 10.6 The design batch, final

23 designed VH variants, one shared light chain and framework, varying in
CDR-H3. All call human germline (IGHV1-69*01 / IGLV1-44*01), so all are in
scope.

**23 of 23 at Level 4**, on scaffold-inherited findings — `cdrh1_len=7` (below
p1 of 29,147 human antibodies, above p0.1: rare, not unprecedented) and
`cdrh3_len=25` (above p95). The parent is a Level 4 antibody, so every child
is. With the 11 shared findings set aside:

| introduced-findings level | n |
|---|---|
| L1 nothing beyond the scaffold | 10 |
| L2 one CDR point repair | 2 |
| L3 redesign | 11 |
| L4 | 0 |

## 10.7 Held-out outcome: the level still does not predict purity QC

`Status` was never read by any code path between profiling and level
assignment. Unsealed after the final run.

| | n | full | introduced-only |
|---|---|---|---|
| Failed purity QC | 7 | all L4 | mean **2.14** |
| Ready | 16 | all L4 | mean **2.00** |

Failures sit marginally *above* successes now, having sat marginally below
before the widening. On 7 versus 16 both directions are noise. The honest
statement is unchanged across all three configurations tried: **the level does
not track purity QC on this batch.**

The one thing that did, across every configuration, is the cysteine-parity
observation from §7 — 4 of 4 odd-cysteine candidates failed, 0 of 16 Ready
candidates had one, Fisher one-sided p = 0.0040. No band touches it, because
parity is a mechanistic observation, so no reference-population change moved it.

This file was never used to tune anything. No threshold, cut point, cost
assignment or filter in this repository was chosen by looking at it.

## 10.8 What remains true and unfinished

- **The cysteine result needs an independent batch.** 23 congeneric sequences
  from one parent, examined post-hoc, is an observation. It is the cheapest
  high-value experiment available: one batch from a different parent, scoring
  nothing but VH cysteine parity.
- **Tier 2 is not run.** TAP's four structure metrics still use Raybould's
  published AMBER/RED bands, so the circularity of §2 is unresolved for them —
  by decision, recorded as a decision.
- **"Failed purity QC" is one assay at one step.** It is not developability.
- **Germline identity cannot exceed deviation tier 1**, because 393 antibodies
  cannot assert a 1-in-100 boundary. Out-of-scope species is handled
  mechanistically (§10.1); a *within-scope* humanness excursion tops out at
  Level 3. Fixing it means a larger human clinical reference, which does not
  currently exist in PLAbDab.
