# What evalAB measures

Every finding evalAB reports comes from this list. For each one: what it is,
what tool produces it, why it has anything to do with an antibody, and the
highest level it can drive on its own.

Nothing here is a prediction. These are descriptions of a sequence, located
against a stated reference population. Whether any of them predicts a wet-lab
outcome in your programme is a question your bench answers, not this table —
see [calibration.md](calibration.md).

---

## How to read this document

Findings come in two families, and the difference is not cosmetic.

**Measured metrics** (`RangeMetric`) are numbers with a distribution. Their
deviation tier is read off the reference population: p5–p95 is typical, and
the tiers above it are p95–p99, p99–p99.9, and beyond p99.9. Nobody's opinion
enters; the cut point is a percentile.

**Mechanistic observations** (`BinaryObservation`) are yes/no facts with no
distribution worth fitting. A conserved fold anchor is present or it is not.
Fitting a percentile to that would be theatre — the answer does not depend on
how often it happens, it depends on whether the immunoglobulin fold survives
it. Their deviation is asserted from mechanism and stated in the table.

Both carry a **repair cost**, which is a property of the finding and is fixed
at the moment of measurement:

| cost | meaning |
|---|---|
| 0 | nothing to undo |
| 1 | one substitution in framework; binding untouched |
| 2 | one substitution in a CDR; affinity must be re-measured |
| 3 | a loop, a charge distribution, a patch — not one residue |
| 4 | no substitution fixes it; the molecule is not this molecule |

Level is then `MATRIX[deviation][cost]`, and the candidate's level is the
**maximum** over its findings. See [range-triage-design.md](range-triage-design.md) §5.

### Why some metrics cannot reach Level 5: `max_tier`

A percentile is a claim about a tail, and a tail needs antibodies in it. evalAB
requires at least **5 real observations** beyond a cut point before it will use
that cut point (`bands.MIN_TAIL_OBS`):

| tier | cut point | tail fraction | antibodies needed |
|---|---|---|---|
| 1 | p95 | 5% | 100 |
| 2 | p99 | 1% | 500 |
| 3 | p99.9 | 0.1% | 5,000 |

`Band.max_tier` records the finest tier that band's population can support, and
`Band.tier()` caps every reading at it. So with n=393 clinical antibodies,
p99.9 is 0.4 antibodies out — a number interpolated from nothing. The band is
capped at tier 1 and the p99.9 column is decoration.

This is why the repertoire population (n=29,147) is used wherever the two
populations agree: the larger n buys tier-3 resolution for free. And it is why
germline identity, which must stay anchored to the clinical population, tops
out at Level 3 no matter how far from germline a candidate sits.

---

## Measured metrics

All fitted on PLAbDab, human-only, cleaned (n=29,147) except where noted.
"Two-sided" means unusually low is flagged as well as unusually high.

### Loop lengths

| metric | tail | cost | max level |
|---|---|---|---|
| `cdrh1_len` `cdrh2_len` `cdrh3_len` `cdrl1_len` `cdrl2_len` `cdrl3_len` | two-sided | 3 | **L5** |

**What:** number of residues in each CDR, by IMGT region boundaries.
**Tool:** ANARCI, IMGT numbering scheme (`numbering.py`).
**Biology:** the CDRs are the six loops that form the antigen-binding site.
CDR-H3 in particular is the product of V(D)J recombination and is the most
variable — it usually contributes the most binding surface. Loop length is a
structural parameter: a loop far outside the natural range may not pack against
the framework the way the fold expects, and long CDR-H3s are associated with
conformational flexibility.
**Why two-sided:** a 6-residue CDR-H3 is as unusual as a 28-residue one, and
both are worth a person seeing.
**Why cost 3:** you cannot change a loop's length with a point substitution.
Changing it is a redesign, and the binding surface changes with it.

### Germline identity

| metric | tail | cost | max level |
|---|---|---|---|
| `vh_germline_identity` `vl_germline_identity` | two-sided | 3 | **L3** (capped) |

**Fitted on TheraSAbDab clinical, n=393** — the only metrics that are. See
[range-triage-design.md](range-triage-design.md) §9.1 for the measurement, and
§10.2 for why the `--human-only` filter does not make this redundant.

**What:** fraction of residues identical to the closest human germline V gene.
**Tool:** ANARCI with `assign_germline=True`, which aligns against the IMGT
germline V gene database and returns the nearest gene plus its identity.
**Biology:** a B cell starts from a germline V gene written in your DNA — a
sequence your immune system has been tolerised to since before birth — and then
accumulates mutations through somatic hypermutation to raise affinity. Identity
measures how far the molecule has travelled from that starting point. The
further it is, the more likely a patient's immune system treats it as foreign,
which is the anti-drug antibody (ADA) risk. This is a sequence-level proxy for
immunogenicity, not a prediction of it.
**Why two-sided:** 100% identity across both chains is also unusual, and it
usually means the CDRs were never engineered.
**Important limitation:** "human" here is ANARCI's nearest-germline call, not a
provenance record. A humanised antibody carrying murine back-mutations passes
this and every human-only filter in the pipeline.

### Charge and isoelectric point

| metric | tail | cost | max level |
|---|---|---|---|
| `fv_pi` `fv_net_charge` `cdr_net_charge` | two-sided | 3 | **L5** |

**What:** `fv_pi` is the pH at which the Fv carries zero net charge.
`fv_net_charge` is net charge at pH 7.4 over the whole Fv; `cdr_net_charge` is
the same over CDR residues only.
**Tool:** `developability.py`, Henderson-Hasselbalch over the ionisable side
chains (D, E, C, Y negative; K, R, H positive) plus both termini; pI by
bisection on that same model. pI and net charge deliberately share one pKa
table, so the two can never contradict each other merely because they were
computed from different constants.
**Biology:** at pH equal to pI the molecule has no net charge and therefore no
electrostatic repulsion from its neighbours, which is where aggregation and
precipitation are most likely — so a formulation buffer must sit away from pI.
pI also determines ion-exchange purification behaviour, and a very high pI
(strongly positive at physiological pH) is associated with non-specific binding
to negatively charged cell surfaces and faster clearance. `cdr_net_charge` is
the sequence-level stand-in for TAP's charged-patch metrics: the patches
themselves need a structure, but the charge that forms them is already visible
in the loops.
**Caveat:** this is a sequence estimate. Real pKa values shift with burial,
three-dimensional context and neighbouring charges, so a measured cIEF value
will differ. Use the band for relative position, not as a replacement for the
assay.

### Aggregation propensity

| metric | tail | cost | max level |
|---|---|---|---|
| `vh_aggregation` `vl_aggregation` | two-sided | 3 | **L5** |

**What:** AGGRESCAN Na4vSS — the whole-chain aggregate of windowed
aggregation-propensity values above the hot-spot threshold.
**Tool:** `developability.py`, a port of the pure-Python AGGRESCAN
implementation (Conchillo-Sole et al. 2007 a3v values, 5-residue window).
**Biology:** stretches of hydrophobic, aggregation-prone residues on the
surface can drive self-association, which shows up as low SEC monomer content,
viscosity at high concentration, and shelf-life problems.
**Why the aggregate and not hot-spot counts:** nearly every real protein has a
handful of individual AGGRESCAN hot spots. That is normal, so counting them
flags everything. The whole-sequence magnitude is the comparison that carries
information.

### Counts of repairable liabilities

| metric | tail | cost | max level |
|---|---|---|---|
| `cdr_ptm_motifs` `cdr_oxidation_sites` | **one-sided (high)** | 2 | **L4** |

**What:** how many post-translational-modification motifs, and how many
oxidation-prone residues, fall inside CDRs.
**Tool:** `profile.py` motif scan over IMGT-assigned regions.
**Why one-sided:** nobody needs telling that zero deamidation motifs is
unusually few. Counts are flagged only from above.
**Note:** each individual motif is *also* reported as its own observation with
its own location — the count exists so that a candidate carrying many of them
is distinguishable from one carrying a single motif.

### Cysteine count

| metric | tail | cost | max level |
|---|---|---|---|
| `fv_cys_count` | two-sided | 3 | **L5** |

**What:** total cysteines in the Fv.
**Biology:** the immunoglobulin fold has one canonical intra-domain disulfide
per domain (IMGT 23–104), so four cysteines across a paired Fv is the norm.
Extra cysteines mean free thiols or non-canonical disulfides, which drive
covalent aggregation and heterogeneity during manufacture.
**Deliberately kept in the reference population:** `odd_cysteine` is *not* a
reference-disqualifying filter (see `fit_bands.REFERENCE_DISQUALIFYING`), both
because it is unusual biology rather than a malformed record and because
excluding it would bias this very band.

---

## Mechanistic observations

Deviation is asserted from mechanism, not fitted. Level is fixed.

| observation | deviation | cost | level | what it means |
|---|---|---|---|---|
| `non_standard_aa` | 3 | 4 | **L5** | a character that is not one of the 20 amino acids — nothing downstream can interpret it |
| `anchor_missing` / `anchor_substituted` at IMGT 23, 41, 104 | 3 | 4 | **L5** | a conserved fold anchor is gone: the two cysteines forming the intra-domain disulfide, or the buried tryptophan that packs the core. These hold the immunoglobulin fold together |
| `anchor_missing` / `anchor_substituted` at IMGT 118 | 2 | 3 | **L4** | the J-region anchor. Usually a truncated paste rather than a broken molecule, so it is scored one tier lower |
| `out_of_scope_species` | 2 | 3 | **L4** | closest germline is not human. Every band was fitted on human antibodies, so the candidate is not "unusual" — it is unmeasured. A person decides |
| `odd_cysteine` | 2 | 2 | **L3** | an odd cysteine count guarantees at least one thiol that cannot pair canonically. Sequence alone cannot tell a defect from a designed non-canonical disulfide, so it is not a rejection |
| `poly_residue_run` | 2 | 2 | **L3** | ≥5 W, ≥4 G, ≥5 Y or ≥5 S consecutively in a CDR (AbSci Origin-1 Supplementary Table 7, Critical tier). Occurs zero times in PLAbDab's therapeutic arm, which is why it has no distribution to fit — it is a generative-sampling signature |
| `ptm_motif` | 0 | 1 FR / 2 CDR | **L1 / L2** | deamidation (NG, NS, NT, NH) or isomerization (DG, DS). Chemical degradation during storage; in a CDR it can cost binding |
| `n_glycosylation` | 0 | 1 FR / 2 CDR | **L1 / L2** | N-X-S/T sequon, X≠P. An unintended glycan adds heterogeneity and can block the binding site |
| `oxidation_site` | 0 | 1 | **L1** | a solvent-exposed methionine. Oxidises during storage; the M→L swap is routine, so it is a repair opportunity, not a problem |

### One thing deliberately *not* reported

**Tryptophan is not emitted as an oxidation site.** It is oxidation-prone, but
it is also large, aromatic, frequently load-bearing for the fold or for binding
— so there is no cheap substitution — and entirely typical, since every V domain
has several.

A finding that is both *typical* and *not cheaply repairable* is not a finding:
there is no action it could prompt. Emitting it put 98% of clinical-stage
therapeutics into Level 3 in the first calibration run. The fix belongs in the
measurement (do not emit it), not in the matrix (do not route it), so that the
matrix stays monotone in both arguments.

---

## Metrics with no band

If a metric is measured but the loaded band set has no band for it, it appears
in the report under **unmeasured**, and `BandSet.tier()` returns `None` rather
than `0`.

This distinction is load-bearing. Tier 0 means "looked at it, it is typical".
`None` means "nothing looked". Collapsing the second into the first would let a
silent gap in the band file read as a clean bill of health.
