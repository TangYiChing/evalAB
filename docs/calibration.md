# Calibrating evalAB to your own programme

evalAB ships bands fitted on a public antibody population. Those bands answer
"where does this candidate sit relative to antibodies that exist". They cannot
answer "will this candidate express, purify and stay monomeric in my hands",
because no wet-lab outcome from your programme is in them.

This document is how you add that — and, just as importantly, how to avoid the
two ways this kind of calibration usually goes wrong.

---

## 1. Refitting the bands to a different population

```bash
# Fetch the public reference (~11 MB, checksummed in data/reference/README.md)
python -c "from antibody_prescreen.calibration import fetch_paired_sequences as f; print(f())"

# Refit. --clean and --human-only are reference-population filters only;
# they never drop a candidate at screening time.
python -m antibody_prescreen.calibration.fit_bands \
    --reference plabdab --human-only --clean \
    --clinical-bands antibody_prescreen/data/bands_clinical.json \
    --out antibody_prescreen/data/bands_human_repertoire.json
```

The script splits the reference set deterministically in half. **Cut points
come from the derivation half only; every rate it prints is measured on the
validation half.** A rejection rate measured on the same antibodies the cut
points were fitted to is a restatement of where you drew the lines, not a
result.

Point `--bands` at your file to use it, or pass `band_set=` to `screen_batch`.
Every report prints the reference, split sizes and cut points in its header,
so a result can always be traced to the population that produced it.

### Which population, and what that choice imports

| Reference | Means | Inherits |
|---|---|---|
| `plabdab` (all pairings) | "the range of things that are antibodies" | no development filter at all; includes murine and chimeric sequences |
| `therasabdab` | "the range of things that survived development" | survivorship — these were made, purified and dosed |

Neither is correct in the abstract; they answer different questions. Germline
identity is the one metric that stays anchored to the clinical population,
because a humanness band fitted on a corpus containing murine sequences widens
exactly far enough to stop catching them — it would calibrate away its own
purpose.

---

## 2. Campaign-specific calibration against your own results

This is the part that makes evalAB yours, and the part with the sharpest
failure mode.

**Pre-register before you look.** Write down, before the batch runs:

- the endpoints — purity, expression yield (mg/L), SEC monomer %, Tm, and
  whatever else you actually measure;
- the thresholds that count as success on each;
- a holdout set that is **independent of this batch** and stays sealed.

Then run. Then compare.

**The failure mode this prevents:** adjusting thresholds until the numbers look
good on the same batch that produced them, and reporting the result as a
validation. It is not one. It is the fit. If you moved a threshold after seeing
the outcomes, the only honest next step is a fresh batch, sealed in advance.

Keep a separate holdout per batch. One holdout reused across campaigns is
consumed the first time you look at it.

---

## 3. Ranking within a level — the part evalAB deliberately does not do

The levels are a triage, not a ranking. Within L1–L3, evalAB has nothing
further to say, and it will not invent a tiebreak by re-weighting sequence
flags into a composite score. That composite is exactly the thing this design
rejects; reintroducing it at the last step would undo the point.

Rank inside a level on data that is actually about your project:

- affinity
- expression yield
- SEC monomer fraction
- Tm / thermal stability
- polyreactivity
- activity on the target cells

A sequence flag is a statement about deviation and repair cost. It is not a
proxy for any of the above, and averaging several flags together does not make
it one.

---

## 4. Generative designs need their own reference

The bands shipped here are fitted on natural, patented and clinical
antibodies. In that population, a poly-residue run or an odd cysteine count is
rare, and its rarity is informative because it reflects what B-cell repertoires
and development pipelines produce.

A de novo generative model has no such process behind it. The same observation
in a generated sequence may be an artefact of the sampler, a benign quirk, or a
genuine liability — the natural-population frequency does not tell you which.
Transferring the interpretation across is an assumption, not a measurement.

If you screen generative output at scale, fit a reference population **from
your own generator's sequences** and compare a design against its own kind.
`fit_bands` takes any population; the reference is a parameter, not a fact.

---

## 5. Tier 2 (TAP) is an upgrade check, not a replacement

`--run-structure` builds an ABodyBuilder2 model and runs TAP. It is off by
default and deliberately so:

- it needs a structural model, which costs minutes per candidate, not
  milliseconds;
- the range-based Tier 1 workflow has no independent outcome validation of
  TAP's contribution to *it* — TAP's own thresholds were drawn around the
  distribution of known therapeutics, so measuring them on known therapeutics
  measures the same ruler twice;
- consequently, a Tier 2 result is reported as an **annotation and never
  changes a Tier 1 level**. `StructureAnnotation.available=False` records a
  failure rather than silently degrading the batch.

Use it as a closer look at a short list you already care about. Do not use it
as a cheaper substitute for the wet-lab endpoints in §2.
