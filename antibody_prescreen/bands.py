"""What "outside the common range" means, and who gets to decide.

## The circularity this module exists to avoid

TAP reports AMBER and RED. Where did those come from? They were drawn around
the distribution of known therapeutic antibodies. So when this repository
measured TAP's likelihood ratio on a population of known therapeutic
antibodies, it measured how much that distribution resembles itself, and got
1.00 to 2.43. That is not a weak signal. That is the same ruler applied twice.

The fix is not a better statistic. It is to stop inheriting other people's
lines and draw our own, from a stated population, with the population written
down next to every threshold. A band produced by this module always carries
its provenance: which reference set, how many antibodies, which half of it.

## What a band is

For one measured quantity, four percentile cut points on the reference
population, and the deviation tier they define:

    tier 0  inside p5-p95      typical; this is what antibodies look like
    tier 1  p95-p99            unusual — 1 in 20 known antibodies is out here
    tier 2  p99-p99.9          strongly unusual — 1 in 100
    tier 3  beyond p99.9       outside the range the reference population spans

Two-sided by default, because a CDR-H3 of 4 residues is as unusual as one of
28 and a person should see both. Counts (`tail="high"`) are one-sided: nobody
needs telling that zero deamidation motifs is unusually few.

## Derivation and validation halves, and why both exist

Fitting the cut points on a population and then reporting the rejection rate
on the same population reports where the lines were drawn, not what they do.
That was measurement (b) in the triage plan's critique of the old thresholds
and it is still true here. So `fit_bands` splits the reference set
deterministically in half: cut points come from the derivation half only, and
every rate quoted in the calibration report is measured on the validation half,
which the cut points have never seen.

## The reference population is a parameter, not a fact

The default is TheraSAbDab — antibodies that reached the clinic. That choice
imports a selection: these molecules were made, purified and dosed, so "common
range" here means "the range of things that survived development".

The alternative is the full PLAbDab table (~170k paired sequences, mostly
patent text), which means "the range of things that are antibodies" and carries
no development filter at all. Neither is right in the abstract. They answer
different questions, and a band file records which one it answered.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path

# The percentile cut points. Changing these changes what "unusual" means, so
# they live in one place with the reasoning attached.
#
# p95/p99/p99.9 rather than something tighter: the union budget is what has to
# clear 5%, not each metric (the lesson from the D3 external validation, where
# eight gates at ~1% each summed to 5.53%). With ~17 metrics, a tier-3 rate of
# 0.2% each has room to be independent and still land near 3%. Whether it
# actually does is measured, not assumed — see the calibration report.
CUTS = (5.0, 95.0, 99.0, 99.9)

# How many observations a tail must actually contain before a cut point placed
# in it is allowed to mean anything.
#
# This is the rule that stopped the first design run from being nonsense. Every
# one of 23 designed candidates came back Level 5 on `cdrh1_len=7`, because the
# TheraSAbDab reference has p0.1 = p1 = p5 = 8 for that metric: 534 clinical
# antibodies, every single one with an IMGT CDR-H1 of 8 or more. The p99.9 cut
# was therefore the minimum of the sample wearing a percentile's name, and
# "below p0.1" was being read as "beyond a 1-in-1000 boundary" when the sample
# cannot see 1 in 1000 at all. 534 x 0.001 = 0.53 antibodies.
#
# So a tier is assertable only when its tail would hold at least this many
# reference antibodies. Anything finer is reported at the finest tier the
# population CAN support, and the band says so.
#
# The consequence is worth stating plainly because it constrains the whole
# system: with a reference of ~1,000 antibodies, NO range metric can reach
# deviation tier 3. Tier 3 becomes reachable only from mechanism — a
# non-standard amino acid, a broken fold anchor — which is the honest position.
# Wanting tier 3 from a distribution means wanting a reference population of
# ~5,000 or more, which is a reason to prefer full PLAbDab over TheraSAbDab,
# not a reason to pretend 534 is enough.
MIN_TAIL_OBS = 5
TIER_TAIL_FRACTION = {1: 5.0, 2: 1.0, 3: 0.1}


def supported_tier(n: int) -> int:
    """Finest deviation tier a reference population of size `n` can assert."""
    best = 0
    for tier, frac in sorted(TIER_TAIL_FRACTION.items()):
        if n * frac / 100.0 >= MIN_TAIL_OBS:
            best = tier
    return best

TIER_NAMES = {
    0: "typical",
    1: "unusual",
    2: "strongly unusual",
    3: "outside the reference range",
}


@dataclass
class Band:
    metric: str
    tail: str
    # Six numbers: p0.1, p1, p5, p95, p99, p99.9 of the reference population.
    low_extreme: float
    low_strong: float
    low_mild: float
    high_mild: float
    high_strong: float
    high_extreme: float
    n: int
    # The finest tier this band's sample size can assert. See MIN_TAIL_OBS.
    max_tier: int = 3
    # Which reference population THIS band came from. Per-band, not per-file,
    # because a band set can legitimately mix them — see `compose`.
    source: str = ""

    def raw_tier(self, value: float) -> int:
        """Which cut point `value` crosses, before the resolution cap."""
        if value > self.high_extreme or (self.tail == "both" and value < self.low_extreme):
            return 3
        if value > self.high_strong or (self.tail == "both" and value < self.low_strong):
            return 2
        if value > self.high_mild or (self.tail == "both" and value < self.low_mild):
            return 1
        return 0

    def tier(self, value: float) -> int:
        """Deviation tier, capped at what this population can support."""
        return min(self.raw_tier(value), self.max_tier)

    def degenerate_side(self, value: float) -> bool:
        """True when the tail this value sits in has collapsed to one number.

        Discrete metrics do this constantly — CDR-L2 is 3 residues in
        essentially every antibody, so p5 and p95 are both 3 and any other
        value is "outside the band". That is real information, but it is not a
        percentile, and the message must not pretend otherwise.
        """
        if value > self.high_mild:
            return self.high_mild == self.high_extreme
        return self.low_mild == self.low_extreme

    def describe(self, value: float) -> str:
        side = "above" if value > self.high_mild else "below"
        raw, capped = self.raw_tier(value), self.tier(value)
        text = (
            f"{self.metric}={value:g} — {side} the "
            f"p{CUTS[0]:g}-p{CUTS[1]:g} range [{self.low_mild:g}, {self.high_mild:g}] "
            f"of {self.source or 'the reference population'} (n={self.n})"
        )
        if self.degenerate_side(value):
            text += (" [that tail is a single value in the reference — "
                     "read this as 'never seen', not as a percentile]")
        if raw > capped:
            text += (f" [crossed the p{CUTS[raw]:g} cut, but n={self.n} cannot "
                     f"resolve past tier {self.max_tier}]")
        return text


@dataclass
class BandSet:
    reference: str            # e.g. "PLAbDab/TheraSAbDab"
    n_derivation: int
    n_validation: int
    cuts: tuple
    bands: dict = field(default_factory=dict)
    notes: str = ""

    def tier(self, metric: str, value: float) -> int | None:
        """None means "no band for this metric" — which must NOT be read as 0.

        A metric with no reference distribution has not been shown to be
        typical; it has not been looked at. Callers surface it as unmeasured.
        """
        band = self.bands.get(metric)
        return None if band is None else band.tier(value)

    def band(self, metric: str) -> Band | None:
        return self.bands.get(metric)

    def to_json(self) -> str:
        return json.dumps(
            {
                "reference": self.reference,
                "n_derivation": self.n_derivation,
                "n_validation": self.n_validation,
                "cuts": list(self.cuts),
                "notes": self.notes,
                "bands": {k: vars(v) for k, v in self.bands.items()},
            },
            indent=2,
        )

    @classmethod
    def load(cls, path: Path | str) -> "BandSet":
        d = json.loads(Path(path).read_text())
        return cls(
            reference=d["reference"],
            n_derivation=d["n_derivation"],
            n_validation=d["n_validation"],
            cuts=tuple(d["cuts"]),
            bands={k: Band(**v) for k, v in d["bands"].items()},
            notes=d.get("notes", ""),
        )


def percentile(sorted_values: list[float], p: float) -> float:
    """Linear-interpolated percentile. No numpy dependency in the Tier 1 path."""
    if not sorted_values:
        raise ValueError("empty")
    if len(sorted_values) == 1:
        return sorted_values[0]
    k = (len(sorted_values) - 1) * p / 100.0
    lo = int(k)
    hi = min(lo + 1, len(sorted_values) - 1)
    return sorted_values[lo] + (sorted_values[hi] - sorted_values[lo]) * (k - lo)


def fit_bands(
    observations: dict[str, list[float]],
    tails: dict[str, str],
    reference: str,
    n_derivation: int,
    n_validation: int,
    notes: str = "",
) -> BandSet:
    """Fit one band per metric from DERIVATION-half values only.

    `observations` must already be the derivation half. This function does not
    split anything — the split belongs to the caller, which is also what writes
    down how it split, so the two cannot drift apart.
    """
    bands = {}
    for metric, values in observations.items():
        values = sorted(v for v in values if v is not None)
        # A p99.9 fitted on fewer than ~200 observations is the maximum of the
        # sample wearing a percentile's name. Refusing is better than pretending.
        if len(values) < 200:
            continue
        bands[metric] = Band(
            metric=metric,
            tail=tails.get(metric, "both"),
            low_extreme=percentile(values, 0.1),
            low_strong=percentile(values, 1.0),
            low_mild=percentile(values, CUTS[0]),
            high_mild=percentile(values, CUTS[1]),
            high_strong=percentile(values, CUTS[2]),
            high_extreme=percentile(values, CUTS[3]),
            n=len(values),
            max_tier=supported_tier(len(values)),
            source=reference,
        )
    return BandSet(
        reference=reference,
        n_derivation=n_derivation,
        n_validation=n_validation,
        cuts=CUTS,
        bands=bands,
        notes=notes,
    )


# --- Mixing reference populations ----------------------------------------
#
# Not every metric belongs to the same population, and pretending otherwise was
# a real error rather than a simplification.
#
# The base population is the cleaned human PLAbDab set (see
# `calibration.fit_bands.clean_reference`). It is the right default for almost
# everything, for two reasons that happen to point the same way:
#
#   (a) Loop lengths and cysteine count are properties of the repertoire, not
#       of what survived a development pipeline. Changing a CDR length is a
#       routine engineering operation, and framing it as an excursion against
#       534 clinical antibodies asks "what do clinical drugs look like" and
#       then reads the answer as "what do antibodies look like". CDR-H3's
#       p99.9 is 24.4 clinically and 63 in the repertoire.
#
#   (b) Surface properties - pI, net charge, aggregation propensity - turned
#       out to be nearly IDENTICAL in the two populations. `vh_aggregation`
#       has a p5 of 0.154 in both, to three significant figures. So the
#       population choice is nearly free for these, and the larger one wins on
#       resolution alone: 534 antibodies cannot assert a 1-in-1000 boundary and
#       ~30,000 can. That these bands barely move is itself a finding - on
#       sequence-level measures, development does not measurably narrow the
#       distribution relative to the patent corpus.
#
# Germline identity is the exception, and it is not a close call. PLAbDab's
# patent corpus contains murine and chimeric sequences. A humanness band fitted
# on a population that includes them widens precisely far enough to stop
# catching them - the check would calibrate away its own purpose. Even after
# the human-only filter, "human" is ANARCI's germline call rather than a
# provenance record, and a humanised antibody with murine back-mutations passes
# it. So this one stays anchored to the clinical population, where "human
# enough to dose a person" is what the distribution actually encodes.
CLINICAL_METRICS = (
    "vh_germline_identity",
    "vl_germline_identity",
)


def compose(base: BandSet, override: BandSet, metrics) -> BandSet:
    """Take `metrics` from `override`, everything else from `base`."""
    bands = dict(base.bands)
    taken = []
    for m in metrics:
        if m in override.bands:
            bands[m] = override.bands[m]
            taken.append(m)
    return BandSet(
        reference=f"{base.reference} + {override.reference} for {', '.join(taken)}",
        n_derivation=base.n_derivation,
        n_validation=base.n_validation,
        cuts=base.cuts,
        bands=bands,
        notes=(base.notes + f" Metrics {taken} were re-fitted on "
               f"{override.reference} (n={override.n_derivation}); see bands.compose."),
    )
