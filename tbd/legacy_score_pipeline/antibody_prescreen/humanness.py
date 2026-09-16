"""Germline identity as a humanness proxy.

## What this is, and what it is NOT

This is **not OASis**. OASis (BioPhi) scores humanness by looking up every
9-mer peptide of a sequence in the Observed Antibody Space database and
reporting what fraction are seen in human repertoires. It is the better
measure, and AbSci used BioPhi for exactly this purpose in Origin-1 lead
optimisation.

It is not implemented here because it needs its own conda environment
(python 3.9 + HMMER + abnumber) and a multi-gigabyte OAS peptide database.
That is a deliberate deferral, recorded in the plan, not an oversight —
and this module must never be described as OASis or as a substitute for it.

What this *is*: the percent identity between the query V domain and its closest
human germline V gene, which ANARCI already computes as a side effect of the
numbering call the pipeline makes anyway. Zero new dependencies, zero new
downloads.

## What it can and cannot tell you

It catches the coarse case — a murine or chimeric framework reads as low
identity — which is the historical HAMA failure mode. It cannot catch a
humanised antibody carrying a handful of non-human back-mutations at positions
that happen to be T-cell epitopes; OASis can, because it works on peptides
rather than on an aggregate percentage.

## An important limitation on calibrating it

The Step 7 population was filtered to human-framework antibodies, so germline
identity has almost no variance there and its likelihood ratio cannot be
estimated on it. It earns its place by mechanism, not by a measured LR on that
population.
"""

from dataclasses import dataclass

from .checks import Flag
from .numbering import NumberedChain

# Clinical-stage humanised and fully human antibodies typically sit above ~80%
# identity to their closest human germline V gene. Below that the framework is
# carrying substantial non-human content.
#
# PROVISIONAL: chosen from the conventional humanisation literature, not fitted
# here — see the calibration limitation in the module docstring.
GERMLINE_IDENTITY_LOW = 0.80
GERMLINE_IDENTITY_VERY_LOW = 0.70

LOW_IDENTITY_WEIGHT = 2.0
VERY_LOW_IDENTITY_WEIGHT = 5.0


@dataclass
class HumannessResult:
    available: bool
    v_gene: str | None = None
    species: str | None = None
    identity: float | None = None
    note: str = ""


def germline_identity(chain: NumberedChain) -> HumannessResult:
    """Identity to the closest HUMAN germline V gene.

    Uses `germline.germline_profile`, whose reference species is the patient's
    (human by default) rather than the antibody's — so a murine chain is
    measured against human germline and correctly reads as low identity,
    instead of scoring high against its own species.
    """
    from .germline import REFERENCE_SPECIES, germline_profile

    profile = germline_profile(chain, reference_species=REFERENCE_SPECIES)
    if not profile.available or profile.identity is None:
        return HumannessResult(
            available=False,
            note="no germline assignment — humanness not scored for this chain",
        )
    return HumannessResult(
        available=True,
        v_gene=profile.v_gene,
        species=profile.species,
        identity=float(profile.identity),
    )


def humanness_flags(result: HumannessResult, chain_name: str) -> list[Flag]:
    if not result.available or result.identity is None:
        return []
    if result.identity >= GERMLINE_IDENTITY_LOW:
        return []

    very_low = result.identity < GERMLINE_IDENTITY_VERY_LOW
    return [
        Flag(
            check="humanness",
            severity="soft",
            region="chain",
            weight=VERY_LOW_IDENTITY_WEIGHT if very_low else LOW_IDENTITY_WEIGHT,
            message=(
                f"{chain_name} is {100 * result.identity:.1f}% identical to its closest "
                f"human germline V gene ({result.v_gene}) — below the "
                f"{100 * GERMLINE_IDENTITY_LOW:.0f}% humanised range. "
                "Germline identity only; NOT an OASis score"
            ),
        )
    ]
