"""Triage levels: what happens to this candidate, and what it would cost.

## The level is the WORST finding, not the sum of findings

An emergency department triages a patient with a broken finger and chest pain
on the chest pain. It does not average them. Same here:

    level = max(severity of every finding)

A candidate reaches Level 1 only because nothing is worse than Level 1. This
is what makes the level answer a question a person can act on — "if I send
this to the bench, what is the largest cost I could pay?" — and it is also why
there are no weights anywhere in this module. A max needs no weights; only a
sum does.

## Numbering

Level 1 is the candidate you act on first; Level 5 is the one you do not
pursue. Note this is the INVERSE of CTAS/ESI severity, where Level 1 is the
resuscitation patient. The two domains allocate in opposite directions: an
emergency department spends its scarce resource on the sickest, drug discovery
spends its scarce resource (bench time) on the cleanest. You cannot preserve
both mappings, and this one preserves "level number = order you get to it".

| Level | Name | Condition |
|---|---|---|
| 1 | Ready | No findings, or only framework-located repairs. One round of point mutagenesis at most; binding is untouched, so no re-measurement. |
| 2 | Rework | A repairable liability sits in a CDR or CDR vicinity. Fixable, but the fix may take the binding with it — affinity must be re-measured afterwards. |
| 3 | Immune risk | A liability with no repair locus that bears on immunogenicity: a substantially non-human framework, or Fv-wide charge asymmetry. Not a substitution away from fixed; a decision about whether to re-engineer. |
| 4 | Hold | A structural finding that is serious but unresolved — an unpaired cysteine in the predicted fold. Not repairable by substitution, not disqualifying either. A person must decide. |
| 5 | Out of range | Outside the range spanned by every known therapeutic antibody, or not a foldable V domain at all. The flag names which boundary and by how much. |

## Two axes decide which level a check routes to

**Axis 1 — does the finding have a locus?** If it sits on a specific residue
or segment, its position decides its cost, because position decides whether
fixing it risks the binding you already have. Framework -> Level 1, CDR ->
Level 2. PTM motifs, N-glycosylation sites, aggregation hot spots and the four
CDR-region TAP metrics all live here.

**Axis 2 — if it has no locus, does it bear on immunogenicity?** A whole-
molecule property cannot be assigned to a CDR or a framework, and changing it
is a redesign rather than a repair. Those are classified by whether they feed
the anti-drug-antibody pathway: germline identity and Fv-wide charge symmetry
(SFvCSP) do, and route to Level 3.

PTM deliberately stays on Axis 1 even though there is a real mechanistic link
to immunogenicity — deamidation and isomerization create isoaspartate, a
neo-epitope tolerance never covered. It stays because it HAS a locus, and the
locus is what determines its repair cost. Moving it to Axis 2 would discard
the CDR-versus-framework distinction that makes it actionable. The neo-epitope
risk is noted on CDR-located PTM flags instead of relocating the check.

## What Level 3 rests on

Mechanism, not measurement, and the distinction matters.

A check justified by DATA has a measured likelihood ratio: we counted how
often it fires on 150 known-good and 150 background antibodies, so the claim
is falsifiable with other data. A check justified by MECHANISM has a reason it
must matter — no IMGT-104 cysteine means no intradomain disulfide means no
immunoglobulin fold — which does not need statistics but also does not tell
you how often it happens.

Level 3 rests on mechanism. The anti-drug-antibody pathway is well documented,
but the measured likelihood ratio of the MHC-II prediction on this population
is 0.500 — chance. Germline identity could not be measured at all, because the
calibration set was filtered to human frameworks and so has no variance in it.
Anyone with outcome data should expect to revise this level, and should know
that is what they are revising.
"""

from dataclasses import dataclass, field

from .checks import Flag

LEVEL_NAMES = {
    1: "Ready",
    2: "Rework",
    3: "Immune risk",
    4: "Hold",
    5: "Out of range",
}

LEVEL_ACTIONS = {
    1: "send to the bench — at most one round of framework mutagenesis first",
    2: "repairable, but the fix is in a CDR — re-measure affinity afterwards",
    3: "decide whether to re-engineer — this is not a substitution away from fixed",
    4: "a person must look before any bench time is committed",
    5: "do not pursue — outside the range of every known therapeutic antibody",
}

# --- Level 5: out of range ------------------------------------------------
#
# Eligibility is the trauma-triage standard: the UPPER 95% bound of a check's
# false-rejection rate on known-good antibodies must clear 5%.
#
# Measured on 150 clinical-stage therapeutics (derivation only):
#
#   sequence_sanity      0.0% (0/150), 95% upper 2.5%   LR+ 11.00
#   v_domain_integrity   0.7% (1/150), 95% upper 3.7%   mechanism: no fold
#   tap_psh RED          0.0% (0/150), 95% upper 2.5%   LR+  7.00
#   tap_pnc RED          0.0% (0/150), 95% upper 2.5%   LR+  5.00
#   tap_cdr_length RED   0.0% (0/150)   no RED events in either population
#   tap_ppc RED          0.0% (0/150)   no RED events in either population
#   tap_sfvcsp RED       0.0% (0/150)   no RED events in either population
#   poly_residue_run     mechanism and precedent only, see below
#
# The three TAP metrics with no observed RED events stay eligible: a RED means
# "outside the range spanned by every known therapeutic" regardless of which
# axis it is on, and costing nothing on natural antibodies is exactly what a
# gate aimed at generative-model output should do.
#
# poly_residue_run is here on MECHANISM AND PRECEDENT. Its motifs occur
# essentially never in natural antibodies (five consecutive tryptophans: zero
# of 300), so no natural population can estimate its likelihood ratio at any
# sample size reachable from PLAbDab. They are generative-model sampling
# artefacts, and AbSci's Origin-1 pipeline filters its Critical tier outright
# for the same reason.
REJECTING_CHECKS = {
    "sequence_sanity",
    "v_domain_integrity",
    "poly_residue_run",
    "tap_cdr_length",
    "tap_psh",
    "tap_ppc",
    "tap_pnc",
    "tap_sfvcsp",
}

# Checks that emit two severities under one name, where only the severe tier
# rejects. Recorded as a minimum weight because that is set at emission.
#
# Found by measurement, not review: an earlier version listed "tap" wholesale,
# so a TAP AMBER rejected. 46/150 (30.7%) of clinical-stage therapeutics were
# auto-rejected, 44 on AMBER and zero on RED, against a 5% budget.
REJECTING_MIN_WEIGHT = {
    "tap_cdr_length": 6.0,
    "tap_psh": 6.0,
    "tap_ppc": 6.0,
    "tap_pnc": 6.0,
    "tap_sfvcsp": 6.0,
    "poly_residue_run": 6.0,
}

# Hard gates predate this module and used to reject unconditionally, bypassing
# the budget. Measured for the first time on the same 300 antibodies:
#
#   non-standard amino acid     0.0% on known-good, 95% upper 2.5%  keep
#   CDR-H3 length out of range  0.0% on known-good, 95% upper 2.5%  keep
#   N-glyc motif inside a CDR   2.7% on known-good, 95% upper 6.7%  DEMOTED
#
# The third was rejecting Ispectamab, Zanolimumab and Puxitatug, three named
# clinical-stage therapeutics. Its point estimate clears 5% but its upper bound
# does not, and the rule is the upper bound.
HARD_GATE_CHECKS = {"sequence_sanity"}

# --- Level 4: hold --------------------------------------------------------
#
# Deliberately short. If every review check could route here, almost everything
# would, and the level would mean nothing.
HOLD_CHECKS = {"disulfide_pairing"}

# --- Level 3: immune risk -------------------------------------------------
#
# Axis 2: no repair locus, bears on immunogenicity.
IMMUNE_RISK_CHECKS = {"humanness", "tap_sfvcsp"}

# --- Shown but never routing ---------------------------------------------
#
# Displayed next to the candidate; changes no level on its own. This is where
# the alert-fatigue remedy lives — a finding with a likelihood ratio near 1
# that fires on every candidate does not merely waste attention, it trains the
# reader to ignore the ones that matter.
#
#   immunogenicity        LR 0.500 measured — chance
#   cysteine_pairing      5.3% on known-good, over the rejecting budget
#   tap_* at AMBER        LR 1.00-2.43, and not repairable by substitution
#   aggregation, charge   LR 1.00, fire on 100% of both populations
REVIEW_CHECKS = {
    "immunogenicity",
    "immunogenicity_observed",
    "cysteine_pairing",
    "n_glycosylation",
    "aggregation",
    "charge",
    "tap_cdr_length",
    "tap_psh",
    "tap_ppc",
    "tap_pnc",
}


@dataclass
class TriageResult:
    level: int
    name: str
    action: str
    reasons: list[str] = field(default_factory=list)
    cdr_repairs: int = 0
    framework_repairs: int = 0
    review_flags: list[Flag] = field(default_factory=list)

    # There is deliberately no GO / CONDITIONAL / NO-GO property. Collapsing
    # five levels onto three words destroyed the distinction the levels exist
    # to carry: Level 1 and Level 2 would both have read "CONDITIONAL", but
    # whether the fix lands in a CDR — and therefore whether affinity has to be
    # re-measured — is the entire reason they are separate.


def _is_rejecting(flag: Flag) -> bool:
    if flag.check not in REJECTING_CHECKS:
        return False
    minimum = REJECTING_MIN_WEIGHT.get(flag.check)
    return True if minimum is None else flag.weight >= minimum


def _is_immune_risk(flag: Flag) -> bool:
    if flag.check == "humanness":
        return True
    # SFvCSP only lands here at AMBER; a RED is a Level 5 excursion.
    if flag.check == "tap_sfvcsp":
        return flag.weight < REJECTING_MIN_WEIGHT["tap_sfvcsp"]
    return False


def _result(level: int, reasons: list[str], **kwargs) -> "TriageResult":
    return TriageResult(
        level=level,
        name=LEVEL_NAMES[level],
        action=LEVEL_ACTIONS[level],
        reasons=reasons,
        **kwargs,
    )


def triage(flags: list[Flag]) -> TriageResult:
    """Assign a triage level. The level is the worst finding, not their sum."""
    review = [f for f in flags if f.check in REVIEW_CHECKS]

    rejecting = [
        f for f in flags if f.severity == "hard_gate" and f.check in HARD_GATE_CHECKS
    ] + [f for f in flags if _is_rejecting(f)]
    if rejecting:
        return _result(5, [f.message for f in rejecting], review_flags=review)

    holds = [f for f in flags if f.check in HOLD_CHECKS]
    if holds:
        return _result(4, [f.message for f in holds], review_flags=review)

    immune = [f for f in flags if _is_immune_risk(f)]
    if immune:
        return _result(3, [f.message for f in immune[:2]], review_flags=review)

    # Repair cost counts only flags explicitly marked repairable — one point
    # mutation fixes it and the rest of the molecule is unchanged. See the
    # comment on Flag.repairable in checks.py for why this is a whitelist.
    repairs = [f for f in flags if f.repairable]
    cdr = sum(1 for f in repairs if f.region.startswith("CDR"))
    framework = len(repairs) - cdr

    if cdr:
        return _result(
            2,
            [f"{cdr} repair(s) in a CDR — affinity must be re-measured"]
            + ([f"{framework} framework repair(s)"] if framework else []),
            cdr_repairs=cdr,
            framework_repairs=framework,
            review_flags=review,
        )

    reasons = (
        [f"{framework} framework repair(s) — one round, binding untouched"]
        if framework
        else ["no repairs needed"]
    )
    return _result(
        1,
        reasons,
        cdr_repairs=0,
        framework_repairs=framework,
        review_flags=review,
    )
