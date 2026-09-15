"""Triage levels: what happens to this candidate, not what it scored.

## Why levels and not a score

The fused score cut at two thresholds could not support the decision it was
being asked to support. Measured on 150 clinical-stage therapeutics: the score
distribution is a single smooth unimodal hill with no gap between GO and
NO-GO, 17% of known-good antibodies sit within +/-3 points of the upper
threshold, and the "10% false rejection" figure was a restatement of where the
p90 line was drawn rather than a measurement.

So this module routes on the same principle the Emergency Severity Index uses:
**a level is defined by the action it triggers.** Level 1 in an emergency
department does not mean "score above X", it means "goes to resuscitation now".

## Numbering follows ESI: 1 is the one you act on first, 5 is the one you do not

This is the direction the Emergency Severity Index uses, and getting it
backwards is easy — an earlier version of this module did. There, Level 1 meant
"reject", which is the opposite of what Level 1 means to anyone who has seen
ESI, where Level 1 is the resuscitation patient. The symptom was visible in the
distribution and went unnoticed: 18 candidates at "Level 1" and zero at
"Level 5", when ESI Level 1 is by nature rare and 4/5 are the bulk.

    LEVEL 1  no repairs needed          -> straight to the bench
    LEVEL 2  framework repairs only     -> one round, no re-measurement
    LEVEL 3  repairs land in a CDR      -> re-measure affinity afterwards
    LEVEL 4  a human must adjudicate    -> serious, but not disqualifying
    LEVEL 5  crossed a boundary         -> do not pursue

Read it as "how much does this cost me, and does it cost me anything at all":
1 costs nothing, 5 costs everything because you never get it back.

## The questions, asked in order

1. Did this cross a boundary no known-good antibody crosses?  -> LEVEL 5
2. Is there a finding a repair count cannot capture?          -> LEVEL 4
3. How many rounds of engineering before it reaches a bench?  -> LEVEL 3/2/1

Only question 1 rejects. The rest order and annotate.

## The repair-cost axis

Levels 3/4/5 count resources the way ESI counts labs and consults. The axis
that matters is NOT how many positions must change — it is:

    **does fixing this risk the binding you already have?**

A framework N-glyc motif and a CDR-H3 N-glyc motif are both one point
mutation. But the CDR fix may destroy the binder, so it costs an extra round
of affinity measurement. That is the real resource difference, and it is
commensurable in a way that summing weighted flags is not.

Conserved IMGT positions are excluded from the count upstream in checks.py:
a "liability" at a structurally required residue is not repairable, so
counting it inflates the estimate with work nobody can do.
"""

from dataclasses import dataclass, field

from .checks import Flag

# Checks eligible to reject a candidate outright (Level 5).
#
# Entry requirement is the trauma-triage standard: measured false-rejection
# rate on known-good antibodies below 5%, ideally with the upper confidence
# bound clearing it too.
#
# Measured on 150 clinical-stage therapeutics (derivation only, not yet
# externally validated — docs/triage_implementation_plan.md Phase 2):
#
#   sequence_sanity      0.0%  (0/150)   LR+ 11.00
#   tap (any RED)        0.0%  (0/150)   LR+ 11.00
#   v_domain_integrity   0.7%  (1/150)   mechanism: no fold without these residues
#   poly_residue_run     see below
#
# poly_residue_run is here on MECHANISM AND PRECEDENT, not on a measured LR.
# Its motifs occur essentially never in natural antibodies (5xW: zero of 300),
# so no natural-antibody population can estimate its LR at any sample size we
# can reach. They are generative-model sampling artefacts, and AbSci's Origin-1
# pipeline filters its Critical tier outright for the same reason.
REJECTING_CHECKS = {
    "sequence_sanity",
    "tap",
    "v_domain_integrity",
    "poly_residue_run",
}

# Several checks emit flags at two severities under one check name, and only
# the severe tier may reject. The tiers are distinguished by weight at the
# point of emission, so the minimum rejecting weight is recorded per check
# rather than inferred from the message text.
#
# This was found by measurement, not by review: the first version put "tap" in
# the rejecting set wholesale, so a TAP AMBER rejected. Result: 46/150 (30.7%) of
# clinical-stage therapeutics auto-rejected, 44 of them on AMBER alone and
# zero on RED. Against a 5% budget. A check that fires on a third of known-good
# candidates cannot gate, and the budget is what catches it.
REJECTING_MIN_WEIGHT = {
    "tap": 6.0,               # tap_runner.RED_WEIGHT; AMBER is 2.0
    "poly_residue_run": 6.0,  # checks.POLY_RUN_CRITICAL_WEIGHT; "Other" is 1.0
}

# Checks worth showing next to a candidate that never change its level on their
# own. Passive display — this is where the alert-fatigue remedy lives.
#
#   n_glycosylation      3.3% on known-good, LR+ 2.64
#   cysteine_pairing     5.3% on known-good, LR+ 2.06  <- over budget, cannot gate
#   disulfide_pairing    structural; LR unmeasured
#   humanness            germline identity; no variance in the human-framework
#                        calibration set, so LR unmeasurable there
#   immunogenicity       LR+ measured at 0.500 — chance — so it informs a human
#                        and never moves a level
REVIEW_CHECKS = {
    "n_glycosylation",
    "cysteine_pairing",
    "disulfide_pairing",
    "humanness",
    "immunogenicity",
    "immunogenicity_observed",
    "poly_residue_run",
}

# Checks that fire on essentially every antibody and therefore cannot
# discriminate, but do estimate work. Measured LR+ = 1.00 for both: they fire
# on 100% of clinical-stage therapeutics AND 100% of patent-text antibodies.
RESOURCE_CHECKS = {"ptm_liability", "developability"}

LEVEL_ACTIONS = {
    1: "no repairs needed — straight to the bench",
    2: "fixable in framework — one round of mutagenesis, no re-measurement",
    3: "fixable, but the fix is in a CDR — re-measure affinity after",
    4: "needs a human to adjudicate before committing bench time",
    5: "do not pursue — crossed a boundary no known-good antibody crosses",
}

# Level 4 is for findings a repair count cannot express: the problem is not
# "how many substitutions" but "somebody has to decide whether this is
# acceptable at all". Deliberately a short list — if every Level 2 check could
# send a candidate here, almost everything would land at 4 and the level would
# mean nothing.
#
#   disulfide_pairing  an unpaired cysteine in the predicted fold is a real
#                      covalent-aggregation risk, and it is also real biology
#                      in some approved molecules. Not repairable by one
#                      substitution, not disqualifying either.
#   humanness          only at the very-low tier, i.e. a framework that is
#                      substantially non-human. Humanising it is a project,
#                      not a repair.
ADJUDICATION_CHECKS = {"disulfide_pairing", "humanness"}
HUMANNESS_ADJUDICATION_MIN_WEIGHT = 5.0  # humanness.VERY_LOW_IDENTITY_WEIGHT


@dataclass
class TriageResult:
    level: int
    action: str
    reasons: list[str] = field(default_factory=list)
    cdr_repairs: int = 0
    framework_repairs: int = 0
    review_flags: list[Flag] = field(default_factory=list)

    # There is deliberately no GO / CONDITIONAL / NO-GO property. Collapsing
    # five levels onto three words destroyed the distinction the levels exist
    # to carry: Level 2 and Level 3 both used to read "CONDITIONAL", but the
    # difference between them — does the fix land in a CDR, so does affinity
    # have to be re-measured — is the entire point of separating them. Worse,
    # framework-repair candidates read as "GO", which invites someone to send
    # them to the bench unrepaired.


def _is_rejecting(flag: Flag) -> bool:
    if flag.check not in REJECTING_CHECKS:
        return False
    minimum = REJECTING_MIN_WEIGHT.get(flag.check)
    return True if minimum is None else flag.weight >= minimum


# Hard gates predate this module and used to reject unconditionally, bypassing
# the false-rejection budget entirely. Measured on the same 300 antibodies:
#
#   non-standard amino acid   0.0% on known-good, 95% upper 2.5%  -> keep
#   CDR-H3 length out of range 0.0% on known-good, 95% upper 2.5% -> keep
#   N-glyc motif inside a CDR  2.7% on known-good, 95% upper 6.7% -> DEMOTED
#
# The third one was rejecting three named clinical-stage therapeutics. It is
# now a weighted, repairable soft flag in checks.py. The surviving two are
# listed here so that "what can reject a candidate" is one visible list rather
# than a severity string scattered across modules.
HARD_GATE_CHECKS = {"sequence_sanity"}


def triage(flags: list[Flag]) -> TriageResult:
    """Assign a triage level from a candidate's flags."""
    hard_gates = [
        f
        for f in flags
        if f.severity == "hard_gate" and f.check in HARD_GATE_CHECKS
    ]
    rejecting = hard_gates + [f for f in flags if _is_rejecting(f)]
    if rejecting:
        return TriageResult(
            level=5,
            action=LEVEL_ACTIONS[5],
            reasons=[f.message for f in rejecting],
        )

    # Repair cost counts only flags explicitly marked repairable — one point
    # mutation fixes it and the rest of the molecule is unchanged. See the
    # comment on Flag.repairable in checks.py for why this is a whitelist and
    # not "everything a resource check emitted".
    repairs = [f for f in flags if f.repairable]
    cdr = sum(1 for f in repairs if f.region.startswith("CDR"))
    framework = len(repairs) - cdr

    adjudication = [f for f in flags if _needs_adjudication(f)]
    review = [f for f in flags if f.check in REVIEW_CHECKS]

    if adjudication:
        level = 4
    elif cdr:
        level = 3
    elif framework:
        level = 2
    else:
        level = 1

    reasons = []
    if adjudication:
        reasons.extend(f.message for f in adjudication[:2])
    if cdr:
        reasons.append(f"{cdr} repair(s) in a CDR — affinity must be re-measured")
    if framework:
        reasons.append(f"{framework} framework repair(s)")
    if not reasons:
        reasons.extend(f.message for f in review[:2])

    return TriageResult(
        level=level,
        action=LEVEL_ACTIONS[level],
        reasons=reasons,
        cdr_repairs=cdr,
        framework_repairs=framework,
        review_flags=review,
    )


def _needs_adjudication(flag: Flag) -> bool:
    if flag.check not in ADJUDICATION_CHECKS:
        return False
    if flag.check == "humanness":
        return flag.weight >= HUMANNESS_ADJUDICATION_MIN_WEIGHT
    return True
