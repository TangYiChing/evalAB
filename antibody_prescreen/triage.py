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

## The three questions, asked in order

1. Did this cross a boundary no known-good antibody crosses? -> LEVEL 1
2. Does anything here need a human to look at it?            -> LEVEL 2
3. How many rounds of engineering before it reaches a bench? -> LEVEL 3/4/5

Only question 1 rejects. Questions 2 and 3 order and annotate.

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

# Checks eligible to reject a candidate outright.
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
LEVEL_1_CHECKS = {
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
# LEVEL_1_CHECKS wholesale, so a TAP AMBER rejected. Result: 46/150 (30.7%) of
# clinical-stage therapeutics auto-rejected, 44 of them on AMBER alone and
# zero on RED. Against a 5% budget. A check that fires on a third of known-good
# candidates cannot gate, and the budget is what catches it.
REJECTING_MIN_WEIGHT = {
    "tap": 6.0,               # tap_runner.RED_WEIGHT; AMBER is 2.0
    "poly_residue_run": 6.0,  # checks.POLY_RUN_CRITICAL_WEIGHT; "Other" is 1.0
}

# Checks worth a human's attention that may not reject on their own.
#
#   n_glycosylation      3.3% on known-good, LR+ 2.64
#   cysteine_pairing     5.3% on known-good, LR+ 2.06  <- over budget, cannot gate
#   disulfide_pairing    structural; LR unmeasured
#   humanness            germline identity; no variance in the human-framework
#                        calibration set, so LR unmeasurable there
#   immunogenicity       LR+ measured at 0.500 — chance — so it informs a human
#                        and never moves a level
LEVEL_2_CHECKS = {
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
    1: "do not pursue — crossed a boundary no known-good antibody crosses",
    2: "needs a human to look before committing bench time",
    3: "fixable, but the fix is in a CDR — re-measure affinity after",
    4: "fixable in framework — one round of mutagenesis, no re-measurement",
    5: "no repairs needed — straight to the bench",
}


@dataclass
class TriageResult:
    level: int
    action: str
    reasons: list[str] = field(default_factory=list)
    cdr_repairs: int = 0
    framework_repairs: int = 0
    level2_flags: list[Flag] = field(default_factory=list)

    @property
    def verdict(self) -> str:
        """Backwards-compatible mapping onto the old vocabulary."""
        return {1: "NO-GO", 2: "CONDITIONAL", 3: "CONDITIONAL"}.get(self.level, "GO")


def _is_rejecting(flag: Flag) -> bool:
    if flag.check not in LEVEL_1_CHECKS:
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
    if hard_gates:
        return TriageResult(
            level=1,
            action=LEVEL_ACTIONS[1],
            reasons=[f.message for f in hard_gates],
        )

    rejecting = [f for f in flags if _is_rejecting(f)]
    if rejecting:
        return TriageResult(
            level=1,
            action=LEVEL_ACTIONS[1],
            reasons=[f.message for f in rejecting],
        )

    # Repair cost counts only flags explicitly marked repairable — one point
    # mutation fixes it and the rest of the molecule is unchanged. See the
    # comment on Flag.repairable in checks.py for why this is a whitelist and
    # not "everything a resource check emitted".
    repairs = [f for f in flags if f.repairable]
    cdr = sum(1 for f in repairs if f.region.startswith("CDR"))
    framework = len(repairs) - cdr

    level2 = [f for f in flags if f.check in LEVEL_2_CHECKS]

    if cdr:
        level = 3
    elif framework:
        level = 4
    else:
        level = 5

    # A Level 2 finding cannot make a candidate worse than "needs a look", but
    # it must not let one be waved through as Level 5 either.
    if level == 5 and level2:
        level = 4

    reasons = []
    if cdr:
        reasons.append(f"{cdr} repair(s) in a CDR — affinity must be re-measured")
    if framework:
        reasons.append(f"{framework} framework repair(s)")
    reasons.extend(f.message for f in level2[:2])

    return TriageResult(
        level=level,
        action=LEVEL_ACTIONS[level],
        reasons=reasons,
        cdr_repairs=cdr,
        framework_repairs=framework,
        level2_flags=level2,
    )
