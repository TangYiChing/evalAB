"""The triage matrix: deviation x repair cost -> level.

## Routing on properties, never on check names

The tempting design routes on the identity of the check that fired — this
check means escalate, that one means hold — but every such set is a list of
separate decisions, each of which has to be argued on its own, and none of
which generalises to the next check somebody adds.

This module routes instead on two properties that every finding already
carries:

    how far outside the reference range is it   (deviation, 0-3, from bands.py)
    what does undoing it cost                   (cost, 0-4, from profile.py)

and one table. No check name appears in the routing logic anywhere in this
module. Adding a check therefore requires no routing argument — only answers
to the same two questions every existing check already answers.

## The table

                         cost ->
                   0     1       2        3         4
                 none  FR pt  CDR pt  redesign  no path
    deviation 0    1      1      2        3         4
              1    1      2      2        3         4
              2    2      2      3        4         5
              3    3      3      4        5         5

Read a cell as a sentence. `deviation=3, cost=4` is "outside the range any
known antibody spans, and no substitution brings it back" — Level 5, do not
pursue. `deviation=0, cost=2` is "a perfectly ordinary deamidation motif that
happens to sit in a CDR" — Level 2, fix it and re-measure affinity.

The two corners are the parts worth defending:

**deviation=0, cost>=3 must not be reachable, and that is a design claim, not
a table entry.** A finding that is both entirely typical and not cheaply
repairable is not a finding: there is no action it could prompt, so reporting
it can only produce noise. The first calibration run proved the point by
violating it — an oxidation-prone tryptophan in FR3 was emitted as a
deviation-0, redesign-cost observation, and put 98% of clinical-stage
therapeutics into Level 3. The fix belongs in `profile.py` (do not emit it),
not in the table (do not route it), because the table should stay monotone in
both arguments. See the note on tryptophan there.

**deviation=3, cost<=1 is Level 3, not Level 5.** A candidate whose pI sits
beyond anything in the reference population, but which one framework
substitution brings back, is not a dead candidate. It is a candidate with a
known price. Rejecting it outright is what a pass/fail gate does, and it is
the specific behaviour this matrix exists to prevent: extreme deviation and
cheap repair is an ordinary combination, not a verdict.

## The level is still a max, and there are still no weights

`level = max(level(finding) for finding in findings)`. Two Level-2 findings do
not make a Level 3, and that is deliberate: a person needs to know the largest
price they could end up paying, and a max answers that question while a sum
answers a different one nobody asked. Summing would also smuggle weights back
in through the choice of what counts as one finding.

## A level is a position, not a category

No level is the home of a particular check. Immunogenicity, aggregation and
germline identity are checks like any other, and where each lands depends
entirely on how far outside the range it sits and what repairing it costs. The
same check can produce a Level 1 on one candidate and a Level 4 on the next.
That is what it means for the routing to be a property of the finding rather
than of the check's name.
"""

from dataclasses import dataclass, field

from .bands import TIER_NAMES, BandSet
from .profile import COST_NAMES, Profile

MATRIX = {
    0: {0: 1, 1: 1, 2: 2, 3: 3, 4: 4},
    1: {0: 1, 1: 2, 2: 2, 3: 3, 4: 4},
    2: {0: 2, 1: 2, 2: 3, 3: 4, 4: 5},
    3: {0: 3, 1: 3, 2: 4, 3: 5, 4: 5},
}

LEVEL_NAMES = {
    1: "Ready",
    2: "Point repair",
    3: "Redesign",
    4: "Human decision",
    5: "Out of range",
}

LEVEL_ACTIONS = {
    1: "put it on the wet-lab list as it stands",
    2: "one or two substitutions first; re-measure affinity if any is in a CDR",
    3: "a loop or a surface has to be redesigned — budget a design cycle, not a mutation",
    4: "nothing routes this automatically; a person decides whether it is worth the bench time",
    5: "outside the reference range with no repair path — the flag names which boundary",
}


@dataclass
class Finding:
    source: str
    deviation: int
    cost: int
    level: int
    message: str


@dataclass
class Triage:
    candidate_id: str
    level: int
    name: str
    action: str
    findings: list[Finding] = field(default_factory=list)
    unmeasured: list[str] = field(default_factory=list)
    error: str | None = None

    @property
    def drivers(self) -> list[Finding]:
        """Only the findings that actually set the level. Everything else is context."""
        return [f for f in self.findings if f.level == self.level]

    def repair_counts(self) -> dict:
        c = {"framework_point": 0, "cdr_point": 0, "redesign": 0, "no_path": 0}
        key = {1: "framework_point", 2: "cdr_point", 3: "redesign", 4: "no_path"}
        for f in self.findings:
            if f.cost in key:
                c[key[f.cost]] += 1
        return c


def assign(profile: Profile, band_set: BandSet) -> Triage:
    """Place one profile on the matrix. Nothing here knows any check's name."""
    if profile.error:
        return Triage(profile.candidate_id, 0, "ERROR", profile.error, error=profile.error)

    findings: list[Finding] = []
    unmeasured: list[str] = []

    for m in profile.metrics:
        tier = band_set.tier(m.name, m.value)
        if tier is None:
            unmeasured.append(f"{m.name}={m.value:g} (no reference band)")
            continue
        if tier == 0:
            continue
        band = band_set.band(m.name)
        findings.append(Finding(
            source=m.name, deviation=tier, cost=m.cost,
            level=MATRIX[tier][m.cost],
            message=f"{TIER_NAMES[tier]}: {band.describe(m.value)}"
                    + (f" [{m.note}]" if m.note else ""),
        ))

    for o in profile.observations:
        findings.append(Finding(
            source=o.name, deviation=o.deviation, cost=o.cost,
            level=MATRIX[o.deviation][o.cost],
            message=o.message + f" — {COST_NAMES[o.cost]}",
        ))

    level = max((f.level for f in findings), default=1)
    findings.sort(key=lambda f: (-f.level, -f.deviation, -f.cost))
    return Triage(
        candidate_id=profile.candidate_id,
        level=level,
        name=LEVEL_NAMES[level],
        action=LEVEL_ACTIONS[level],
        findings=findings,
        unmeasured=unmeasured,
    )
