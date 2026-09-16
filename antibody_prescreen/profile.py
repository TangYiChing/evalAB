"""What a candidate IS, before anyone decides whether that is good or bad.

## The question this module refuses to answer

"What is wrong with this antibody" is unanswerable here, and saying so plainly
is the foundation of the whole design. No measurement in this repository can
identify a bad antibody, because there is no wet-lab outcome label anywhere in
it. Building a good/bad classifier on the only labels that do exist — the
provenance contrast "reached the clinic" vs "was merely patented" — measures
an AUC of 0.587 (docs/range-triage-design.md). That is statistically real and
practically nothing, and it is the number that rules out scoring.

So this module answers the question that IS answerable from a reference
population alone:

    Where does this candidate sit relative to antibodies we have seen?

That is descriptive. It needs a reference population and nothing else — no
outcome data, no expert judgement, no weights. `bands.py` supplies the
reference population; this module supplies the coordinates.

## The two kinds of observation, and why they must not be mixed

**MEASURED** (`RangeMetric`) — a number with a distribution. CDR lengths,
isoelectric point, germline identity, aggregation propensity. For these,
"outside the range" is defined by where the reference population actually
sits, and the cut point is a percentile, not an opinion.

**OBSERVED** (`BinaryObservation`) — a mechanistic yes/no with no distribution
worth fitting. A conserved IMGT anchor is present or it is not; an amino acid
is standard or it is not. Fitting a percentile to these would be theatre: the
answer does not depend on how often it happens, it depends on whether the
immunoglobulin fold can survive it.

Holding both in one list, each carrying a weight, is precisely the weighted-sum
thinking this pipeline exists to avoid: it lets a frequency and a mechanism be
traded off against each other as though they were the same currency. They are
kept in separate lists so that they cannot be.

## Repair cost is recorded here, not derived later

Every observation carries the cost of undoing it, on a fixed five-point scale
(see `COST_*`). Cost is a property of the finding, not of the candidate, and it
is knowable at the moment of measurement: a deamidation motif in a framework is
one substitution away from gone, and a CDR-H3 four residues longer than any
clinical antibody is not. Deferring the cost decision to the triage function
loses that locality — the finding is what knows its own repair, so the finding
carries it.
"""

from dataclasses import dataclass, field

from .developability import aggrescan_na4vss, isoelectric_point, net_charge
from .numbering import NumberedChain

# --- Repair cost scale ----------------------------------------------------
#
# The second axis of the triage matrix. Read it as "what a person has to do,
# and what they lose by doing it".
COST_NONE = 0        # nothing to undo
COST_FR_POINT = 1    # one substitution in framework; binding untouched
COST_CDR_POINT = 2   # one substitution in a CDR; affinity must be re-measured
COST_REDESIGN = 3    # a loop, a charge distribution, a patch — not one residue
COST_NO_PATH = 4     # no substitution fixes it; the molecule is not this molecule

COST_NAMES = {
    COST_NONE: "nothing to do",
    COST_FR_POINT: "one framework substitution",
    COST_CDR_POINT: "one CDR substitution + re-measure affinity",
    COST_REDESIGN: "redesign a loop or a surface — not a substitution",
    COST_NO_PATH: "no repair path by substitution",
}


@dataclass
class RangeMetric:
    """A measured number, to be placed against a reference distribution."""
    name: str
    value: float
    cost: int
    # Which tail is a problem. Almost everything here is two-sided: a CDR-H3
    # of 4 residues is as unusual as one of 28, and both are worth seeing.
    # `high` is for counts, where zero is the floor and nobody is alarmed by it.
    tail: str = "both"  # "both" or "high"
    note: str = ""


@dataclass
class BinaryObservation:
    """A mechanistic yes/no. `deviation` is asserted, not fitted — see module docstring."""
    name: str
    deviation: int  # 0..3, on the same scale bands.py produces
    cost: int
    region: str
    message: str


@dataclass
class Profile:
    candidate_id: str
    metrics: list[RangeMetric] = field(default_factory=list)
    observations: list[BinaryObservation] = field(default_factory=list)
    error: str | None = None

    def metric(self, name: str) -> float | None:
        for m in self.metrics:
            if m.name == name:
                return m.value
        return None


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Each anchor carries an explicit cost and deviation rather than a weight.
# 23/41/104 hold the fold together; 118 marks
# the J end, and its usual cause of failure is a truncated paste rather than a
# broken molecule, so it is scored one tier lower.
ANCHORS = {23: ("C", "1st-CYS"), 41: ("W", "CONSERVED-TRP"),
           104: ("C", "2nd-CYS"), 118: ("FW", "J-PHE/J-TRP")}

DEAMIDATION = {("N", "G"), ("N", "S"), ("N", "T"), ("N", "H")}
ISOMERIZATION = {("D", "G"), ("D", "S")}
OXIDATION_PRONE = {"M", "W"}
CONSERVED_POSITIONS = {23, 41, 104, 118}

# AbSci Origin-1 Supplementary Table 7, Critical tier. Left in as an
# observation rather than a metric because a run of five serines has no
# distribution on natural antibodies to fit against — it occurs zero times in
# PLAbDab's therapeutic arm. Its deviation is asserted from the mechanism
# (a generative sampling artefact), which is what BinaryObservation is for.
POLY_RUN_CRITICAL = [("W", 5), ("G", 4), ("Y", 5), ("S", 5)]


def _cdr_residues(chain: NumberedChain):
    return [r for r in chain.residues if r.is_cdr and r.aa != "-"]


def build_profile(candidate_id: str, chains: dict[str, NumberedChain]) -> Profile:
    """Measure a numbered VH/VL pair. No thresholds are applied here."""
    p = Profile(candidate_id=candidate_id)
    vh, vl = chains["VH"], chains["VL"]

    # --- lengths -----------------------------------------------------------
    # Loop length is the one property that genuinely cannot be repaired by
    # substitution: changing it changes the backbone, so cost is REDESIGN even
    # when the excursion is one residue.
    for chain_name, chain in (("h", vh), ("l", vl)):
        for i in (1, 2, 3):
            p.metrics.append(RangeMetric(
                name=f"cdr{chain_name}{i}_len",
                value=float(len(chain.region_sequence(f"CDR{i}"))),
                cost=COST_REDESIGN,
            ))

    # --- germline identity -------------------------------------------------
    # One-sided in spirit (low is the risk) but kept two-sided: a candidate at
    # 100.0% identity to germline across both chains is also unusual, and worth
    # a person seeing, because it usually means the CDRs were not engineered.
    for chain_name, chain in (("vh", vh), ("vl", vl)):
        if chain.v_identity is not None:
            p.metrics.append(RangeMetric(
                name=f"{chain_name}_germline_identity",
                value=float(chain.v_identity),
                cost=COST_REDESIGN,
                note=f"closest human germline {chain.v_gene}",
            ))

    # --- whole-Fv surface properties ---------------------------------------
    fv = vh.full_sequence() + vl.full_sequence()
    p.metrics.append(RangeMetric("fv_pi", isoelectric_point(fv), COST_REDESIGN))
    p.metrics.append(RangeMetric("fv_net_charge", net_charge(fv, 7.4), COST_REDESIGN))

    # CDR-only charge. This is the sequence-level stand-in for TAP's PPC/PNC —
    # the patch metrics need a structure, but the charge that forms them is
    # already visible in the loops.
    cdr_seq = "".join(r.aa for c in (vh, vl) for r in _cdr_residues(c))
    p.metrics.append(RangeMetric("cdr_net_charge", net_charge(cdr_seq, 7.4), COST_REDESIGN))

    for chain_name, chain in (("vh", vh), ("vl", vl)):
        p.metrics.append(RangeMetric(
            f"{chain_name}_aggregation", aggrescan_na4vss(chain.full_sequence()),
            COST_REDESIGN,
            note="AGGRESCAN Na4vSS — a regional property, not one residue",
        ))

    # --- counts ------------------------------------------------------------
    # These ARE repairable per-instance, which is why the same finding appears
    # twice in the profile: once as a count (how unusual is this many?) and once
    # per instance as an observation (what would it cost to remove one?).
    ptm_cdr = ox_cdr = 0
    for chain_name, chain in (("VH", vh), ("VL", vl)):
        seq = chain.full_sequence()
        residues = [r for r in chain.residues if r.aa != "-"]
        for i in range(len(seq) - 1):
            pair = (seq[i], seq[i + 1])
            if pair not in DEAMIDATION and pair not in ISOMERIZATION:
                continue
            r = residues[i]
            kind = "deamidation" if pair in DEAMIDATION else "isomerization"
            in_cdr = r.region.startswith("CDR")
            ptm_cdr += in_cdr
            p.observations.append(BinaryObservation(
                name="ptm_motif", deviation=0,
                cost=COST_CDR_POINT if in_cdr else COST_FR_POINT,
                region=f"{chain_name} {r.region}",
                message=f"{kind} motif '{''.join(pair)}' at {chain_name} {r.region} {r.position}",
            ))
        for r in residues:
            if r.aa not in OXIDATION_PRONE or r.position in CONSERVED_POSITIONS:
                continue
            in_cdr = r.region.startswith("CDR")
            ox_cdr += in_cdr
            # Methionine has a routine M->L swap, so it is a repair
            # opportunity and gets its own observation.
            #
            # Tryptophan does not. It is large, aromatic and frequently
            # load-bearing for the fold or for binding, so there is no cheap
            # substitution. And it is entirely typical — every V domain has
            # several. A finding that is BOTH typical and not cheaply
            # repairable is not a finding at all; there is no action it could
            # prompt. Emitting it anyway is what put 98% of clinical-stage
            # therapeutics into Level 3 on the first calibration run: a
            # tryptophan in FR3 was being read as "redesign needed", on every
            # antibody ever made.
            #
            # Tryptophans still count toward `cdr_oxidation_sites`, which IS
            # banded. That is the right place for them: the question worth
            # asking is not "is there a W in a CDR" (there usually is) but
            # "are there more of them than in any antibody we have seen".
            if r.aa == "M":
                p.observations.append(BinaryObservation(
                    name="oxidation_site", deviation=0,
                    cost=COST_CDR_POINT if in_cdr else COST_FR_POINT,
                    region=f"{chain_name} {r.region}",
                    message=f"oxidation-prone M at {chain_name} {r.region} {r.position}",
                ))

    p.metrics.append(RangeMetric("cdr_ptm_motifs", float(ptm_cdr), COST_CDR_POINT, tail="high"))
    p.metrics.append(RangeMetric("cdr_oxidation_sites", float(ox_cdr), COST_CDR_POINT, tail="high"))
    p.metrics.append(RangeMetric("fv_cys_count", float(fv.count("C")), COST_REDESIGN))

    # --- mechanistic observations -----------------------------------------
    for chain_name, chain in (("VH", vh), ("VL", vl)):
        seq = chain.full_sequence()
        residues = [r for r in chain.residues if r.aa != "-"]

        bad = sorted(set(seq) - STANDARD_AA)
        if bad:
            p.observations.append(BinaryObservation(
                "non_standard_aa", 3, COST_NO_PATH, chain_name,
                f"{chain_name} contains non-standard amino acid(s): {bad}"))

        observed = {r.position: r.aa for r in residues}
        for pos, (expected, label) in ANCHORS.items():
            aa = observed.get(pos)
            deviation, cost = (2, COST_REDESIGN) if pos == 118 else (3, COST_NO_PATH)
            if aa is None:
                p.observations.append(BinaryObservation(
                    "anchor_missing", deviation, cost, chain_name,
                    f"conserved {label} at IMGT {pos} absent in {chain_name} — truncated or incomplete"))
            elif aa not in expected:
                p.observations.append(BinaryObservation(
                    "anchor_substituted", deviation, cost, chain_name,
                    f"conserved {label} at IMGT {pos} is {aa!r}, expected {'/'.join(expected)} ({chain_name})"))

        # --- declared scope -----------------------------------------------
        #
        # This pipeline is scoped to HUMAN antibodies supplied as a paired
        # VH/VL — from an IgG, an scFv, or a paired-chain design. The reference
        # population is filtered to match (see
        # calibration.fit_bands.REFERENCE_DISQUALIFYING and --human-only), so a
        # murine or chimeric candidate is not merely unusual against it, it is
        # being measured against a distribution that was never meant to
        # describe it.
        #
        # Reported as an observation rather than left to `germline_identity`,
        # because that band tops out at deviation tier 1: it is fitted on 393
        # human clinical antibodies, and 393 cannot assert a 1-in-100 boundary.
        # A murine framework would therefore come back "unusual", which is not
        # what is true about it. What is true is that nobody has calibrated
        # anything for it. Level 4 — a person decides — says that; Level 3 does
        # not.
        #
        # This fires on zero antibodies inside the declared scope, so it costs
        # nothing on the reference population by construction.
        germline_species = (chain.v_species or chain.species or "").lower()
        if germline_species and germline_species != "human":
            p.observations.append(BinaryObservation(
                "out_of_scope_species", 2, COST_REDESIGN, chain_name,
                f"{chain_name} closest germline is {germline_species}, not human "
                f"({chain.v_gene}) — outside this pipeline's declared scope; "
                "every band was fitted on human antibodies and does not describe this"))

        # Odd cysteine count. Sequence alone cannot say whether this is a real
        # defect or a genuine non-canonical disulfide, so it is not a rejection.
        # But it is not nothing either: it guarantees at least one cysteine that
        # cannot be canonically paired, which is a free thiol in the expressed
        # protein. Cost is a CDR point mutation because on engineered candidates
        # the unpaired cysteine is almost always one the designer put in a loop.
        if seq.count("C") % 2:
            p.observations.append(BinaryObservation(
                "odd_cysteine", 2, COST_CDR_POINT, chain_name,
                f"{chain_name} has an odd cysteine count ({seq.count('C')}) — "
                "at least one cannot be canonically paired"))

        for i in range(len(seq) - 2):
            if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in "ST":
                r = residues[i]
                in_cdr = r.region.startswith("CDR")
                p.observations.append(BinaryObservation(
                    "n_glycosylation", 0,
                    COST_CDR_POINT if in_cdr else COST_FR_POINT,
                    f"{chain_name} {r.region}",
                    f"N-glycosylation motif '{seq[i:i+3]}' at {chain_name} {r.region} {r.position}"))

        for region in ("CDR1", "CDR2", "CDR3"):
            segment = chain.region_sequence(region)
            for aa, k in POLY_RUN_CRITICAL:
                run = _longest_run(segment, aa)
                if run >= k:
                    p.observations.append(BinaryObservation(
                        "poly_residue_run", 2, COST_CDR_POINT, f"{chain_name} {region}",
                        f"{run} consecutive {aa!r} in {chain_name} {region} "
                        f"(AbSci Origin-1 Critical, threshold {k})"))
                    break

    return p


def _longest_run(segment: str, aa: str) -> int:
    best = run = 0
    for residue in segment:
        run = run + 1 if residue == aa else 0
        best = max(best, run)
    return best
