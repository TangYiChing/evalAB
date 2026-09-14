"""Immunogenicity risk: will this antibody be seen as foreign and attacked?

## What this check is actually asking

An antibody drug fails immunogenically like this:

    an antigen-presenting cell engulfs the drug
    -> chops it into ~15-mers
    -> loads them onto MHC-II and presents them at the surface
    -> a CD4+ helper T cell recognises one
    -> B cells are recruited and make anti-drug antibodies (ADA)
    -> the drug is neutralised and cleared

This belongs next to the developability checks because it is another way the
molecule dies before it helps anyone — and the two are mechanistically linked,
not merely adjacent: aggregates are taken up by APCs far more efficiently than
monomer, so a high Tier 2 hydrophobic patch score (`structure/tap_runner.py`)
feeds directly into step 1 above.

## The trap, and why this module looks the way it does

MHC-II binding is a *necessary* condition, not a sufficient one. Between
presentation and a T-cell response sits **central tolerance**: T cells that
react to self peptides are deleted in the thymus. A germline-encoded framework
peptide can be presented beautifully and provoke nothing, because no T cell
is left to see it.

IEDB cannot supply that. Its prediction API covers MHC-II binding; its query
API records observed T-cell assay outcomes. Neither knows what is self. So the
germline correction is ours to make, and without it this check is worthless:

- Measured on the reference VH: 8 strong binders (rank < 2.0), **0 in CDRs,
  8 in framework**. The strongest carries `YLQMNSLRAEDT` — present in 24.6%
  of PLAbDab heavy chains and 34.1% of clinical-stage therapeutics.
- The filter must work on the **9-mer binding core**, not the 15-mer window.
  Measured: window-level filtering kept 8/8 of those false positives (each
  window happened to contain one mutated flanking residue); core-level
  filtering dropped 8/8. All eight shared the identical, fully germline core
  `YLQMNSLRA`.
- Positive control: grafting influenza HA306-318 into CDR-H3 gives 15 strong
  binders — the 8 germline framework ones are dropped and all 7 windows of the
  grafted epitope (core `YVKQNTLKL`) are kept.

The same correction is applied before the evidence layer, not just the
prediction layer: `YLQMNSLRA` has 16 positive human T-cell assay records in
IEDB. The experimental layer false-positives on germline exactly as the
predictor does.

See docs/iedb-immunogenicity-evaluation.md for the full measurements.

## Network

Two IEDB services, both optional and both degrading gracefully:

- prediction: `tools-cluster-interface.iedb.org` (**https** — the http URL
  308-redirects and urllib will not re-POST, which is what previously made
  this check report itself as unavailable)
- evidence: `query-api.iedb.org` (IQ-API)

Responses are cached on disk, because the input is immutable and calibration
re-runs the same batch repeatedly.
"""

import hashlib
import json
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

from .checks import Flag
from .germline import germline_profile
from .numbering import NumberedChain

IEDB_MHCII_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
IEDB_IQ_TCELL_URL = "https://query-api.iedb.org/tcell_search"

DEFAULT_CACHE_DIR = Path.home() / ".cache" / "evalab" / "iedb"

# A small default panel of common human HLA-DR alleles. This is a starting
# point, not population coverage — it is a real parameter of the check and
# should be chosen deliberately for the population a drug is aimed at.
DEFAULT_ALLELES = [
    "HLA-DRB1*01:01",
    "HLA-DRB1*03:01",
    "HLA-DRB1*04:01",
    "HLA-DRB1*07:01",
    "HLA-DRB1*15:01",
]

METHOD = "netmhciipan_el"
WINDOW_SIZE = 15

# netmhciipan reports a percentile rank; lower = stronger predicted binder.
# < 2.0 is the field's conventional "likely binder" cutoff.
BINDER_PERCENTILE_CUTOFF = 2.0

# Weights for a NON-GERMLINE predicted binder. Germline-core binders are
# dropped entirely before weighting, so these only ever apply to peptides
# tolerance does not cover.
#
# PROVISIONAL, like the TAP weights — not fitted against outcomes.
CDR_BINDER_WEIGHT = 5.0
FRAMEWORK_BINDER_WEIGHT = 2.0
# An observed positive T-cell assay is evidence, not prediction, so it is
# weighted separately and reported separately. It is still not a verdict: the
# record may come from a vaccine study or an unrelated donor context.
OBSERVED_EPITOPE_WEIGHT = 3.0

HTTP_TIMEOUT = 120


@dataclass
class PredictedBinder:
    """One predicted MHC-II binder that survived the germline filter."""

    allele: str
    peptide: str
    core: str
    rank: float
    start: int  # window start, 0-based index into full_sequence()
    end: int  # window end, exclusive
    core_start: int  # binding-core start, the span clustering works on
    core_end: int
    regions: str
    overlaps_cdr: bool
    core_status: list[str]


@dataclass
class ImmunogenicityResult:
    available: bool
    flags: list[Flag]
    note: str = ""
    n_predicted_total: int = 0
    n_dropped_germline: int = 0
    binders: list[PredictedBinder] = field(default_factory=list)


# --- IEDB prediction ------------------------------------------------------


def _cache_path(kind: str, key: str, cache_dir: Path | None) -> Path:
    directory = Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE_DIR
    directory.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(key.encode()).hexdigest()[:20]
    return directory / f"{kind}_{digest}.txt"


def predict_mhcii(
    sequence: str,
    alleles: list[str],
    cache_dir: Path | None = None,
) -> list[dict] | None:
    """Score every 15-mer window of `sequence` against `alleles` in ONE call.

    The API does its own windowing, so a whole chain and a comma-separated
    allele list go in one POST. Measured: a 120-residue VH against 5 alleles
    returns 530 rows in ~3s. Doing this per-window-per-allele instead would be
    530 separate requests for the same data.

    Returns a list of row dicts, or None if the service is unreachable.
    NOTE: rows come back sorted by score descending, NOT by position.
    """
    body = urllib.parse.urlencode(
        {
            "method": METHOD,
            "sequence_text": sequence,
            "allele": ",".join(alleles),
            "length": ",".join([str(WINDOW_SIZE)] * len(alleles)),
        }
    ).encode()

    cache_file = _cache_path("mhcii", f"{METHOD}|{sequence}|{','.join(alleles)}", cache_dir)
    if cache_file.exists() and cache_file.stat().st_size > 0:
        text = cache_file.read_text()
    else:
        try:
            request = urllib.request.Request(IEDB_MHCII_URL, data=body)
            with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT) as response:
                text = response.read().decode()
        except (urllib.error.URLError, TimeoutError, OSError):
            return None
        cache_file.write_text(text)

    lines = [line for line in text.strip().splitlines() if line.strip()]
    if len(lines) < 2:
        return None
    header = lines[0].split("\t")
    if "rank" not in header:
        return None

    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(header):
            continue
        row = dict(zip(header, fields))
        try:
            row["rank"] = float(row["rank"])
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
        except (ValueError, KeyError):
            continue
        rows.append(row)
    return rows or None


# --- IEDB experimental evidence (IQ-API) ----------------------------------


def observed_tcell_records(
    peptide: str, cache_dir: Path | None = None, limit: int = 50
) -> list[dict] | None:
    """Look up published T-cell assay records containing `peptide`.

    Substring match — exact 15-mer matching finds essentially nothing, since
    IEDB entries are whatever length the original study used.

    This is the *evidence* layer: it says a response was observed, not that
    it predicts one. Call it only on peptides that already survived the
    germline filter — each query takes ~0.8s, so scanning every window would
    cost minutes per candidate.
    """
    query = urllib.parse.urlencode(
        {
            "linear_sequence": f"like.*{peptide}*",
            "select": "structure_id,linear_sequence,qualitative_measure,source_organism_name",
            "limit": str(limit),
        }
    )
    cache_file = _cache_path("iq", query, cache_dir)
    if cache_file.exists() and cache_file.stat().st_size > 0:
        text = cache_file.read_text()
    else:
        try:
            with urllib.request.urlopen(
                f"{IEDB_IQ_TCELL_URL}?{query}", timeout=HTTP_TIMEOUT
            ) as response:
                text = response.read().decode()
        except (urllib.error.URLError, TimeoutError, OSError):
            return None
        cache_file.write_text(text)

    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return None


# --- the check itself -----------------------------------------------------


def check_immunogenicity(
    chain: NumberedChain,
    chain_name: str,
    alleles: list[str] | None = None,
    check_observed: bool = True,
    cache_dir: Path | None = None,
) -> ImmunogenicityResult:
    """Scan one chain for non-germline predicted MHC-II binders.

    Args:
        check_observed: also query IEDB for published T-cell assay records on
            the surviving peptides. ~0.8s per distinct core.
    """
    alleles = alleles or DEFAULT_ALLELES
    sequence = chain.full_sequence()
    residues = [r for r in chain.residues if r.aa != "-"]
    profile = germline_profile(chain)

    rows = predict_mhcii(sequence, alleles, cache_dir=cache_dir)
    if rows is None:
        return ImmunogenicityResult(
            available=False,
            flags=[],
            note=(
                f"IEDB MHC-II prediction API unreachable ({IEDB_MHCII_URL}) — "
                "immunogenicity not scored for this chain."
            ),
        )

    strong = [r for r in rows if r["rank"] < BINDER_PERCENTILE_CUTOFF]

    binders: list[PredictedBinder] = []
    dropped_germline = 0
    for row in strong:
        window_start = row["start"] - 1  # API positions are 1-based
        window_end = row["end"]
        core = row.get("core_peptide", "")

        # Locate the 9-mer binding core inside the window. MHC-II binding is
        # determined by the core sitting in the groove; the flanks hang
        # outside it and do not decide whether the peptide is presented. This
        # distinction is the difference between a filter that works and one
        # that keeps everything.
        core_start = sequence.find(core, window_start) if core else -1
        if core_start == -1:
            core_start, core_end = window_start, window_end
        else:
            core_end = core_start + len(core)

        if profile.available and profile.is_self(core_start, core_end):
            dropped_germline += 1
            continue

        window_residues = residues[window_start:window_end]
        overlaps_cdr = any(r.is_cdr for r in window_residues)
        binders.append(
            PredictedBinder(
                allele=row["allele"],
                peptide=row["peptide"],
                core=core,
                rank=row["rank"],
                start=window_start,
                end=window_end,
                core_start=core_start,
                core_end=core_end,
                regions="/".join(sorted({r.region for r in window_residues})),
                overlaps_cdr=overlaps_cdr,
                core_status=profile.window_status(core_start, core_end),
            )
        )

    flags = _flags_for_binders(binders, chain_name, profile.available)

    if check_observed and binders:
        flags.extend(_observed_flags(binders, chain_name, cache_dir))

    note = ""
    if not profile.available:
        note = (
            "no germline assignment for this chain — every predicted binder is "
            "reported, including ones central tolerance would cover. Treat the "
            "count as an upper bound."
        )

    return ImmunogenicityResult(
        available=True,
        flags=flags,
        note=note,
        n_predicted_total=len(strong),
        n_dropped_germline=dropped_germline,
        binders=binders,
    )


def _cluster_binders(binders: list[PredictedBinder]) -> list[list[PredictedBinder]]:
    """Group binders whose CORES overlap in sequence into one cluster each.

    Deduplicating by exact core string is not enough. A single liable stretch
    produces several shifted cores as the window slides across it — measured
    on a real candidate: 'LLISAASSL', 'LISAASSLQ', 'ISAASSLQS' are three
    offsets of one region, and scoring them as three independent findings
    weights that one liability three times over. Cluster by overlap, then
    emit one flag per cluster.
    """
    if not binders:
        return []

    # Core span within the chain, recovered from the window plus the core's
    # offset inside it; binders that share any residue belong together.
    ordered = sorted(binders, key=lambda b: (b.core_start, b.core_end))
    clusters: list[list[PredictedBinder]] = [[ordered[0]]]
    reach = ordered[0].core_end
    for binder in ordered[1:]:
        if binder.core_start < reach:
            clusters[-1].append(binder)
            reach = max(reach, binder.core_end)
        else:
            clusters.append([binder])
            reach = binder.core_end
    return clusters


def _flags_for_binders(
    binders: list[PredictedBinder], chain_name: str, germline_available: bool
) -> list[Flag]:
    """One flag per overlapping cluster of cores — never one per allele or window."""
    flags = []
    for cluster in _cluster_binders(binders):
        best = min(cluster, key=lambda b: b.rank)
        n_alleles = len({b.allele for b in cluster})
        n_cores = len({b.core for b in cluster})
        caveat = "" if germline_available else " (germline unknown)"
        extra = f", {n_cores} overlapping cores" if n_cores > 1 else ""
        flags.append(
            Flag(
                check="immunogenicity",
                severity="soft",
                region=best.regions,
                weight=CDR_BINDER_WEIGHT if best.overlaps_cdr else FRAMEWORK_BINDER_WEIGHT,
                message=(
                    f"predicted non-germline MHC-II binder core '{best.core}' "
                    f"({chain_name} {best.regions}) best rank={best.rank} "
                    f"across {n_alleles} allele(s){extra}"
                    f"{' — overlaps CDR' if best.overlaps_cdr else ''}{caveat}"
                ),
            )
        )
    return flags


def _observed_flags(
    binders: list[PredictedBinder], chain_name: str, cache_dir: Path | None
) -> list[Flag]:
    """Evidence layer: has this core ever been seen to provoke a T-cell response?

    Reported separately from the prediction flags and never summed with them —
    "might bind MHC-II" and "was observed to provoke a response" are different
    claims with different standards of proof.
    """
    flags = []
    # One evidence lookup per cluster, using its strongest core — querying
    # every shifted core of the same region would repeat one finding and cost
    # an extra ~0.8s round trip each time.
    representatives = [
        min(cluster, key=lambda b: b.rank).core for cluster in _cluster_binders(binders)
    ]
    for core in sorted({c for c in representatives if c}):
        records = observed_tcell_records(core, cache_dir=cache_dir)
        if not records:
            continue
        positive = [
            r
            for r in records
            if str(r.get("qualitative_measure") or "").startswith("Positive")
        ]
        if not positive:
            continue
        # IQ-API returns JSON nulls, not missing keys, so .get(k, default)
        # still hands back None here — coerce before sorting.
        organisms = sorted(
            {str(r.get("source_organism_name") or "unknown source") for r in positive}
        )[:2]
        flags.append(
            Flag(
                check="immunogenicity_observed",
                severity="soft",
                region="evidence",
                weight=OBSERVED_EPITOPE_WEIGHT,
                message=(
                    f"core '{core}' ({chain_name}) matches {len(positive)} positive "
                    f"T-cell assay record(s) in IEDB, e.g. {', '.join(organisms)} — "
                    "observed, not predicted; context may be unrelated to therapeutic use"
                ),
            )
        )
    return flags
