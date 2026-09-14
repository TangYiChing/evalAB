"""Own-sequence MHC-II binding scan for immunogenicity risk.

This is the corrected approach discussed in the design phase: score the
antibody's OWN sequence (15-mer sliding window) against common human HLA-DR
alleles, rather than the ToolUniverse skill's `iedb_search_epitopes(target=...)`
call, which searches known epitopes ON the target antigen — a different,
mostly-irrelevant question for drug immunogenicity.

Network dependency: calls IEDB's MHC-II binding prediction API
(tools-cluster-interface.iedb.org). In THIS sandbox that host is not on the
network egress allowlist (confirmed via the agent proxy status endpoint), so
this check degrades gracefully — it returns an "unavailable" result rather
than crashing the pipeline. In a real deployment, allowlist that host (or
swap in a local NetMHCIIpan install) to get real signal here.
"""

from dataclasses import dataclass

import urllib.error
import urllib.request

from .checks import Flag
from .numbering import NumberedChain

IEDB_MHCII_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"

# A small default panel of common human HLA-DR alleles, chosen for broad
# population coverage. Extend this list for a more thorough scan; each
# allele added multiplies the number of API calls.
DEFAULT_ALLELES = [
    "HLA-DRB1*01:01",
    "HLA-DRB1*03:01",
    "HLA-DRB1*04:01",
    "HLA-DRB1*07:01",
    "HLA-DRB1*15:01",
]

WINDOW_SIZE = 15
# netmhciipan reports percentile rank; lower = stronger predicted binder.
# <2.0 is the field's conventional "likely binder" cutoff.
BINDER_PERCENTILE_CUTOFF = 2.0


@dataclass
class ImmunogenicityResult:
    available: bool
    flags: list[Flag]
    note: str = ""


def _windows(seq: str, size: int = WINDOW_SIZE):
    for i in range(len(seq) - size + 1):
        yield i, seq[i : i + size]


def _query_iedb(peptide: str, allele: str) -> float | None:
    """Return predicted binding percentile rank, or None if unavailable."""
    data = f"method=netmhciipan_el&sequence_text={peptide}&allele={allele}&length={len(peptide)}".encode()
    req = urllib.request.Request(IEDB_MHCII_URL, data=data)
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            body = resp.read().decode()
    except (urllib.error.URLError, TimeoutError, OSError):
        return None

    # Response is TSV; percentile rank is typically the last numeric column.
    # Real parsing depends on the live response shape — verify against a real
    # response once this host is allowlisted, this is a best-effort parser.
    lines = [ln for ln in body.strip().splitlines() if ln and not ln.startswith("allele")]
    if not lines:
        return None
    try:
        fields = lines[0].split("\t")
        return float(fields[-1])
    except (ValueError, IndexError):
        return None


def check_immunogenicity(
    chain: NumberedChain, chain_name: str, alleles: list[str] | None = None
) -> ImmunogenicityResult:
    alleles = alleles or DEFAULT_ALLELES
    seq = chain.full_sequence()
    residues = [r for r in chain.residues if r.aa != "-"]
    flags: list[Flag] = []
    any_call_succeeded = False

    for start_idx, peptide in _windows(seq):
        for allele in alleles:
            rank = _query_iedb(peptide, allele)
            if rank is None:
                continue
            any_call_succeeded = True
            if rank < BINDER_PERCENTILE_CUTOFF:
                window_residues = residues[start_idx : start_idx + WINDOW_SIZE]
                overlaps_cdr = any(r.is_cdr for r in window_residues)
                region_label = "/".join(sorted({r.region for r in window_residues}))
                flags.append(
                    Flag(
                        check="immunogenicity",
                        severity="soft",
                        region=region_label,
                        weight=5.0 if overlaps_cdr else 1.0,
                        message=(
                            f"predicted MHC-II binder '{peptide}' ({chain_name} {region_label}) "
                            f"allele={allele} percentile_rank={rank}"
                            f"{' — overlaps CDR' if overlaps_cdr else ''}"
                        ),
                    )
                )

    if not any_call_succeeded:
        return ImmunogenicityResult(
            available=False,
            flags=[],
            note=(
                "IEDB MHC-II API unreachable from this environment (host not on network "
                "egress allowlist) — immunogenicity check skipped, not scored. Allowlist "
                "tools-cluster-interface.iedb.org or install a local NetMHCIIpan to enable."
            ),
        )

    return ImmunogenicityResult(available=True, flags=flags)
