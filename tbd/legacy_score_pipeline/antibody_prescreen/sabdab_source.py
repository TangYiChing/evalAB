"""Load real VH/VL(/antigen) sequences from the ANTIPASTI-curated SAbDab
dataset (kevinmicha/ANTIPASTI on GitHub), as a replacement data source for
the calibration set — TheraSAbDab, SAbDab, RCSB, and PDBe were all
unreachable from this build environment's network egress allowlist, but
this dataset's residue-list files are plain GitHub-hosted .npy arrays and
work directly.

Coverage: 627 unique PDB entries in sabdab_summary_all.tsv have complete
Hchain/Lchain/antigen_chain; 617 (98.4%) have a matching lists_of_residues
npy file. A PDB with multiple antibody copies in the asymmetric unit (e.g.
5w08, 5 Fab copies) has ONE npy file covering exactly one copy — this
module matches the npy's actual chain letters against the TSV rows to find
the correct one, rather than assuming the first TSV row for that PDB ID.

Token format per residue, e.g. 'DH  1 ': [amino acid][chain letter][position].
'START-Ab' / 'END-Ab' mark the antibody (H+L) block; antigen residues follow
END-Ab. Position numbers can have gaps (unresolved residues in the crystal
structure) — sequences are built by residue order in the array, not by
filling numeric gaps.
"""

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class SabdabCandidate:
    pdb_id: str
    heavy_chain_id: str
    light_chain_id: str
    antigen_chain_id: str
    antigen_type: str
    antigen_name: str
    vh_sequence: str
    vl_sequence: str
    antigen_sequence: str


def load_summary_rows(tsv_path: str) -> list[dict]:
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def complete_rows(rows: list[dict]) -> list[dict]:
    """Rows with Hchain, Lchain, and antigen_chain all present (not 'NA'/empty)."""
    return [
        r
        for r in rows
        if r.get("Hchain") not in ("NA", "", None)
        and r.get("Lchain") not in ("NA", "", None)
        and r.get("antigen_chain") not in ("NA", "", None)
    ]


def _parse_residue_array(arr: np.ndarray) -> dict[str, str]:
    """Parse one lists_of_residues/{pdb}.npy array into per-chain sequences.

    Returns {chain_letter: sequence}. Antibody chains (before END-Ab) and
    antigen chains (after END-Ab) are both included — caller matches chain
    letters against the TSV row to pick out VH/VL/antigen.
    """
    chains: dict[str, list[str]] = {}
    for token in arr:
        if token in ("START-Ab", "END-Ab"):
            continue
        aa, chain_letter = token[0], token[1]
        if aa not in STANDARD_AA:
            continue  # skip HETATM-derived or non-standard residue codes
        chains.setdefault(chain_letter, []).append(aa)
    return {chain: "".join(residues) for chain, residues in chains.items()}


def load_candidate(pdb_id: str, row: dict, residues_dir: Path) -> SabdabCandidate | None:
    """Load one candidate given its matched TSV row. Returns None if the npy
    file is missing or doesn't contain the row's expected chain letters
    (e.g. this npy covers a different antibody copy of the same PDB ID).
    """
    npy_path = residues_dir / f"{pdb_id}.npy"
    if not npy_path.exists():
        return None

    arr = np.load(npy_path, allow_pickle=True)
    chains = _parse_residue_array(arr)

    h, l, ag = row["Hchain"], row["Lchain"], row["antigen_chain"]
    if h not in chains or l not in chains or ag not in chains:
        return None

    return SabdabCandidate(
        pdb_id=pdb_id,
        heavy_chain_id=h,
        light_chain_id=l,
        antigen_chain_id=ag,
        antigen_type=row.get("antigen_type", ""),
        antigen_name=row.get("antigen_name", ""),
        vh_sequence=chains[h],
        vl_sequence=chains[l],
        antigen_sequence=chains[ag],
    )


def load_all_candidates(tsv_path: str, residues_dir: str) -> list[SabdabCandidate]:
    """Load every complete, npy-matched candidate. Rows whose npy is missing
    or doesn't match (wrong antibody copy) are silently skipped — call
    complete_rows() + iterate manually if you need to see what was dropped.
    """
    rows = complete_rows(load_summary_rows(tsv_path))
    residues_path = Path(residues_dir)
    candidates = []
    for row in rows:
        candidate = load_candidate(row["pdb"], row, residues_path)
        if candidate is not None:
            candidates.append(candidate)
    return candidates
