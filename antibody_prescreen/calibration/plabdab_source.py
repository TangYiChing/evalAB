"""PLAbDab as the reference population the range bands are fitted to.

A reference population decides what "unusual" means, so its selection bias is
inherited wholesale by every band fitted on it. Choosing one is therefore a
design decision that has to be stated, not a data-loading detail.

PLAbDab is chosen because it is large (~150k paired entries) and because it
carries a *provenance label*. The obvious alternative — antibodies with solved
structures — is a population of a few hundred selected for "someone
crystallised it", which is a structural-biology selection dressed up as a
developability one, and is far too small to assert a distribution tail.

## Do not download the 5 GB archive

PLAbDab's documented download is `plabdab_data.tar.gz`, ~5.05 GB. Almost all
of that is `data-dev/models/*.pdb` — ~65k pre-built ABodyBuilder2 models,
~280 KB each — and those members come FIRST in the tar, so streaming and
aborting early does not reach the sequence CSVs either.

But the individual files are served standalone from the same directory, and
`paired_sequences.csv.gz` is 11 MB. That is the only file needed to calibrate,
so that is the only file this module fetches. The `PLAbDab` python package and
its KA-Search index are NOT dependencies here — they serve sequence/structure
*search*, which is a separate feature with a separate (large) download.

## The label, and what it is not

The `pairing` column records how each entry was assembled. Two of its values
are useful as weak, opposing labels:

  - "TheraSAbDab"  (~1.2k)  clinical-stage therapeutics
  - "Patent text"  (~91k)   merely patented, no developability filter applied

"Reached the clinic" is not "manufacturable", and "appears in a patent" is
not "bad". This is a provenance label, not an assay. It is a better signal
than "a structure exists", and it is still not experimental ground truth.
"""

import os
import urllib.request
from pathlib import Path

PLABDAB_PAIRED_URL = (
    "https://opig.stats.ox.ac.uk/webapps/plabdab/static/downloads/"
    "paired_sequences.csv.gz"
)

# Where the reference download lives. In a clone, that is `data/reference/`,
# so the file sits beside the repository it calibrates and a reader can see
# exactly which public data produced the shipped bands. Outside a clone
# (pip-installed), it falls back to the user cache.
_REPO_REFERENCE_DIR = Path(__file__).resolve().parents[2] / "data" / "reference"
_USER_CACHE_DIR = Path.home() / ".cache" / "evalab" / "plabdab"


def default_cache_path() -> Path:
    """Resolve the reference file location: repo `data/reference/` if present."""
    override = os.environ.get("EVALAB_REFERENCE_DIR")
    if override:
        return Path(override).expanduser() / "paired_sequences.csv.gz"
    if _REPO_REFERENCE_DIR.is_dir():
        return _REPO_REFERENCE_DIR / "paired_sequences.csv.gz"
    return _USER_CACHE_DIR / "paired_sequences.csv.gz"

# `pairing` values used as the two calibration populations.
THERAPEUTIC_PAIRING = "TheraSAbDab"
BACKGROUND_PAIRING = "Patent text"

# Shortest plausible V-domain. Anything below this is a fragment, not a chain,
# and will only waste an ANARCI call.
MIN_CHAIN_LEN = 90


def fetch_paired_sequences(cache_path: Path | str | None = None) -> Path:
    """Download PLAbDab's paired_sequences.csv.gz (~11 MB) if not already cached."""
    path = Path(cache_path) if cache_path is not None else default_cache_path()
    if path.exists() and path.stat().st_size > 0:
        return path

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".partial")
    urllib.request.urlretrieve(PLABDAB_PAIRED_URL, tmp)
    tmp.replace(path)
    return path


def load_paired_sequences(cache_path: Path | str | None = None):
    """Return the full PLAbDab paired table as a pandas DataFrame."""
    import pandas as pd

    return pd.read_csv(fetch_paired_sequences(cache_path))


def sample_pilot(
    n_per_group: int = 50,
    seed: int = 0,
    cache_path: Path | str | None = None,
) -> list[dict]:
    """Build a small balanced pilot set for proving the pipeline runs.

    Returns candidates in the shape `screen_batch` expects, plus a `pairing`
    key recording which population each came from.

    Deduplicates on (heavy, light): PLAbDab's "Same entry" / "Unique source"
    pairings mean one antibody can appear under several IDs, and scoring the
    same sequence twice would silently weight it twice in any distribution.
    """
    df = load_paired_sequences(cache_path)

    df = df[df["pairing"].isin([THERAPEUTIC_PAIRING, BACKGROUND_PAIRING])]
    df = df.dropna(subset=["heavy_sequence", "light_sequence"])
    df = df[
        (df["heavy_sequence"].str.len() >= MIN_CHAIN_LEN)
        & (df["light_sequence"].str.len() >= MIN_CHAIN_LEN)
    ]
    df = df.drop_duplicates(subset=["heavy_sequence", "light_sequence"])

    out = []
    for pairing in (THERAPEUTIC_PAIRING, BACKGROUND_PAIRING):
        group = df[df["pairing"] == pairing]
        take = min(n_per_group, len(group))
        for row in group.sample(n=take, random_state=seed).itertuples():
            out.append(
                {
                    "candidate_id": str(row.ID),
                    "vh_sequence": row.heavy_sequence,
                    "vl_sequence": row.light_sequence,
                    "pairing": pairing,
                }
            )
    return out
