"""ABodyBuilder2 structure prediction, with an on-disk model cache.

This is the expensive step of Tier 2 — everything downstream (TAP) is cheap
by comparison. Calibration re-runs the scoring many times over the same
population, so models are cached on disk and keyed by sequence content, not
by candidate id: two candidates with identical VH/VL share one model, and
renaming a candidate does not invalidate its model.
"""

import hashlib
from pathlib import Path

DEFAULT_CACHE_DIR = Path.home() / ".cache" / "evalab" / "model_cache"


class ModellingError(Exception):
    """Raised when ABodyBuilder2 cannot produce a model for a VH/VL pair.

    Callers should degrade to Tier-1-only for that candidate rather than
    failing the whole batch — same contract as immunogenicity's
    `available=False`.
    """


def model_key(vh_sequence: str, vl_sequence: str) -> str:
    """Content-addressed cache key for a VH/VL pair."""
    digest = hashlib.sha256(
        f"{vh_sequence.strip().upper()}|{vl_sequence.strip().upper()}".encode()
    )
    return digest.hexdigest()[:16]


def build_model(
    vh_sequence: str,
    vl_sequence: str,
    cache_dir: Path | str | None = None,
    predictor=None,
) -> Path:
    """Predict an IMGT-numbered Fv structure and return the path to its .pdb.

    Args:
        vh_sequence: heavy chain variable domain sequence.
        vl_sequence: light chain variable domain sequence.
        cache_dir: where models are cached. Defaults to DEFAULT_CACHE_DIR.
        predictor: an already-constructed ABodyBuilder2 instance. Constructing
            one loads ~100MB of weights, so batch callers should build it once
            and pass it in (see `get_predictor`).

    Returns:
        Path to the model .pdb. ABodyBuilder2 writes IMGT-numbered output with
        chains named H and L, which is exactly what TAP expects.

    Raises:
        ModellingError: if ImmuneBuilder is not installed, or prediction fails.
    """
    cache_dir = Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)

    path = cache_dir / f"{model_key(vh_sequence, vl_sequence)}.pdb"
    if path.exists() and path.stat().st_size > 0:
        return path

    if predictor is None:
        predictor = get_predictor()

    try:
        antibody = predictor.predict({"H": vh_sequence, "L": vl_sequence})
    except Exception as e:  # ImmuneBuilder raises a variety of bare exceptions
        raise ModellingError(f"ABodyBuilder2 prediction failed: {e}") from e

    # Write via a temp name so an interrupted run cannot leave a truncated
    # model in the cache that later runs would happily reuse.
    tmp = path.with_suffix(".pdb.partial")
    try:
        antibody.save(str(tmp))
        tmp.replace(path)
    except Exception as e:
        tmp.unlink(missing_ok=True)
        raise ModellingError(f"Could not write model for VH/VL pair: {e}") from e

    return path


def get_predictor():
    """Construct an ABodyBuilder2 predictor, with a clear error if unavailable.

    Loads model weights (downloaded to ~/.cache on first use), so construct
    once per batch rather than per candidate.
    """
    try:
        from ImmuneBuilder import ABodyBuilder2
    except ImportError as e:
        raise ModellingError(
            "ImmuneBuilder is not installed — Tier 2 structure screening is "
            "unavailable. Install it with `pip install ImmuneBuilder`, or run "
            "with run_structure=False."
        ) from e

    try:
        return ABodyBuilder2()
    except Exception as e:
        raise ModellingError(f"Could not initialise ABodyBuilder2: {e}") from e
