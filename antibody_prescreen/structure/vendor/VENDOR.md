# Vendored third-party code

## `tap/` — TAP (Therapeutic Antibody Profiler) implementation

- **Upstream**: https://github.com/Exscientia/ab-characterisation
- **Path**: `src/ab_characterisation/developability_tools/tap/`
- **Pinned commit**: `46ccd6452ec22e31c6c4a97327740b7a88ab0b83`
- **Licence**: MIT (Copyright (c) 2024 Exscientia) — full text in `tap/LICENSE.upstream`

TAP itself is from Raybould et al. (2019), "Five computational developability
guidelines for therapeutic antibody profiling", PNAS 116(10):4025-4030.

### Why vendored rather than pip-installed

`ab-characterisation` has no PyPI release, and installing the whole package
from git pulls in its Rosetta / ChimeraX / mpi4py pipeline stages, which this
project does not use. The TAP subtree only needs `biopython` + `loguru` plus
the bundled `psa` binary, so vendoring one MIT-licensed subtree at a recorded
commit is the smaller commitment.

### Local edits

Only two, both mechanical:

1. **Import rewrites.** Every `from ab_characterisation.developability_tools.tap.X`
   absolute import became the equivalent relative import. No logic changed.
2. **`outputs.py`**: upstream's `write_output_file` delegates to
   `ab_characterisation.developability_tools.utils.outputs.write_file`, which
   also handles `s3://` destinations. Replaced with a plain local write, marked
   in-file with a `VENDOR EDIT` comment, rather than vendoring the S3 helper
   and its dependencies.

Metric definitions, thresholds, and the `psa` binaries are byte-for-byte
upstream.

### `psa_executables/`

`psa` version 2.0 (Feb 2000) computes per-residue solvent-accessible surface
area. Two builds are shipped and `structure_annotation.py` picks between them
on `sys.platform`:

| File | Format | Used on |
|---|---|---|
| `psa` | ELF 64-bit x86-64, Linux | everything except macOS |
| `psa_mac` | Mach-O 64-bit x86_64 | macOS (incl. Apple Silicon via Rosetta 2) |

There is no arm64 build. On Apple Silicon, `psa_mac` runs under Rosetta 2 —
install it once with:

```bash
softwareupdate --install-rosetta --agree-to-license
```

Both files must keep their executable bit. If it is lost (e.g. by a checkout
that drops the mode), TAP fails with `PSAError("psa executable was not found.")`.
Restore with `chmod +x` and `git update-index --chmod=+x`.

### Updating

Re-download the same file list from a newer upstream commit, re-apply the two
edits above, and update the pinned SHA in this file. Do not edit the vendored
files for project-specific behaviour — put that in
`antibody_prescreen/structure/tap_runner.py` instead.
