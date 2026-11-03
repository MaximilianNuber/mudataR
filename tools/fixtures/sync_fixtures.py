"""Copy generated fixtures into the R packages (development-only).

Generators write raw outputs to tools/fixtures/generated/. Zarr stores are
zipped when copied, because R CMD check flags the hidden metadata files of
Zarr v2 stores (anndataR ships its Zarr fixtures zipped for the same reason).

Usage (from tools/fixtures):
    .venv/bin/python -I sync_fixtures.py [path/to/anndataR]
"""

from __future__ import annotations

import shutil
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
GEN = HERE / "generated"
MUDATAR = HERE.parent.parent / "inst" / "extdata"

MUDATAR_FILES = [
    "scirpy/gex_airr.h5mu",
    "scirpy/gex_airr_v2.zarr",
    "scirpy/gex_airr_v3.zarr",
    "scirpy/expected_airr.json",
    "scirpy/manifest.json",
]
ANNDATAR_FILES = [
    "scirpy/airr.h5ad",
    "scirpy/airr_v2.zarr",
    "scirpy/airr_v3.zarr",
    "scirpy/gex_airr.h5mu",
    "scirpy/expected_airr.json",
    "scirpy/manifest.json",
    "awkward/edge_cases.h5ad",
    "awkward/edge_cases_v2.zarr",
    "awkward/edge_cases_v3.zarr",
    "awkward/union.h5ad",
    "awkward/expected_edge_cases.json",
]


def zip_store(src: Path, dest: Path) -> None:
    with zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(src.rglob("*")):
            zf.write(path, path.relative_to(src.parent))


def sync(files: list[str], dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    for rel in files:
        src = GEN / rel
        if src.is_dir():
            out = dest / (src.name + ".zip")
            shutil.rmtree(dest / src.name, ignore_errors=True)
            zip_store(src, out)
        else:
            out = dest / src.name
            shutil.copy2(src, out)
        print(f"{rel} -> {out}")


if __name__ == "__main__":
    sync(MUDATAR_FILES, MUDATAR)
    if len(sys.argv) > 1:
        sync(ANNDATAR_FILES, Path(sys.argv[1]) / "tests" / "testthat" / "fixtures" / "awkward")
