"""Awkward edge-case fixtures not produced by scirpy (development-only).

Covers node classes and null/empty distinctions the scirpy fixtures do not
exercise: ByteMaskedArray, BitMaskedArray, ListArray (non-contiguous),
IndexedArray, UnmaskedArray, RegularArray, tuples, narrow integer and float32
primitives, empty strings vs missing strings vs empty lists vs missing rows.

Usage (from tools/fixtures):
    .venv/bin/python -I generate_awkward_edge_fixtures.py <outdir>
"""

from __future__ import annotations

import json
import shutil
import sys
import warnings
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import awkward as ak
import h5py
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

N = 5


def edge_arrays() -> dict[str, ak.Array]:
    i64 = ak.index.Index64
    out = {}
    # ?var * ?string: missing row, empty list, empty string, missing string
    out["strings"] = ak.Array([["a", ""], [], None, [None, "xyz"], [""]])
    # ByteMaskedArray (valid_when=True) over int64
    out["bytemasked"] = ak.Array(
        ak.contents.ByteMaskedArray(
            ak.index.Index8(np.array([1, 0, 1, 0, 1], dtype=np.int8)),
            ak.contents.NumpyArray(np.array([10, 20, 30, 40, 50], dtype=np.int64)),
            valid_when=True,
        )
    )
    # BitMaskedArray, valid_when=False, lsb_order=True over float32
    out["bitmasked"] = ak.Array(
        ak.contents.BitMaskedArray(
            ak.index.IndexU8(np.packbits(np.array([0, 1, 0, 0, 1, 0, 0, 0], dtype=bool), bitorder="little")),
            ak.contents.NumpyArray(np.array([1.5, 2.5, 3.5, 4.5, 5.5], dtype=np.float32)),
            valid_when=False,
            length=N,
            lsb_order=True,
        )
    )
    # ListArray with non-contiguous, overlapping ranges over int32
    out["listarray"] = ak.Array(
        ak.contents.ListArray(
            i64(np.array([3, 0, 2, 2, 5])),
            i64(np.array([5, 2, 2, 4, 6])),
            ak.contents.NumpyArray(np.array([0, 1, 2, 3, 4, 5], dtype=np.int32)),
        )
    )
    # RegularArray of uint8
    out["regular"] = ak.to_regular(ak.Array(np.arange(10, dtype=np.uint8).reshape(N, 2)), axis=1)
    # IndexedArray (non-option) over strings
    out["indexed"] = ak.Array(
        ak.contents.IndexedArray(
            i64(np.array([2, 0, 1, 2, 0])),
            ak.Array(["x", "yy", ""]).layout,
        )
    )
    # UnmaskedArray over int16
    out["unmasked"] = ak.Array(ak.contents.UnmaskedArray(ak.contents.NumpyArray(np.array([1, -2, 3, -4, 5], dtype=np.int16))))
    # tuple records
    out["tuples"] = ak.Array([(1, "a"), (2, ""), (3, "c"), (4, "d"), (5, "e")])
    # nested records with lists of uint16 and bools
    out["nested"] = ak.Array(
        [
            {"x": [1, 2], "flag": True, "name": "p"},
            {"x": [], "flag": False, "name": ""},
            {"x": [3], "flag": True, "name": None},
            {"x": [4, 5, 6], "flag": False, "name": "q"},
            {"x": [7], "flag": True, "name": "r"},
        ],
        with_name=None,
    )
    return out


def main(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for name in ("edge_cases.h5ad", "union.h5ad"):
        (outdir / name).unlink(missing_ok=True)
    for zf in (2, 3):
        shutil.rmtree(outdir / f"edge_cases_v{zf}.zarr", ignore_errors=True)

    obs = pd.DataFrame(index=[f"cell{i}" for i in range(N)])
    arrays = edge_arrays()
    adata = ad.AnnData(obs=obs, obsm=arrays)
    adata.write_h5ad(outdir / "edge_cases.h5ad")
    for zf in (2, 3):
        ad.settings.zarr_write_format = zf
        adata.write_zarr(outdir / f"edge_cases_v{zf}.zarr")
    ad.settings.zarr_write_format = 3

    union = ad.AnnData(obs=obs, obsm={"union": ak.Array([1, "a", 2, "b", 3])})
    union.write_h5ad(outdir / "union.h5ad")

    forms = {}
    with h5py.File(outdir / "edge_cases.h5ad", "r") as f:
        for k in f["obsm"]:
            forms[k] = json.loads(f["obsm"][k].attrs["form"])
    with h5py.File(outdir / "union.h5ad", "r") as f:
        forms["union"] = json.loads(f["obsm/union"].attrs["form"])

    expected = {
        "obs_names": list(adata.obs_names),
        "values": {k: ak.to_list(v) for k, v in arrays.items()},
        "types": {k: str(v.type) for k, v in arrays.items()},
        "forms": forms,
        "versions": {p: version(p) for p in ["anndata", "awkward", "zarr", "h5py", "numpy"]},
    }
    (outdir / "expected_edge_cases.json").write_text(json.dumps(expected, indent=1))
    print(json.dumps(expected["types"], indent=1))


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "generated" / "awkward")
