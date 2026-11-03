"""Independently verify files written by R (development-only oracle).

Reads R-written H5AD / Zarr / H5MU / MuData-Zarr files with the Python
reference implementations (anndata, mudata, awkward, scirpy) and compares them
with the expected JSON produced by the fixture generators.

Usage (from tools/fixtures):
    .venv/bin/python -I verify_r_outputs.py anndata <path> <expected.json> [--rows 1,3,5] [--edge]
    .venv/bin/python -I verify_r_outputs.py mudata  <path> <expected_airr.json> [--rows ...]

`--rows` gives the zero-based observation subset that R applied. Exits with a
non-zero status on the first mismatch.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import anndata as ad
import awkward as ak
import mudata as md
import numpy as np

warnings.filterwarnings("ignore")


def fail(msg: str) -> None:
    print(f"FAIL: {msg}")
    sys.exit(1)


def read_any_anndata(path: str) -> ad.AnnData:
    return ad.read_zarr(path) if path.rstrip("/").endswith(".zarr") else ad.read_h5ad(path)


def jsonable(x):
    """ak.to_list returns tuples for tuple records; JSON fixtures hold lists."""
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, dict):
        return {k: jsonable(v) for k, v in x.items()}
    return x


def check_awkward(adata: ad.AnnData, key: str, expected_values, expected_type: str | None, rows) -> None:
    if key not in adata.obsm:
        fail(f"obsm[{key!r}] missing")
    arr = adata.obsm[key]
    if not isinstance(arr, ak.Array):
        fail(f"obsm[{key!r}] is {type(arr)}, not an awkward Array")
    exp = expected_values if rows is None else [expected_values[i] for i in rows]
    got = jsonable(ak.to_list(arr))
    if got != exp:
        for i, (g, e) in enumerate(zip(got, exp)):
            if g != e:
                fail(f"obsm[{key!r}] row {i}: {g!r} != {e!r}")
        fail(f"obsm[{key!r}]: length {len(got)} != {len(exp)}")
    if expected_type is not None:
        n_exp = len(exp)
        want = f"{n_exp} * " + expected_type.split(" * ", 1)[1]
        if str(arr.type) != want:
            fail(f"obsm[{key!r}] type {arr.type} != {want}")


def check_scirpy(adata: ad.AnnData) -> None:
    import scirpy as ir

    # exercises scirpy's own accessors on the R-written awkward data
    ir.tl.chain_qc(adata)
    loci = ir.get.airr(adata, "locus")
    if loci.shape[0] != adata.n_obs:
        fail("scirpy.get.airr returned the wrong number of cells")


def verify_anndata(args) -> None:
    expected = json.loads(Path(args.expected).read_text())
    adata = read_any_anndata(args.path)
    rows = [int(r) for r in args.rows.split(",")] if args.rows else None
    obs_exp = expected["obs_names"] if rows is None else [expected["obs_names"][i] for i in rows]
    if list(adata.obs_names) != obs_exp:
        fail(f"obs_names {list(adata.obs_names)} != {obs_exp}")
    if args.edge:
        for key, values in expected["values"].items():
            check_awkward(adata, key, values, expected["types"][key], rows)
    else:
        check_awkward(adata, "airr", expected["airr"], expected["airr_type"], rows)
        check_awkward(adata, "chain_indices", expected["chain_indices"], expected["chain_indices_type"], rows)
        check_scirpy(adata)
    print(f"OK: {args.path}")


def verify_mudata(args) -> None:
    expected = json.loads(Path(args.expected).read_text())
    path = args.path
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        mdata = md.read_zarr(path) if path.rstrip("/").endswith(".zarr") else md.read_h5mu(path)
    for w in caught:
        if "not created by muon/mudata" in str(w.message) or "mudata" in str(w.message).lower() and "created" in str(w.message):
            fail(f"mudata warned: {w.message}")
    em = expected["mudata"]
    if list(mdata.mod.keys()) != em["mod_order"]:
        fail(f"modality order {list(mdata.mod.keys())} != {em['mod_order']}")
    if list(mdata.obs_names) != em["obs_names"]:
        fail(f"obs_names {list(mdata.obs_names)} != {em['obs_names']}")
    if list(mdata.var_names) != em["var_names"]:
        fail(f"var_names {list(mdata.var_names)} != {em['var_names']}")
    for k, v in em["obsmap"].items():
        if np.asarray(mdata.obsmap[k]).tolist() != v:
            fail(f"obsmap[{k}] differs")
    for k, v in em["varmap"].items():
        if np.asarray(mdata.varmap[k]).tolist() != v:
            fail(f"varmap[{k}] differs")
    airr = mdata.mod["airr"]
    if list(airr.obs_names) != expected["obs_names"]:
        fail("airr modality obs_names differ")
    check_awkward(airr, "airr", expected["airr"], expected["airr_type"], None)
    check_awkward(airr, "chain_indices", expected["chain_indices"], expected["chain_indices_type"], None)
    check_scirpy(airr)
    # scirpy also accepts the MuData directly
    import scirpy as ir

    ir.tl.chain_qc(mdata)
    print(f"OK: {path}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("kind", choices=["anndata", "mudata"])
    p.add_argument("path")
    p.add_argument("expected")
    p.add_argument("--rows", default=None)
    p.add_argument("--edge", action="store_true")
    args = p.parse_args()
    (verify_anndata if args.kind == "anndata" else verify_mudata)(args)


if __name__ == "__main__":
    main()
