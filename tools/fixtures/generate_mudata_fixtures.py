"""Differential oracle cases for MuData index/map construction (development-only).

Each case builds a MuData with Python mudata from modality obs/var names, runs
optional steps (replace a modality, update) and records the resulting global
names, maps and global obs columns. R tests rebuild the same objects natively
and compare.

Usage (from tools/fixtures):
    .venv/bin/python -I generate_mudata_fixtures.py ../../tests/testthat/fixtures/expected_update.json
"""

from __future__ import annotations

import json
import sys
import warnings
from importlib.metadata import version

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


def adata(obs, var):
    return ad.AnnData(obs=pd.DataFrame(index=obs), var=pd.DataFrame(index=var))


def snapshot(mdata):
    return {
        "obs_names": list(map(str, mdata.obs_names)),
        "var_names": list(map(str, mdata.var_names)),
        "obsmap": {k: np.asarray(v).astype(int).tolist() for k, v in mdata.obsmap.items()},
        "varmap": {k: np.asarray(v).astype(int).tolist() for k, v in mdata.varmap.items()},
        "obs_columns": {c: [None if pd.isna(x) else x for x in mdata.obs[c].tolist()] for c in mdata.obs.columns},
    }


CASES = [
    {"name": "axis0_overlap", "axis": 0, "mods": {"rna": [["a", "b", "c"], ["g1", "g2"]], "prot": [["b", "c", "d"], ["p1"]]}},
    {"name": "axis0_disjoint", "axis": 0, "mods": {"rna": [["a", "b"], ["g1"]], "prot": [["c", "d"], ["p1", "p2"]]}},
    {
        "name": "axis0_three",
        "axis": 0,
        "mods": {"rna": [["a", "b", "c"], ["g1", "g2"]], "atac": [["c", "e"], ["r1"]], "prot": [["e", "a"], ["p1"]]},
    },
    {"name": "axis0_shared_feature_names", "axis": 0, "mods": {"rna": [["a", "b"], ["g1", "x"]], "prot": [["b", "c"], ["x", "p"]]}},
    {"name": "axis0_empty_vars", "axis": 0, "mods": {"gex": [["a", "b"], ["g1"]], "airr": [["b", "c"], []]}},
    {"name": "axis1_overlap", "axis": 1, "mods": {"rna": [["a", "b"], ["g1", "g2", "g3"]], "prot": [["c", "d", "e"], ["g2", "g3", "g4"]]}},
    {"name": "axism1_union", "axis": -1, "mods": {"x": [["a", "b"], ["g1", "g2"]], "y": [["b", "c"], ["g2", "g3"]]}},
    {
        "name": "axis0_filter_reorder_update",
        "axis": 0,
        "mods": {"rna": [["a", "b", "c", "d"], ["g1"]], "prot": [["b", "c"], ["p1"]]},
        "obs_column": {"score": [1.0, 2.0, 3.0, 4.0]},
        "replace": {"rna": [["c", "a", "b"], ["g1"]]},
    },
    {
        "name": "axis0_add_modality_update",
        "axis": 0,
        "mods": {"rna": [["a", "b"], ["g1"]]},
        "obs_column": {"score": [1.0, 2.0]},
        "replace": {"prot": [["b", "z"], ["p1"]]},
    },
]


def main(out):
    results = []
    for case in CASES:
        mods = {k: adata(*v) for k, v in case["mods"].items()}
        mdata = md.MuData(mods, axis=case["axis"])
        initial = snapshot(mdata)
        final = None
        if "obs_column" in case:
            for col, values in case["obs_column"].items():
                mdata.obs[col] = values
        if "replace" in case:
            for k, v in case["replace"].items():
                mdata.mod[k] = adata(*v)
            mdata.update()
            final = snapshot(mdata)
        results.append({**case, "initial": initial, "final": final})
    payload = {"versions": {p: version(p) for p in ["mudata", "anndata", "pandas"]}, "cases": results}
    with open(out, "w") as f:
        json.dump(payload, f, indent=1)
    print(f"wrote {len(results)} cases to {out}")


if __name__ == "__main__":
    main(sys.argv[1])
