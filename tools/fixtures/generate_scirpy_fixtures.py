"""Generate deterministic scirpy / MuData oracle fixtures for mudataR and anndataR.

Development-only. The R packages never import Python at runtime; these files are
committed and read natively by the R test suites.

Usage (from tools/fixtures):
    uv run python -I generate_scirpy_fixtures.py <outdir>

Writes into <outdir>:
    airr.h5ad, airr_v2.zarr, airr_v3.zarr       standalone AIRR AnnData
    gex_airr.h5mu, gex_airr_v2.zarr, gex_airr_v3.zarr   paired GEX + AIRR MuData
    expected_airr.json                          expected ragged values (Python view)
    manifest.json                               versions + awkward form JSON per element
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
import mudata as md
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scirpy as ir

warnings.filterwarnings("ignore")

SEED = 20260914


def chain(locus: str, **fields) -> dict:
    c = ir.io.AirrCell.empty_chain_dict()
    c.update(
        locus=locus,
        productive=True,
        consensus_count=1,
        duplicate_count=1,
        umi_count=1,
    )
    c.update(fields)
    return c


def make_airr_cells() -> list[ir.io.AirrCell]:
    specs: list[tuple[str, list[dict]]] = [
        # zero chains: the cell exists in the AIRR modality but has no chains
        ("c0", []),
        # single VDJ chain
        ("c1", [chain("TRB", junction_aa="CASSLGQAYEQYF", v_call="TRBV7-9", umi_count=12, my_score=0.5)]),
        # paired alpha/beta
        (
            "c2",
            [
                chain("TRA", junction_aa="CAVRDSNYQLIW", v_call="TRAV1-2", umi_count=7, my_score=1.25),
                chain("TRB", junction_aa="CASSPGTGELFF", v_call="TRBV20-1", umi_count=30, my_score=-2.0),
            ],
        ),
        # multichain: >2 of each
        (
            "c3",
            [
                chain("TRA", junction_aa="CAAA", v_call="TRAV12-1", umi_count=3),
                chain("TRA", junction_aa="CAAB", v_call="TRAV12-2", umi_count=2),
                chain("TRA", junction_aa="CAAC", v_call="TRAV12-3", umi_count=1),
                chain("TRB", junction_aa="CBBA", v_call="TRBV5-1", umi_count=9),
                chain("TRB", junction_aa="CBBB", v_call="TRBV5-2", umi_count=8),
            ],
        ),
        # BCR heavy + kappa
        (
            "c4",
            [
                chain("IGH", junction_aa="CARDYW", v_call="IGHV3-23", c_call="IGHM", umi_count=40),
                chain("IGK", junction_aa="CQQYNSYPLTF", v_call="IGKV1-39", umi_count=22),
            ],
        ),
        # missing (None) and empty-string fields, non-productive chain
        (
            "c5",
            [
                chain("TRB", junction_aa=None, v_call="", productive=False, umi_count=0, my_score=None),
                chain("TRA", junction_aa="", v_call=None, umi_count=5),
            ],
        ),
        # single lambda light chain, non-ASCII-free but long string
        ("c6", [chain("IGL", junction_aa="CQSYDSSLSGSVF", v_call="IGLV1-40", umi_count=4, my_score=3.0)]),
    ]
    cells = []
    for cell_id, chains in specs:
        cell = ir.io.AirrCell(cell_id)
        for c in chains:
            cell.add_chain(c)
        cells.append(cell)
    return cells


def make_gex(obs_names: list[str]) -> ad.AnnData:
    rng = np.random.default_rng(SEED)
    n_obs, n_vars = len(obs_names), 5
    X = sp.random(n_obs, n_vars, density=0.5, format="csr", random_state=rng, dtype=np.float32)
    obs = pd.DataFrame(
        {"sample": pd.Categorical(["s1", "s2"] * (n_obs // 2) + ["s1"] * (n_obs % 2), categories=["s2", "s1"], ordered=False)},
        index=obs_names,
    )
    var = pd.DataFrame(index=[f"gene{i}" for i in range(n_vars)])
    return ad.AnnData(X=X, obs=obs, var=var)


def ak_to_jsonable(x: ak.Array):
    return ak.to_list(x)


def awkward_forms(store_path: Path, prefix: str) -> dict:
    """Read the on-disk form JSON for every awkward-array element under prefix/obsm."""
    out = {}
    if store_path.suffix in (".h5ad", ".h5mu"):
        with h5py.File(store_path, "r") as f:
            g = f[prefix + "obsm"] if prefix + "obsm" in f else {}
            for k in g:
                attrs = g[k].attrs
                if attrs.get("encoding-type") == "awkward-array":
                    out[prefix + "obsm/" + k] = {
                        "encoding-version": attrs["encoding-version"],
                        "length": int(attrs["length"]),
                        "form": json.loads(attrs["form"]),
                        "buffers": sorted(g[k].keys()),
                    }
    else:
        import zarr

        root = zarr.open_group(str(store_path), mode="r")
        g = root[prefix + "obsm"] if prefix + "obsm" in root else None
        if g is not None:
            for k in g.keys():
                attrs = dict(g[k].attrs)
                if attrs.get("encoding-type") == "awkward-array":
                    form = attrs["form"]
                    out[prefix + "obsm/" + k] = {
                        "encoding-version": attrs["encoding-version"],
                        "length": int(attrs["length"]),
                        "form": json.loads(form) if isinstance(form, str) else form,
                        "buffers": sorted(g[k].keys()),
                    }
    return out


def main(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for p in outdir.glob("*"):
        if p.name.endswith((".zarr",)):
            shutil.rmtree(p)
        elif p.suffix in (".h5ad", ".h5mu", ".json"):
            p.unlink()

    airr = ir.io.from_airr_cells(make_airr_cells())
    ir.pp.index_chains(airr)  # adds obsm["chain_indices"]

    # GEX overlaps the AIRR cells partially: c0..c4 shared, g7/g8 GEX-only, c5/c6 AIRR-only
    gex = make_gex(["c0", "c1", "c2", "c3", "c4", "g7", "g8"])
    mdata = md.MuData({"gex": gex, "airr": airr}, axis=0)
    mdata.uns["fixture"] = "scirpy-gex-airr"

    # standalone AIRR AnnData (for anndataR-level codec tests)
    airr.write_h5ad(outdir / "airr.h5ad")
    mdata.write_h5mu(outdir / "gex_airr.h5mu")
    for zf in (2, 3):
        ad.settings.zarr_write_format = zf
        airr.write_zarr(outdir / f"airr_v{zf}.zarr")
        mdata.write_zarr(outdir / f"gex_airr_v{zf}.zarr")
    ad.settings.zarr_write_format = 3

    expected = {
        "obs_names": list(airr.obs_names),
        "airr": ak_to_jsonable(airr.obsm["airr"]),
        "chain_indices": ak_to_jsonable(airr.obsm["chain_indices"]),
        "airr_type": str(airr.obsm["airr"].type),
        "chain_indices_type": str(airr.obsm["chain_indices"].type),
        "mudata": {
            "obs_names": list(mdata.obs_names),
            "var_names": list(mdata.var_names),
            "obsmap": {k: np.asarray(v).tolist() for k, v in mdata.obsmap.items()},
            "varmap": {k: np.asarray(v).tolist() for k, v in mdata.varmap.items()},
            "mod_order": list(mdata.mod.keys()),
        },
    }
    (outdir / "expected_airr.json").write_text(json.dumps(expected, indent=1))

    manifest = {
        "generator": "tools/fixtures/generate_scirpy_fixtures.py",
        "seed": SEED,
        "versions": {
            p: version(p) for p in ["mudata", "anndata", "scirpy", "awkward", "zarr", "h5py", "numpy", "pandas", "scipy"]
        },
        "forms": {
            "airr.h5ad": awkward_forms(outdir / "airr.h5ad", ""),
            "gex_airr.h5mu": awkward_forms(outdir / "gex_airr.h5mu", "mod/airr/"),
            "airr_v2.zarr": awkward_forms(outdir / "airr_v2.zarr", ""),
            "airr_v3.zarr": awkward_forms(outdir / "airr_v3.zarr", ""),
        },
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"wrote fixtures to {outdir}")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "generated" / "scirpy")
