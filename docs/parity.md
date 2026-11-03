# mudataR ↔ Python mudata parity

This is the contract for what "parity" means in mudataR. Every public member of
`mudata.MuData` and every public module-level function of the pinned Python
version has a row. A row is **complete** only when it works natively (no
Python) for every backend listed and is covered by tests.

## Targets

| Component | Version / commit | Source |
| --- | --- | --- |
| Python `mudata` (reference) | 0.4.1 | <https://mudata.readthedocs.io/stable/generated/mudata.MuData.html> |
| MuData on-disk format | RFC, `encoding-version` 0.1.0 | <https://mudata.readthedocs.io/stable/io/spec.html> |
| Python `anndata` | 0.13.4 | <https://anndata.readthedocs.io/en/stable/fileformat-prose.html> |
| Python `awkward` | 2.14.0 | `awkward-array` encoding 0.1.0 |
| Python `scirpy` | 0.25.1 | <https://scirpy.readthedocs.io/en/latest/data-structure.html> |
| `anndataR` | ≥ 1.3.3 (`devel` @ `577fdc5`), unmodified | <https://github.com/scverse/anndataR> |

The member list was produced by introspecting `mudata.MuData` 0.4.1
(`inspect.getattr_static`); re-run that introspection when the target version
moves (see `tools/fixtures/`).

Status legend: **complete**, **equivalent** (same behaviour, R naming or
idiom), **partial** (what works is listed), **planned**, **deprecated
upstream**, **unsupported** (intentionally, with reason).

Backends: M = `InMemoryMuData`, H = `HDF5MuData`, Z = `ZarrMuData`,
R = `ReticulateMuData`, V = `MuDataView`.

## `mudata.MuData` properties

| Python | mudataR | Status | Notes |
| --- | --- | --- | --- |
| `axis` | `$axis` | complete (M, H, Z) | Read-only after construction. Missing `axis` in old files reads as `0`. |
| `filename` | `$filename` | partial (M, H, Z) | Read-only. Python allows assigning to re-bind/move the backing file; mudataR errors and points to `write_h5mu()`/`read_h5mu()`. |
| `is_view` | `$is_view` | partial | Always `FALSE` until `MuDataView` exists. |
| `isbacked` | `$isbacked` | complete (M, H, Z) | |
| `mod` | `$mod` | complete (M, H, Z) | A guarded `ModalityMapping`: `names()`, `length()`, `[[`, `$`, `as.list()`; assignment validates and rebuilds maps. Returns `InMemoryAnnData` (M), read-only `HDF5AnnData` (H) or `ZarrAnnData` (Z); see *Modality classes per backend*. |
| `mod_names` | `$mod_names()` | equivalent | Method rather than property (anndataR convention). |
| `n_mod` | `$n_mod()` | equivalent | |
| `n_obs` | `$n_obs()`, `nrow()` | equivalent | |
| `n_vars` / `n_var` | `$n_vars()`, `ncol()` | equivalent | `n_var` is an alias in Python. |
| `obs` | `$obs` | complete (M, H, Z) | Row names are added on access when `obs_names` are unique. |
| `obs_names` | `$obs_names`, `rownames()` | partial | Settable on M (same length); read-only on H/Z. |
| `obsm` | `$obsm` | complete (M, H, Z) | Includes the per-modality presence masks Python stores (`obsm[[mod]]`). |
| `obsmap` | `$obsmap` | complete (M, H, Z) | Integer maps, `0` = absent. Written as `uint32` like Python. |
| `obsp` | `$obsp` | complete (M, H, Z) | |
| `shape` | `$shape()`, `dim()` | equivalent | |
| `uns` | `$uns` | complete (M, H, Z) | Any element encoding supported by anndataR. |
| `var` | `$var` | complete (M, H, Z) | |
| `var_names` | `$var_names`, `colnames()` | partial | As `obs_names`. |
| `varm` | `$varm` | complete (M, H, Z) | |
| `varmap` | `$varmap` | complete (M, H, Z) | |
| `varp` | `$varp` | complete (M, H, Z) | |

## `mudata.MuData` methods

| Python | mudataR | Status | Notes |
| --- | --- | --- | --- |
| `MuData(data, ..., axis)` | `MuData()` / `InMemoryMuData$new()` | partial | Named list of anndataR AnnData objects; supplied maps are validated, not recomputed. Not yet: splitting a single AnnData by `feature_types`, nested MuData modalities. |
| `copy(filename=None)` | `$copy()` | partial | In-memory copy; `filename` (copy to a backed file) is planned. |
| `obs_keys()` / `var_keys()` | `$obs_keys()` / `$var_keys()` | complete | |
| `obsm_keys()` / `varm_keys()` | `$obsm_keys()` / `$varm_keys()` | complete | Also `$obsp_keys()`, `$varp_keys()`. |
| `uns_keys()` | `$uns_keys()` | complete | |
| `update()` | `$update()` | partial | Port of `_update_attr` main case, differential-tested against Python (`test-update-differential.R`): axes 0/1/-1, overlapping, disjoint and duplicated-across-modality names, reorder-after-filter, new modalities, global column reindexing, `obsm[[mod]]` masks. Not yet: duplicated names *within* a modality (errors with `mudatar_unsupported`; anndataR cannot hold such names, scverse/anndataR#256). |
| `update_obs()` / `update_var()` | `$update_obs()` / `$update_var()` | partial | As `update()`. Like mudata ≥ 0.3 they never pull modality columns. |
| `pull_obs()` / `pull_var()` | `$pull_obs()` / `$pull_var()` | partial (R) | Available on `ReticulateMuData` (runs Python, full signature). Native backends: planned, with the full signature (`columns`, `mods`, `common`, `join_common`, `nonunique`, `join_nonunique`, `unique`, `prefix_unique`, `drop`, `only_drop`). |
| `push_obs()` / `push_var()` | `$push_obs()` / `$push_var()` | partial (R) | Available on `ReticulateMuData` (runs Python, full signature). Native backends: planned, with the full signature (`columns`, `mods`, `common`, `prefixed`, `drop`, `only_drop`). |
| `obs_names_make_unique()` / `var_names_make_unique()` | — | planned | |
| `strings_to_categoricals(df)` | — | planned | |
| `to_anndata()` | — | planned | Also module-level `mudata.to_anndata()`. |
| `obs_vector()` / `var_vector()` | — | deprecated upstream | Deprecated in mudata 0.3.4; not planned. |
| `write()` / `write_h5mu()` | `$write_h5mu()`, `write_h5mu()` | complete (M, H, Z sources) | Native; temp sibling + rename; `mode = "w-"` default. |
| `write_zarr()` | `$write_zarr()`, `write_mudata_zarr()` | complete (M, H, Z sources) | Zarr v2 and v3. |
| `__getitem__(str)` / `__setitem__` | `[[`, `[[<-` | complete (M); partial (H, Z) | H/Z (`mode = "r+"`) can add or replace modalities; removing one is planned. |
| `__getitem__(index)` (views) | `[` | planned | `MuDataView` (lazy, read-only, zero-sentinel translation). `[` currently errors with `mudatar_unsupported`. |
| `__repr__` | `print()`, `format()` | equivalent | |

## R-only conversions

| mudataR | Status | Notes |
| --- | --- | --- |
| `$as_SingleCellExperiment(main, mod_args, join)` | complete (M, H, Z) | `main` (default `"rna"`) is the main experiment, other modalities are `altExp()`s. Each modality goes through `anndataR::as_SingleCellExperiment()` with its `mod_args` entry. `join = "intersect"` (shared cells) or `"main"` (all main cells, altExps padded with `NA` plus a `present` column). Global `obs` → main `colData` (`mudata:` prefix on clashes), global `obsm` matrices → main `reducedDims`, global `uns`/`obsm` → `metadata(sce)$mudata`, Awkward Arrays → `metadata(<exp>)$awkward_obsm`. Requires `axis` 0 or -1. |
| `$as_MultiAssayExperiment(mod_args)` | complete (M, H, Z) | Modalities → experiments (via anndataR), global `obs` → `colData`, `obsmap` → `sampleMap`; global `uns`/`obsm`/`obsp` → `metadata(mae)$mudata`. Requires `axis` 0 or -1. |

## Modality classes per backend

| MuData backend | Modality class |
| --- | --- |
| `InMemoryMuData` | always `anndataR::InMemoryAnnData` (backed, Python-backed and view AnnData objects are materialised on construction and on `[[<-`) |
| `HDF5MuData` | read-only `anndataR::HDF5AnnData` reading `/mod/<name>` through a temporary H5AD of HDF5 external links (released anndataR cannot open a non-root group); replace a modality with `mdata[[name]] <- adata`. Awkward Arrays are left out; use `get_awkward()`. |
| `ZarrMuData` | `anndataR::ZarrAnnData` opened at `<store>/mod/<name>` with the MuData's `mode`; read-only (through a symlink shim) when the modality holds Awkward Arrays, which `get_awkward()` returns. |
| `ReticulateMuData` | `anndataR::ReticulateAnnData` (nested MuData modalities as `ReticulateMuData`) |

## Module-level functions

| Python | mudataR | Status | Notes |
| --- | --- | --- | --- |
| `read_h5mu(filename, backed)` | `read_h5mu(path, as, mode, backed, modalities)` | complete | `as = "InMemoryMuData"` or `"HDF5MuData"`; `backed = TRUE` gives `DelayedArray` matrices and lazily loaded Awkward buffers (`get_awkward()`); `modalities` exposes a subset without touching the rest. |
| `read_zarr(store)` | `read_mudata_zarr(path, as, mode, modalities)` | complete | Named to avoid masking `anndataR::read_zarr()`. Zarr v2 and v3. |
| `read(path)` | — | planned | Extension-dispatching reader. |
| `write_h5mu(filename, mdata)` | `write_h5mu(x, path, mode, compression)` | complete | Argument order follows R (object first). Writes the 512-byte user block. |
| `write_zarr(store, data)` | `write_mudata_zarr(x, path, mode, compression, zarr_format)` | complete | |
| `write(path, data)` | — | planned | |
| `read_h5ad(filename, mod)` / `write_h5ad(filename, mod, data)` | `read_h5mu(as = "HDF5MuData")[[mod]]`, `$set_mod()` on `mode = "r+"` | equivalent | Single-modality access to an H5MU file. |
| `read_anndata` / `write_anndata` | `anndataR::read_h5ad()` / `anndataR::write_h5ad()` | equivalent | Owned by anndataR. |
| `concat()` | — | planned | |
| `to_anndata()` / `to_mudata()` | — | planned | |
| `register_mudata_namespace`, `ExtensionNamespace`, `acc` | — | unsupported | Python accessor machinery; R uses ordinary functions. |
| `settings`, `suppress`, `utils` | — | unsupported | Python-specific. `pull_on_update` is always `FALSE` (the mudata ≥ 0.4 default). |

## Backends

| Backend | Status | Notes |
| --- | --- | --- |
| `InMemoryMuData` | partial | See tables above. |
| `HDF5MuData` | partial | Lazy: construction reads root attributes only. Modalities are read-only `HDF5AnnData` objects (see above); with `mode = "r+"` global slots are written in place and modalities are replaced with `[[<-`. Lifecycle: closing a modality never affects the parent or siblings; closing the parent stops access through it, but modality objects obtained earlier stay usable (released anndataR has no hook to invalidate them). |
| `ZarrMuData` | partial | As `HDF5MuData`, with `ZarrAnnData` modalities that are editable in place when they hold no Awkward Arrays. After in-place edits a Python-written consolidated metadata file (`.zmetadata`) is not updated. `backed` is not available (anndataR's Zarr backed mode waits on ZarrArray, scverse/anndataR#507). |
| `ReticulateMuData` | complete | Refactor of the MofaflexR prototype onto `AbstractMuData`. Wraps a Python `mudata.MuData` (`convert = FALSE`); modalities are `anndataR::ReticulateAnnData`. `update()`, `pull_*()`, `push_*()` and `copy()` run in Python; `py_mudata()` and `py_call()` are kept. `as_ReticulateMuData()` converts any MuData. Python is optional and never needed by native code. Tested against mudata 0.4.1 / anndata 0.13.4 when `MUDATAR_ORACLE_PYTHON` is set. |
| `MuDataView` | planned | |

## On-disk format

| Feature | Status | Notes |
| --- | --- | --- |
| Root `encoding-type = "MuData"`, `encoding-version = "0.1.0"`, `axis`, `encoder`, `encoder-version` | complete | |
| `/mod` `mod-order`, honoured only when every listed modality exists (RFC) | complete | Always written as an array, also for one modality. |
| `/mod/<name>` as AnnData `0.1.0` groups | complete | Read through anndataR via shims (HDF5 external links / Zarr sub-paths); written by anndataR to a temporary H5AD copied in with `H5Ocopy` (HDF5) or directly to `<store>/mod/<name>` (Zarr). |
| `obsmap`/`varmap` as `uint32` arrays | complete | |
| H5MU 512-byte user block with `MuData (format-version=…)` marker | complete | Created with `H5Pset_userblock()` (no compiled code needed); verified byte-level and by Python mudata reading without warning. |
| Unknown element encodings | complete | anndataR's "could not read/write element" warnings are raised as `mudatar_anndatar_error` errors, so elements are never dropped silently. |
| Safe writes | complete | Validation before writing; temp sibling + rename. Directories (Zarr) are removed then replaced, which is not atomic. |
| Duplicated global names on write | unsupported (for now) | Errors; needs anndataR data frames with duplicate indices (scverse/anndataR#256). |

## Awkward Arrays (`awkward-array` 0.1.0, in mudataR)

anndataR does not handle this encoding; mudataR reads and writes it itself and never shows it to anndataR's readers or writers. In in-memory modalities aligned Awkward Arrays are `AwkwardFrame`s, data frames of row positions that anndataR subsets with its standard rules; `get_awkward()` returns them (also for backed modalities) as `AwkwardArrowArray`s.


| Node / feature | Read | Write | Notes |
| --- | --- | --- | --- |
| `NumpyArray` bool, (u)int8/16/32/64, float32/64 | yes | yes | int64 values are exact up to 2^53 in R; larger values are kept as `integer64` when read from HDF5 and refused (not rounded) when writing. |
| `NumpyArray` with `inner_shape` | yes | via `RegularArray` | |
| `NumpyArray` datetime64/timedelta64, float16, complex | no | no | Explicit `mudatar_unsupported_encoding` error. |
| `EmptyArray` | yes | yes | |
| `ListOffsetArray` (i32/u32/i64) | yes | yes | |
| `ListArray` | yes | as `ListOffsetArray` | Python's writer never stores it (it packs first); tested with R-built forms. |
| `RegularArray` | yes | yes | |
| `RecordArray` (records and tuples) | yes | yes | Field order preserved. |
| `IndexedArray` | yes | yes | |
| `IndexedOptionArray` | yes | yes | |
| `ByteMaskedArray` | yes | yes | |
| `BitMaskedArray` | yes | as `ByteMaskedArray` | Unmodified arrays are written back with their original form. |
| `UnmaskedArray` | yes | yes | |
| `UnionArray` | no | no | Explicit error naming the element and form node (fixture `union.h5ad`). |
| `__array__ = "string"` / `"bytestring"` | yes | yes | Arrow `large_string`/`string`, `large_binary`/`binary`. |
| Arrow conversion | yes | yes | Option nodes ↔ nullable fields; non-option ↔ non-nullable, so `?string` vs `string` survives Arrow. |
| Lossless pass-through | — | yes | Unmodified arrays keep form and buffers. |
| Lazy backed reads (H5MU) | yes | — | Buffers load on first use. Row-range reads of buffers are planned. |
| Zarr int64 buffers beyond the int32 range | read: as Rarr returns them | refused | Rarr writes `<i8` only from R integers. |
| Python `awkward.Array` via reticulate | yes | yes | `r_to_py()` / `py_to_r()` with `ak.from_buffers()` / `ak.to_buffers()`, lossless (int64 exact to 2^53). |
| scirpy `obsm["airr"]`, `obsm["chain_indices"]` | yes | yes | H5AD, H5MU, Zarr v2/v3; verified with Python anndata, mudata and `scirpy.tl.chain_qc`. |
