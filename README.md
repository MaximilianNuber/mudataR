# mudataR

Native [MuData](https://mudata.readthedocs.io) for R: read and write H5MU
files and MuData Zarr stores **without Python**, with an
[anndataR](https://github.com/scverse/anndataR)-style backend API and
lossless support for [scirpy](https://scirpy.scverse.org) immune repertoire
data (Awkward Arrays) as nested Arrow.

> Status: early development (0.0.0.9000). mudataR works with **unmodified
> anndataR** (>= 1.3.3, currently Bioconductor devel) through its public API
> only. Everything anndataR does not offer yet (MuData's nested AnnData
> groups, single-element I/O, Awkward Arrays) is handled inside mudataR; see
> the [roadmap](#roadmap) for what that means today and how it is meant to
> change. [`docs/parity.md`](docs/parity.md) lists exactly what is
> implemented.

## Installation

```r
# install.packages("BiocManager")
BiocManager::install(c("scverse/anndataR", "Huber-group-EMBL/Rarr@devel", "rhdf5"))
BiocManager::install("MaximilianNuber/mudataR")
# For Awkward Arrays / Arrow: install.packages(c("nanoarrow", "bit64", "arrow"))
```

## Example

```r
library(mudataR)

# Read a Python-written H5MU into memory (no Python involved)
mdata <- read_h5mu("gex_airr.h5mu")
mdata
mdata$obsmap

# Modalities are anndataR AnnData objects; scirpy's AIRR data are an
# AwkwardFrame that anndataR subsets with its standard rules
airr <- mdata[["airr"]]
airr$obsm$airr                       # <AwkwardFrame> 7 * var * {...}
sub <- airr[c(2, 4), ]               # an anndataR view, AIRR rows included

# The Awkward data as nested Arrow, or one row per chain
awk <- as_AwkwardArrowArray(sub$obsm$airr)
arrow::as_arrow_table(awk, name = "airr")
chains <- as.data.frame(awk, unnest = TRUE)

# Open lazily instead; Awkward Arrays of backed modalities via get_awkward()
backed <- read_h5mu("gex_airr.h5mu", as = "HDF5MuData")
get_awkward(backed[["airr"]])$airr   # <AwkwardArrowArray> 7 * var * {...}

# Write back; Python mudata/scirpy read the result
write_h5mu(mdata, "out.h5mu")
write_mudata_zarr(mdata, "out.zarr", zarr_format = 3)

# Build a MuData in R
rna <- anndataR::AnnData(
  obs = data.frame(row.names = c("c1", "c2", "c3")),
  var = data.frame(row.names = c("g1", "g2"))
)
prot <- anndataR::AnnData(
  obs = data.frame(row.names = c("c2", "c3", "c4")),
  var = data.frame(row.names = "p1")
)
m <- MuData(list(rna = rna, prot = prot))
m$obsmap
```

## Roadmap

mudataR is meant to be a thin MuData layer on top of anndataR. Today it
bridges three gaps in released anndataR itself; each gap closes when the
matching (small, general-purpose) feature lands in anndataR, without changing
mudataR's user-facing API.

### Awkward Arrays (scirpy)

| | Now | Later |
| --- | --- | --- |
| Reading and writing `awkward-array` elements | mudataR's own codec (rhdf5 / Rarr), lossless; unchanged arrays keep their original form and buffers | Same codec, plugged into anndataR's readers and writers through an encoding-registration hook, so anndataR objects can hold Awkward Arrays directly |
| In-memory modalities | `obsm`/`varm` hold an `AwkwardFrame`: a data frame of row positions carrying the `AwkwardArrowArray`. anndataR subsets it with its normal data-frame rules (views, `adata[i, ]`, `as_InMemoryAnnData()`) | Same, or the `AwkwardArrowArray` itself once anndataR views dispatch subsetting to a generic instead of passing unknown classes through unchanged |
| Backed modalities (`HDF5MuData`, `ZarrMuData`) | anndataR opens them **without** their Awkward Arrays; `get_awkward(mdata[["airr"]])` returns them (lazily with `backed = TRUE`) | `mdata[["airr"]]$obsm$airr` directly, via the registration hook |
| Arrow / data frames | `as_AwkwardArrowArray()`, `arrow::as_arrow_table()`, `as.data.frame(unnest = TRUE)` | Unchanged |
| Python (`ReticulateMuData`) | `awkward.Array` ↔ `AwkwardArrowArray` via `ak.to_buffers()`/`ak.from_buffers()`; R ↔ Python AnnData goes through a temporary H5AD | Direct object conversion once anndataR's `ReticulateAnnData` handles Awkward Arrays, `layers` with anndata ≥ 0.13 and zero-length axes |
| SingleCellExperiment / MultiAssayExperiment | Kept in `metadata(<experiment>)$awkward_obsm` | Unchanged |

### Loading and writing

| | Now | Later |
| --- | --- | --- |
| H5MU modalities (`/mod/<name>`) | Read through a temporary H5AD whose entries are HDF5 external links into the H5MU; written by anndataR to a temporary H5AD that is copied in with `H5Ocopy` | anndataR opens and writes `root = "/mod/<name>"` directly (rooted stores) |
| Zarr modalities | anndataR opens `<store>/mod/<name>` directly (a Zarr group is a directory) | Unchanged |
| Global slots (`obs`, `obsm`, `uns`, …) | anndataR's decoders/encoders on a single element through a temporary AnnData (`uns/value`) | anndataR's `read_elem()`/`write_elem()` |
| `obsmap`/`varmap`, user block, root attributes | mudataR (rhdf5 / Rarr) | Unchanged |
| Elements anndataR cannot read or write | Raised as errors (`mudatar_anndatar_error`), never silently dropped | Unchanged |
| Cost | Writing a modality rewrites the whole modality, including `X` | Elements are written individually |

### Backed MuData objects

| | Now | Later |
| --- | --- | --- |
| `HDF5MuData` modalities | Read-only `HDF5AnnData` objects (`backed = TRUE` gives `DelayedArray`s); change a modality by replacing it: `mdata[["rna"]] <- adata` | Writable `HDF5AnnData` objects rooted in the H5MU (`mdata[["rna"]]$obs <- ...` writes in place) |
| `ZarrMuData` modalities | `ZarrAnnData` objects editable in place (with `mode = "r+"`); read-only when the modality holds Awkward Arrays | Always editable in place |
| Global slots with `mode = "r+"` | Written in place | Unchanged |
| Closing | `mdata$close()` stops access through the MuData; modality objects obtained earlier stay usable | Closing the MuData invalidates its modalities with a clear error |
| Python-written consolidated Zarr metadata (`.zmetadata`) | Not updated by in-place edits | Updated or removed on write |
| Lazy `MuDataView` (`mdata[i, j]`) | Planned | Planned |

### What closes the gaps in anndataR

Each of these is useful to anndataR on its own and is tracked for upstream
contribution:

1. **Rooted stores**: `read_h5ad(root =)`, `write_h5ad(root =)` (and Zarr),
   plus `read_elem()`/`write_elem()`. Removes the HDF5 shims and makes H5MU
   modalities writable in place.
2. **Extension points**: let view subsetting dispatch to a generic for
   classes anndataR does not know, and let packages register readers and
   writers for extra encodings such as `awkward-array`
   ([scverse/anndataR#404](https://github.com/scverse/anndataR/issues/404)).
3. **`ReticulateAnnData` fixes**: `layers` with Python anndata ≥ 0.13 (the
   `None` key for `X`) and AnnData objects with zero observations or
   variables.

## Design

* `AbstractMuData` → `InMemoryMuData`, `HDF5MuData`, `ZarrMuData`,
  `ReticulateMuData` (R6, active bindings for slots, S3 `dim()`, `[[`,
  `[[<-`, `print()`).
* Modalities are real anndataR objects of the matching backend
  (`InMemoryAnnData`, `HDF5AnnData`, `ZarrAnnData`, `ReticulateAnnData`).
* Element encodings are anndataR's; mudataR adds the MuData layer (root and
  `mod` attributes, `uint32` maps, the 512-byte H5MU user block) and the
  Awkward Array codec. `R/anndatar-compat.R` is the only place that works
  around anndataR's current API.
* Global indices and maps follow Python mudata 0.4.1, checked by
  differential tests against Python-generated oracle cases.
* Writes validate the whole object first and go through a temporary sibling
  path.

## Development

Python is used only to create fixtures and to independently verify R output
(`tools/fixtures/`, pinned with `uv`). Native tests never start Python;
cross-language and `ReticulateMuData` tests run when `MUDATAR_ORACLE_PYTHON`
points at the fixture environment:

```sh
(cd tools/fixtures && uv sync)
MUDATAR_ORACLE_PYTHON=$PWD/tools/fixtures/.venv/bin/python \
  Rscript -e 'devtools::test()'
```
