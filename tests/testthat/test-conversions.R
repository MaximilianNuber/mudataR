skip_if_not_installed("SingleCellExperiment")
skip_if_not_installed("SummarizedExperiment")
skip_if_not_installed("S4Vectors")

conversion_mudata <- function() {
  rna <- anndataR::AnnData(
    X = matrix(as.numeric(1:6), 3, 2),
    layers = list(counts = matrix(as.numeric(11:16), 3, 2)),
    obs = data.frame(row.names = c("a", "b", "c"), ct = c("T", "B", "T")),
    var = data.frame(row.names = c("g1", "g2"))
  )
  prot <- anndataR::AnnData(
    X = matrix(as.numeric(1:3), 3, 1),
    obs = data.frame(row.names = c("c", "b", "d"), ct = c("x", "y", "z")),
    var = data.frame(row.names = "p1")
  )
  m <- MuData(list(prot = prot, rna = rna))
  obs <- m$obs
  obs$batch <- c("b1", "b1", "b2", "b2")
  obs$ct <- c("g", "g", "g", "g")
  m$obs <- obs
  m$obsm <- c(m$obsm, list(X_umap = matrix(as.numeric(1:8), 4, 2)))
  m$uns <- list(note = "global")
  m
}

test_that("as_SingleCellExperiment() puts the main modality in the main experiment and the rest in altExps", {
  m <- conversion_mudata()
  expect_message(
    sce <- m$as_SingleCellExperiment(mod_args = list(rna = list(x_mapping = "logcounts"))),
    "Keeping the 2 of 4"
  )
  expect_s4_class(sce, "SingleCellExperiment")
  # "rna" is the default main experiment even though it is not first
  expect_equal(SingleCellExperiment::altExpNames(sce), "prot")
  expect_equal(colnames(sce), c("b", "c"))
  expect_equal(SummarizedExperiment::assayNames(sce), c("logcounts", "counts"))
  expect_equal(unname(SummarizedExperiment::assay(sce, "logcounts")[1, ]), c(2, 3))
  alt <- SingleCellExperiment::altExp(sce, "prot")
  expect_equal(colnames(alt), c("b", "c"))
  # prot columns are reordered through obsmap (prot stores c, b, d)
  expect_equal(unname(SummarizedExperiment::assay(alt, "X")[1, ]), c(2, 1))
  expect_equal(alt$ct, c("y", "x"))
  # global obs: new columns added, clashing ones prefixed
  # global order is c, b, d, a (prot first), so b and c are both in batch b1
  expect_equal(sce$batch, c("b1", "b1"))
  expect_equal(sce$ct, c("B", "T"))
  expect_equal(sce[["mudata:ct"]], c("g", "g"))
  expect_equal(SingleCellExperiment::reducedDimNames(sce), "X_umap")
  meta <- S4Vectors::metadata(sce)$mudata
  expect_equal(meta$main, "rna")
  expect_equal(meta$uns$note, "global")
})

test_that("join = \"main\" keeps every main cell and pads altExps", {
  m <- conversion_mudata()
  sce <- m$as_SingleCellExperiment(main = "rna", join = "main")
  expect_equal(colnames(sce), c("a", "b", "c"))
  alt <- SingleCellExperiment::altExp(sce, "prot")
  expect_equal(alt$present, c(FALSE, TRUE, TRUE))
  expect_equal(unname(SummarizedExperiment::assay(alt, "X")[1, ]), c(NA, 2, 1))
  expect_equal(alt$ct, c(NA, "y", "x"))
})

test_that("conversion arguments are validated", {
  m <- conversion_mudata()
  expect_error(m$as_SingleCellExperiment(main = "atac"), "No modality")
  expect_error(m$as_SingleCellExperiment(mod_args = list(atac = list())), "unknown modalities")
  expect_error(m$as_SingleCellExperiment(mod_args = list(rna = list(xmapping = "a"))), "unknown arguments")
  expect_error(m$as_SingleCellExperiment(mod_args = list(rna = "counts")), "named list")
  shared_var <- MuData(
    list(a = make_adata("x", c("g1", "g2")), b = make_adata("y", c("g1", "g2"))),
    axis = 1
  )
  expect_error(shared_var$as_SingleCellExperiment(), class = "mudatar_unsupported")
  expect_error(shared_var$as_MultiAssayExperiment(), class = "mudatar_unsupported")
})

test_that("scirpy AIRR data stay available in the converted altExp", {
  skip_if_no_hdf5()
  skip_if_no_awkward()
  m <- read_h5mu(fixture("gex_airr.h5mu"), as = "HDF5MuData")
  sce <- suppressMessages(m$as_SingleCellExperiment(main = "gex"))
  alt <- SingleCellExperiment::altExp(sce, "airr")
  airr <- S4Vectors::metadata(alt)$awkward_obsm$airr
  expect_s3_class(airr, "AwkwardArrowArray")
  expect_length(airr, ncol(sce))
  expected <- scirpy_expected()
  chains <- as.data.frame(airr, unnest = TRUE)
  # the 5 shared cells c0..c4 are the first 5 AIRR cells
  n_chains <- lengths(expected$airr[1:5])
  expect_equal(nrow(chains), sum(n_chains))

  padded <- m$as_SingleCellExperiment(main = "gex", join = "main")
  alt <- SingleCellExperiment::altExp(padded, "airr")
  expect_equal(alt$present, c(rep(TRUE, 5), FALSE, FALSE))
  expect_equal(S4Vectors::metadata(alt)$awkward_obsm_columns, paste0("c", 0:4))
})

test_that("as_MultiAssayExperiment() maps modalities, obs and obsmap", {
  skip_if_not_installed("MultiAssayExperiment")
  m <- conversion_mudata()
  mae <- m$as_MultiAssayExperiment(mod_args = list(rna = list(x_mapping = "logcounts")))
  expect_s4_class(mae, "MultiAssayExperiment")
  expect_equal(names(mae), c("prot", "rna"))
  expect_s4_class(mae[["rna"]], "SingleCellExperiment")
  expect_equal(SummarizedExperiment::assayNames(mae[["rna"]]), c("logcounts", "counts"))
  # colData is the MuData obs, in MuData's global order (prot listed first)
  expect_equal(rownames(MultiAssayExperiment::colData(mae)), c("c", "b", "d", "a"))
  expect_equal(mae$batch, c("b1", "b1", "b2", "b2"))
  sm <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
  expect_equal(nrow(sm), 6L)
  prot <- sm[sm$assay == "prot", ]
  expect_equal(prot$primary, c("c", "b", "d"))
  expect_equal(prot$colname, c("c", "b", "d"))
  rna <- sm[sm$assay == "rna", ]
  expect_equal(rna$primary, c("c", "b", "a"))
  expect_equal(rna$colname, c("c", "b", "a"))
  expect_equal(S4Vectors::metadata(mae)$mudata$uns$note, "global")
  expect_true("X_umap" %in% names(S4Vectors::metadata(mae)$mudata$obsm))

  skip_if_no_hdf5()
  h5 <- read_h5mu(fixture("gex_airr.h5mu"), as = "HDF5MuData")
  mae <- h5$as_MultiAssayExperiment()
  expect_equal(ncol(mae[["airr"]]), 7L)
  expect_length(S4Vectors::metadata(mae[["airr"]])$awkward_obsm$airr, 7)
})

test_that("each MuData backend holds the matching anndataR class", {
  skip_if_no_hdf5()
  h5 <- read_h5mu(fixture("gex_airr.h5mu"), as = "HDF5MuData")
  for (mod in as.list(h5$mod)) expect_s3_class(mod, "HDF5AnnData")
  mem <- read_h5mu(fixture("gex_airr.h5mu"))
  for (mod in as.list(mem$mod)) expect_s3_class(mod, "InMemoryAnnData")
  # backed AnnData objects and views given to MuData() are materialised
  mixed <- MuData(list(gex = h5[["gex"]], view = mem[["gex"]][1:2, ]))
  for (mod in as.list(mixed$mod)) expect_s3_class(mod, "InMemoryAnnData")
  mixed[["airr"]] <- h5[["airr"]]
  expect_s3_class(mixed[["airr"]], "InMemoryAnnData")
  skip_if_no_zarr()
  z <- read_mudata_zarr(fixture_zarr("gex_airr_v3.zarr"), as = "ZarrMuData")
  for (mod in as.list(z$mod)) expect_s3_class(mod, "ZarrAnnData")
})
