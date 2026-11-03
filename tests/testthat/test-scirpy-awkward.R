skip_if_no_hdf5()
skip_if_no_awkward()

expected <- scirpy_expected()
h5mu <- fixture("gex_airr.h5mu")

test_that("scirpy AIRR data in an H5MU become nested Arrow without Python", {
  skip_if_not_installed("arrow")
  m <- read_h5mu(h5mu, as = "HDF5MuData")
  airr <- get_awkward(m[["airr"]])$airr
  expect_s3_class(airr, "AwkwardArrowArray")
  expect_length(airr, m[["airr"]]$n_obs())
  tab <- arrow::as_arrow_table(airr, name = "airr")
  expect_equal(names(tab), "airr")
  chains <- as.data.frame(airr, unnest = TRUE)
  exp_rows <- unlist(lapply(seq_along(expected$airr), function(i) rep(i, length(expected$airr[[i]]))))
  expect_equal(chains$.row, exp_rows)
  exp_locus <- unlist(lapply(expected$airr, function(cell) vapply(cell, `[[`, "", "locus")))
  expect_equal(chains$locus, exp_locus)
  exp_junction <- unlist(lapply(expected$airr, function(cell) vapply(cell, function(ch) ch$junction_aa %||% NA_character_, "")))
  expect_equal(chains$junction_aa, exp_junction)

  ci <- as.data.frame(arrow::as_arrow_table(get_awkward(m[["airr"]])$chain_indices))
  expect_equal(ci$multichain, vapply(expected$chain_indices, `[[`, TRUE, "multichain"))
})

test_that("AIRR data subset along observations keep their ragged rows", {
  m <- read_h5mu(h5mu, as = "HDF5MuData")
  airr <- get_awkward(m[["airr"]])$airr
  sub <- airr[c(4, 1)]
  rows <- as.data.frame(sub, unnest = TRUE)
  expect_equal(rows$.row, rep(1, 5))
  expect_equal(rows$junction_aa, c("CAAA", "CAAB", "CAAC", "CBBA", "CBBB"))
  expect_length(sub, 2)
})

test_that("unchanged AIRR data are written back with their original form", {
  path <- withr::local_tempfile(fileext = ".h5mu")
  src <- read_h5mu(h5mu)
  write_h5mu(src, path)
  form_src <- rhdf5::h5readAttributes(h5mu, "/mod/airr/obsm/airr")$form
  form_out <- rhdf5::h5readAttributes(path, "/mod/airr/obsm/airr")$form
  expect_identical(form_out, form_src)
  back <- read_h5mu(path, as = "HDF5MuData")
  expect_equal(
    as.data.frame(get_awkward(back[["airr"]])$airr, unnest = TRUE),
    as.data.frame(as_AwkwardArrowArray(src[["airr"]]$obsm$airr), unnest = TRUE)
  )
})

test_that("backed modalities are opened without their Awkward Arrays", {
  m <- read_h5mu(h5mu, as = "HDF5MuData")
  airr <- m[["airr"]]
  expect_s3_class(airr, "HDF5AnnData")
  expect_length(airr$obsm, 0)
  expect_named(get_awkward(airr), c("airr", "chain_indices"))
  expect_length(get_awkward(airr, "varm"), 0)
  expect_length(get_awkward(m[["gex"]]), 0)
})

test_that("in-memory modalities hold AIRR data as AwkwardFrame", {
  mem <- read_h5mu(h5mu)
  airr <- mem[["airr"]]
  expect_s3_class(airr, "InMemoryAnnData")
  frame <- airr$obsm$airr
  expect_s3_class(frame, "AwkwardFrame")
  expect_equal(nrow(frame), airr$n_obs())
  expect_s3_class(as_AwkwardArrowArray(frame), "AwkwardArrowArray")
  expect_output(print(frame), "<AwkwardFrame> 7 \\* var \\*")
  expect_equal(names(get_awkward(airr)), c("airr", "chain_indices"))
})

test_that("anndataR subsets AwkwardFrames with its standard rules", {
  airr <- read_h5mu(h5mu)[["airr"]]
  full <- as_AwkwardArrowArray(airr$obsm$airr)
  for (idx in list(c(4L, 1L), airr$obs_names %in% c("c0", "c3"), c("c3", "c0"))) {
    view <- airr[idx, ]
    expect_s3_class(view$obsm$airr, "AwkwardFrame")
    expect_equal(nrow(view$obsm$airr), view$n_obs())
    sub <- view$as_InMemoryAnnData()
    expect_s3_class(sub$obsm$airr, "AwkwardFrame")
  }
  sub <- airr[c(4, 1), ]$as_InMemoryAnnData()
  expect_equal(
    as.data.frame(as_AwkwardArrowArray(sub$obsm$airr), unnest = TRUE),
    as.data.frame(full[c(4, 1)], unnest = TRUE)
  )
  # subset modalities are written with the selected rows only
  m <- MuData(list(airr = sub))
  path <- withr::local_tempfile(fileext = ".h5mu")
  write_h5mu(m, path)
  back <- get_awkward(read_h5mu(path, as = "HDF5MuData")[["airr"]])$airr
  expect_length(back, 2)
  expect_equal(
    as.data.frame(back, unnest = TRUE)$junction_aa,
    c("CAAA", "CAAB", "CAAC", "CBBA", "CBBB")
  )
})

test_that("an anndataR SingleCellExperiment conversion skips AwkwardFrames", {
  skip_if_not_installed("SingleCellExperiment")
  airr <- read_h5mu(h5mu)[["airr"]]
  sce <- airr$as_SingleCellExperiment()
  expect_false(any(c("airr", "chain_indices") %in% SingleCellExperiment::reducedDimNames(sce)))
})

test_that("AwkwardFrames are created from Awkward and Arrow data", {
  skip_if_not_installed("arrow")
  awk <- as_AwkwardArrowArray(nanoarrow::as_nanoarrow_array(data.frame(a = 1:3, b = c("x", "y", "z"))))
  frame <- as_AwkwardFrame(awk)
  expect_true(is_AwkwardFrame(frame))
  expect_identical(as_AwkwardFrame(frame), frame)
  expect_identical(as_AwkwardArrowArray(frame), awk)
  adata <- anndataR::AnnData(
    obs = data.frame(row.names = c("c1", "c2", "c3")),
    var = data.frame(row.names = "g1"),
    obsm = list(chains = frame)
  )
  sub <- adata[c(3, 1), ]$obsm$chains
  expect_equal(as.data.frame(arrow::as_arrow_table(sub))$b, c("z", "x"))
  expect_error(frame[1, 1], "rows only")
  expect_error(as_AwkwardArrowArray(frame[c(1, NA), ]), class = "mudatar_unsupported")
})

test_that("Zarr stores with AIRR data round trip, keeping the form", {
  skip_if_no_zarr()
  for (zarr_format in c(2L, 3L)) {
    src <- read_mudata_zarr(fixture_zarr(sprintf("gex_airr_v%d.zarr", zarr_format)))
    expect_s3_class(src[["airr"]]$obsm$airr, "AwkwardFrame")
    path <- withr::local_tempfile(fileext = ".zarr")
    write_mudata_zarr(src, path, zarr_format = zarr_format)
    form_src <- Rarr::read_zarr_attributes(file.path(fixture_zarr(sprintf("gex_airr_v%d.zarr", zarr_format)), "mod", "airr", "obsm", "airr"))$form
    form_out <- Rarr::read_zarr_attributes(file.path(path, "mod", "airr", "obsm", "airr"))$form
    expect_identical(form_out, form_src)
    z <- read_mudata_zarr(path, as = "ZarrMuData")
    expect_equal(
      as.data.frame(get_awkward(z[["airr"]])$airr, unnest = TRUE),
      as.data.frame(as_AwkwardArrowArray(src[["airr"]]$obsm$airr), unnest = TRUE)
    )
  }
})
