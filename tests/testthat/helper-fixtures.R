fixture <- function(name) {
  path <- system.file("extdata", name, package = "mudataR")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "extdata", name)
  }
  path
}

# Zarr fixtures ship zipped (Zarr v2 metadata files are hidden files, which
# R CMD check flags); unpack each once per session
.fixture_cache <- new.env()
fixture_zarr <- function(name) {
  if (is.null(.fixture_cache[[name]])) {
    dir <- tempfile("mudatar-fixture-")
    dir.create(dir)
    utils::unzip(fixture(paste0(name, ".zip")), exdir = dir)
    .fixture_cache[[name]] <- file.path(dir, name)
  }
  .fixture_cache[[name]]
}

scirpy_expected <- function() {
  jsonlite::fromJSON(fixture("expected_airr.json"), simplifyVector = FALSE)
}

make_adata <- function(obs, var, ...) {
  anndataR::AnnData(
    obs = data.frame(row.names = obs %||% character(0)),
    var = data.frame(row.names = var %||% character(0)),
    ...
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

skip_if_no_hdf5 <- function() {
  testthat::skip_if_not_installed("rhdf5")
}

skip_if_no_zarr <- function() {
  testthat::skip_if_not_installed("Rarr")
}

skip_if_no_awkward <- function() {
  testthat::skip_if_not_installed("nanoarrow")
  testthat::skip_if_not_installed("bit64")
  testthat::skip_if_not_installed("jsonlite")
}

# Ragged values as nested R lists in the style of Python's ak.to_list(), via
# nanoarrow (public API only)
awk_rows <- function(x) {
  df <- as.data.frame(x, unnest = TRUE)
  split(df[setdiff(names(df), ".row")], factor(df$.row, levels = seq_along(x)))
}
