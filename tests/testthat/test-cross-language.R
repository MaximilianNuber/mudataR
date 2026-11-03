# Optional oracle job: Python mudata/scirpy read files written by R.
# Enabled by pointing MUDATAR_ORACLE_PYTHON at the interpreter of the
# environment in tools/fixtures (never needed for native use).

oracle_python <- function() {
  py <- Sys.getenv("MUDATAR_ORACLE_PYTHON")
  if (!nzchar(py) || !file.exists(py)) {
    skip("MUDATAR_ORACLE_PYTHON is not set")
  }
  py
}

verifier <- function() {
  path <- test_path("..", "..", "tools", "fixtures", "verify_r_outputs.py")
  if (!file.exists(path)) skip("verify_r_outputs.py not available")
  normalizePath(path)
}

run_verifier <- function(...) {
  res <- system2(oracle_python(), c("-I", verifier(), ...), stdout = TRUE, stderr = TRUE)
  status <- attr(res, "status") %||% 0L
  list(ok = status == 0L && any(grepl("^OK:", res)), output = paste(res, collapse = "\n"))
}

test_that("Python mudata and scirpy read R-written H5MU and Zarr stores", {
  skip_if_no_hdf5()
  skip_if_no_zarr()
  py <- oracle_python()
  expected <- fixture("expected_airr.json")
  src <- read_h5mu(fixture("gex_airr.h5mu"), as = "HDF5MuData")
  outputs <- c(
    h5mu = withr::local_tempfile(fileext = ".h5mu"),
    v2 = withr::local_tempfile(fileext = ".zarr"),
    v3 = withr::local_tempfile(fileext = ".zarr")
  )
  write_h5mu(src, outputs[["h5mu"]])
  write_mudata_zarr(src, outputs[["v2"]], zarr_format = 2L)
  write_mudata_zarr(src, outputs[["v3"]], zarr_format = 3L)
  for (out in outputs) {
    res <- run_verifier("mudata", out, expected)
    expect_true(res$ok, info = res$output)
  }
})
