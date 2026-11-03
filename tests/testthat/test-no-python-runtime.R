test_that("native reads and writes never initialise Python", {
  skip_if_no_hdf5()
  skip_if_no_awkward()
  skip_on_cran()
  script <- withr::local_tempfile(fileext = ".R")
  out <- withr::local_tempfile(fileext = ".h5mu")
  writeLines(
    c(
      "Sys.setenv(RETICULATE_PYTHON = '/nonexistent/python')",
      sprintf(
        "if (requireNamespace('mudataR', quietly = TRUE)) library(mudataR) else pkgload::load_all('%s', quiet = TRUE)",
        normalizePath(test_path("..", ".."))
      ),
      sprintf("m <- read_h5mu('%s', as = 'HDF5MuData')", fixture("gex_airr.h5mu")),
      "x <- get_awkward(m[['airr']])$airr",
      "invisible(as.data.frame(x, unnest = TRUE))",
      sprintf("write_h5mu(read_h5mu('%s'), '%s')", fixture("gex_airr.h5mu"), out),
      "loaded <- 'reticulate' %in% loadedNamespaces()",
      "init <- loaded && reticulate::py_available(initialize = FALSE)",
      "cat(if (init) 'PYTHON_INITIALISED' else 'NO_PYTHON', '\\n')"
    ),
    script
  )
  res <- system2(file.path(R.home("bin"), "Rscript"), script, stdout = TRUE, stderr = TRUE)
  expect_true(any(grepl("NO_PYTHON", res)), info = paste(res, collapse = "\n"))
  expect_true(file.exists(out))
})
