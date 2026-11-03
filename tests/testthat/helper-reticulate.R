# ReticulateMuData tests need Python with mudata, anndata, awkward. They use
# the fixture environment named by MUDATAR_ORACLE_PYTHON and are skipped
# otherwise (native tests never need Python).
skip_if_no_reticulate_mudata <- function() {
  testthat::skip_if_not_installed("reticulate")
  py <- Sys.getenv("MUDATAR_ORACLE_PYTHON")
  if (!nzchar(py) || !file.exists(py)) {
    testthat::skip("MUDATAR_ORACLE_PYTHON is not set")
  }
  if (reticulate::py_available(initialize = FALSE)) {
    if (normalizePath(reticulate::py_config()$python) != normalizePath(py)) {
      testthat::skip("reticulate is already bound to a different Python")
    }
  } else {
    Sys.setenv(RETICULATE_PYTHON = py)
  }
  if (!reticulate::py_module_available("mudata")) {
    testthat::skip("Python mudata is not available")
  }
}

# anndataR warns once that Python anndata >= 0.13 is not fully supported
quietly <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("anndata >= 0.13|Awkward Arrays is currently experimental", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
