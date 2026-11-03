#' Standard R methods for MuData objects
#'
#' `dim()`, `nrow()`, `ncol()`, `dimnames()`, `rownames()` and `colnames()`
#' describe the global observations and variables. `mdata[["name"]]`
#' extracts a modality (the backend-matching anndataR AnnData object) and
#' `mdata[["name"]] <- adata` adds or replaces one and rebuilds the global
#' indices (assign `NULL` to remove it).
#'
#' @param x A MuData object
#' @param i Modality name or position
#' @param j Unused (subsetting is not supported yet)
#' @param j Unused (subsetting is not supported yet)
#' @param value An anndataR AnnData object or `NULL`
#' @param ... Unused
#'
#' @name MuData-s3methods
NULL

#' @rdname MuData-s3methods
#' @export
dim.AbstractMuData <- function(x) {
  x$shape()
}

#' @rdname MuData-s3methods
#' @export
dimnames.AbstractMuData <- function(x) {
  list(x$obs_names, x$var_names)
}

#' @rdname MuData-s3methods
#' @export
`[[.AbstractMuData` <- function(x, i, ...) {
  x$get_mod(resolve_mod_name(x$mod_names(), i))
}

#' @rdname MuData-s3methods
#' @export
`[[<-.AbstractMuData` <- function(x, i, value) {
  if (is.numeric(i)) {
    i <- resolve_mod_name(x$mod_names(), i)
  }
  x$set_mod(i, value)
  x
}

#' @rdname MuData-s3methods
#' @export
`[.AbstractMuData` <- function(x, i, j, ...) {
  mudatar_abort_unsupported(
    "Subsetting MuData objects with {.code [}",
    "i" = "Lazy MuData views are planned; see {.file docs/parity.md}."
  )
}

#' Test for MuData objects
#'
#' @param x Any object
#' @return `TRUE` or `FALSE`
#' @name is_mudata
#' @export
is_mudata <- function(x) {
  inherits(x, "AbstractMuData")
}

#' @rdname is_mudata
#' @export
is_InMemoryMuData <- function(x) {
  # nolint start: object_name_linter
  inherits(x, "InMemoryMuData")
  # nolint end: object_name_linter
}

#' @rdname is_mudata
#' @export
is_HDF5MuData <- function(x) {
  # nolint start: object_name_linter
  inherits(x, "HDF5MuData")
  # nolint end: object_name_linter
}

#' @rdname is_mudata
#' @export
is_ZarrMuData <- function(x) {
  # nolint start: object_name_linter
  inherits(x, "ZarrMuData")
  # nolint end: object_name_linter
}
