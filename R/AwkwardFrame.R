#' Awkward Arrays inside anndataR objects
#'
#' @description
#' anndataR does not know Awkward Arrays, but it subsets every `data.frame` in
#' `obsm`/`varm` with its standard rules (views, `adata[i, ]`,
#' `as_InMemoryAnnData()`). An `AwkwardFrame` uses exactly that: it is a
#' `data.frame` with one row per observation (or variable) and a single
#' column, `.awkward_row`, holding the position of that row in the
#' [AwkwardArrowArray] it carries. anndataR subsets the row positions; the
#' Awkward data are only touched when you ask for them.
#'
#' * `as_AwkwardArrowArray(x)` returns the selected rows as an
#'   [AwkwardArrowArray]. When no rows were dropped or reordered it is the
#'   original array, so unchanged data are written back with their original
#'   form and buffers.
#' * [arrow::as_arrow_table()] and [nanoarrow::as_nanoarrow_array()] work on
#'   an `AwkwardFrame` directly.
#' * mudataR writes `AwkwardFrame`s as `awkward-array` elements; anndataR's own
#'   writers never see them.
#'
#' In-memory modalities read by mudataR hold scirpy data such as
#' `obsm[["airr"]]` as `AwkwardFrame`s. For backed modalities use
#' [get_awkward()].
#'
#' @param x An [AwkwardArrowArray], an Arrow or nanoarrow object (for
#'   `as_AwkwardFrame()`), or an `AwkwardFrame`
#' @param i,j,drop Row selection, as for data frames; `j` must be missing
#' @param ... Passed on
#'
#' @return `as_AwkwardFrame()` an `AwkwardFrame`; `is_AwkwardFrame()` a flag
#' @examples
#' if (requireNamespace("nanoarrow", quietly = TRUE)) {
#'   awk <- as_AwkwardArrowArray(
#'     nanoarrow::as_nanoarrow_array(data.frame(a = 1:3, b = c("x", "y", "z")))
#'   )
#'   adata <- anndataR::AnnData(
#'     obs = data.frame(row.names = c("c1", "c2", "c3")),
#'     var = data.frame(row.names = "g1"),
#'     obsm = list(chains = as_AwkwardFrame(awk))
#'   )
#'   # anndataR subsets the Awkward rows with its standard rules
#'   as_AwkwardArrowArray(adata[c(3, 1), ]$obsm$chains)
#' }
#' @name AwkwardFrame
NULL

new_awkward_frame <- function(awkward, rows) {
  x <- data.frame(.awkward_row = as.integer(rows))
  # set attributes one by one: structure() would expand the compact
  # automatic row names, which anndataR then treats as real row names
  attr(x, "awkward") <- awkward
  class(x) <- c("AwkwardFrame", "data.frame")
  x
}

#' @rdname AwkwardFrame
#' @export
# nolint start: object_name_linter
as_AwkwardFrame <- function(x) {
  # nolint end: object_name_linter
  if (is_AwkwardFrame(x)) {
    return(x)
  }
  if (!inherits(x, "AwkwardArrowArray")) {
    x <- as_AwkwardArrowArray(x)
  }
  n <- length(x)
  if (n > .Machine$integer.max) {
    mudatar_abort_unsupported("Awkward Arrays with more than 2^31 - 1 rows inside anndataR objects")
  }
  new_awkward_frame(x, seq_len(n))
}

#' @rdname AwkwardFrame
#' @export
# nolint start: object_name_linter
is_AwkwardFrame <- function(x) {
  # nolint end: object_name_linter
  inherits(x, "AwkwardFrame")
}

awkward_frame_rows <- function(x) {
  rows <- .subset2(x, ".awkward_row")
  names(rows) <- NULL
  rows
}

#' @rdname AwkwardFrame
#' @export
`[.AwkwardFrame` <- function(x, i, j, drop = FALSE) {
  if (!missing(j)) {
    mudatar_abort("An {.cls AwkwardFrame} is indexed along its rows only: use {.code x[i, ]}")
  }
  rows <- awkward_frame_rows(x)
  if (!missing(i)) {
    if (is.character(i)) {
      i <- match(i, rownames(x))
    }
    rows <- rows[i]
  }
  new_awkward_frame(attr(x, "awkward"), rows)
}

#' @rdname AwkwardFrame
#' @export
# nolint start: object_name_linter
as_AwkwardArrowArray.AwkwardFrame <- function(x, ...) {
  # nolint end: object_name_linter
  awk <- attr(x, "awkward")
  rows <- awkward_frame_rows(x)
  if (anyNA(rows)) {
    mudatar_abort_unsupported("Missing (padding) rows in an Awkward Array")
  }
  if (length(rows) == length(awk) && identical(rows, seq_along(rows))) {
    return(awk)
  }
  awk[rows]
}

#' @export
as.matrix.AwkwardFrame <- function(x, ...) {
  # Not numeric, so anndataR never mistakes it for an embedding when
  # converting to a SingleCellExperiment
  matrix(rep("<awkward>", nrow(x)), ncol = 1L)
}

#' @export
format.AwkwardFrame <- function(x, ...) {
  awk <- attr(x, "awkward")
  paste0(
    "<AwkwardFrame> ",
    nrow(x),
    " * ",
    awkward_type_string(awk_field(awk, "form")),
    " (rows of an AwkwardArrowArray; see ?AwkwardFrame)"
  )
}

#' @export
print.AwkwardFrame <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
str.AwkwardFrame <- function(object, ...) {
  cat(" ", format(object), "\n", sep = "")
}

#' @exportS3Method nanoarrow::as_nanoarrow_array
as_nanoarrow_array.AwkwardFrame <- function(x, ..., schema = NULL) {
  as_nanoarrow_array.AwkwardArrowArray(as_AwkwardArrowArray(x), ..., schema = schema)
}

#' @exportS3Method arrow::as_arrow_table
as_arrow_table.AwkwardFrame <- function(x, ..., schema = NULL, name = "value") {
  as_arrow_table.AwkwardArrowArray(as_AwkwardArrowArray(x), ..., schema = schema, name = name)
}

#' @exportS3Method arrow::as_arrow_array
as_arrow_array.AwkwardFrame <- function(x, ..., type = NULL) {
  as_arrow_array.AwkwardArrowArray(as_AwkwardArrowArray(x), ..., type = type)
}

#' Is this an Awkward value (an array, or rows of one inside anndataR)?
#' @noRd
is_awkward_value <- function(x) {
  inherits(x, c("AwkwardArrowArray", "AwkwardFrame"))
}
