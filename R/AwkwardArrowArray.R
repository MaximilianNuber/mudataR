#' Awkward Arrays as Arrow data
#'
#' @description
#' `AwkwardArrowArray` is the R representation of an AnnData element encoded
#' as `awkward-array` (version `0.1.0`), such as the AIRR data that
#' [scirpy](https://scirpy.scverse.org) stores in `obsm[["airr"]]` and
#' `obsm[["chain_indices"]]`. It holds the Awkward *form*, *length* and typed
#' *buffers* exactly as stored on disk and converts them to nested Arrow data
#' on demand. No Python is involved.
#'
#' * **Lossless:** an array read from disk and not modified is written back
#'   with its original form and buffers.
#' * **Ragged-aware:** `x[i]` selects top-level elements (e.g. cells) without
#'   flattening, and the result keeps the Awkward node classes of the input.
#' * **Arrow:** [arrow::as_arrow_table()] returns an Arrow Table (one column
#'   per field for a top-level record, otherwise a single column named by
#'   `name`), [arrow::as_arrow_array()] an Arrow Array and
#'   [nanoarrow::as_nanoarrow_array()] the underlying nanoarrow array. Option
#'   types become nullable Arrow fields; non-option types are non-nullable, so
#'   `?string` and `string` stay distinct. Missing values, empty lists and
#'   empty strings are kept distinct.
#'
#' Supported Awkward node classes: `NumpyArray`, `EmptyArray`,
#' `ListOffsetArray`, `ListArray`, `RegularArray`, `RecordArray`,
#' `IndexedArray`, `IndexedOptionArray`, `ByteMaskedArray`, `BitMaskedArray`
#' and `UnmaskedArray`. Other forms (e.g. `UnionArray`, datetime primitives)
#' raise an `mudatar_unsupported_encoding` error naming the element and form
#' node.
#'
#' @param x An `AwkwardArrowArray`, or for `as_AwkwardArrowArray()` an Arrow
#'   or nanoarrow object
#' @param i Row selector: positive or negative positions or a logical vector
#' @param name Column name used by `arrow::as_arrow_table()` when the
#'   top-level type is not a record
#' @param schema,type Must be `NULL`
#' @param unnest For `as.data.frame()`, whether to return one row per element
#'   of a top-level list (e.g. one row per chain) instead of one row per cell
#' @param ... Unused
#'
#' @return `as_AwkwardArrowArray()` returns an `AwkwardArrowArray`;
#'   `arrow::as_arrow_table()` an Arrow Table; `length()` the top-level
#'   length.
#'
#' @examples
#' if (requireNamespace("nanoarrow", quietly = TRUE)) {
#'   x <- as_AwkwardArrowArray(
#'     nanoarrow::as_nanoarrow_array(data.frame(a = 1:3, b = c("x", "y", "z")))
#'   )
#'   x
#'   length(x[2:3])
#' }
#' @name AwkwardArrowArray
NULL

new_awkward_arrow_array <- function(
  form,
  form_json,
  length,
  buffers = NULL,
  loader = NULL,
  layout = NULL,
  dirty = FALSE,
  source = NULL
) {
  state <- new.env(parent = emptyenv())
  state$buffers <- buffers
  state$loader <- loader
  state$layout <- layout
  state$nanoarrow <- NULL
  structure(
    list(
      form = form,
      form_json = form_json,
      length = as.numeric(length),
      dirty = dirty,
      source = source,
      state = state
    ),
    class = "AwkwardArrowArray"
  )
}

awk_field <- function(x, field) {
  .subset2(x, field)
}

#' Buffers of an `AwkwardArrowArray`, loading lazily if needed
#' @noRd
awkward_buffers <- function(x) {
  state <- awk_field(x, "state")
  if (is.null(state$buffers) && !is.null(state$loader)) {
    state$buffers <- state$loader()
    state$loader <- NULL
  }
  state$buffers %||% list()
}

awkward_layout <- function(x) {
  state <- awk_field(x, "state")
  if (is.null(state$layout)) {
    state$layout <- awkward_decode(
      awk_field(x, "form"),
      awkward_buffers(x),
      awk_field(x, "length"),
      path = awk_field(x, "source")$path
    )
  }
  state$layout
}

#' Build an `AwkwardArrowArray` from stored form, length and buffers
#'
#' @param form_json The `form` attribute (JSON string)
#' @param length The `length` attribute
#' @param buffers Named list of buffers, or `NULL` when `loader` is given
#' @param loader Function returning the buffers (for lazy/backed reads)
#' @param source List describing where the element was read from
#' @noRd
awkward_from_store <- function(form_json, length, buffers = NULL, loader = NULL, source = NULL) {
  awkward_check_requires("Reading Awkward Arrays", c("jsonlite", "nanoarrow", "bit64"))
  form <- awkward_parse_form(form_json, path = source$path)
  specs <- awkward_buffer_specs(form)
  if (!is.null(buffers)) {
    missing <- setdiff(specs$key, names(buffers))
    if (length(missing)) {
      awkward_abort(
        "Missing Awkward buffers {.val {missing}}",
        path = source$path
      )
    }
  }
  new_awkward_arrow_array(
    form = form,
    form_json = form_json,
    length = length,
    buffers = buffers,
    loader = loader,
    source = source
  )
}

awkward_from_layout <- function(layout, source = NULL) {
  enc <- awkward_encode(layout)
  new_awkward_arrow_array(
    form = enc$form,
    form_json = enc$form_json,
    length = enc$length,
    buffers = enc$buffers,
    layout = layout,
    dirty = TRUE,
    source = source
  )
}

#' @rdname AwkwardArrowArray
#' @export
# nolint start: object_name_linter
as_AwkwardArrowArray <- function(x, ...) {
  # nolint end: object_name_linter
  UseMethod("as_AwkwardArrowArray")
}

#' @export
# nolint start: object_name_linter
as_AwkwardArrowArray.AwkwardArrowArray <- function(x, ...) {
  # nolint end: object_name_linter
  x
}

#' @export
# nolint start: object_name_linter
as_AwkwardArrowArray.default <- function(x, ...) {
  # nolint end: object_name_linter
  awkward_check_requires("Creating Awkward Arrays", c("jsonlite", "nanoarrow"))
  if (inherits(x, c("Table", "RecordBatch"))) {
    awkward_check_requires("Creating Awkward Arrays from Arrow objects", "arrow")
    # Columns are converted one by one: the field flags say whether each is
    # an option type (arrow's own RecordBatch conversion drops them)
    fields <- x$schema$fields
    layouts <- lapply(seq_along(fields), function(i) {
      arrow_column_layout(x[[i]], fields[[i]]$nullable, fields[[i]]$name)
    })
    if (length(fields) == 1L) {
      return(awkward_from_layout(layouts[[1]]))
    }
    return(awkward_from_layout(list(
      class = "RecordArray",
      fields = vapply(fields, function(f) f$name, character(1)),
      parameters = list(),
      length = x$num_rows,
      contents = layouts
    )))
  }
  if (inherits(x, c("ChunkedArray", "Array"))) {
    awkward_check_requires("Creating Awkward Arrays from Arrow objects", "arrow")
    return(awkward_from_layout(arrow_column_layout(x, TRUE, "value")))
  }
  array <- nanoarrow::as_nanoarrow_array(x)
  awkward_from_layout(awkward_layout_from_nanoarrow(array))
}

arrow_column_layout <- function(column, nullable, name) {
  if (inherits(column, "ChunkedArray")) {
    chunks <- column$chunks
    column <- if (length(chunks) == 1L) {
      chunks[[1]]
    } else if (length(chunks) == 0L) {
      arrow::concat_arrays(type = column$type)
    } else {
      # no `type`: casting nested lists while concatenating is unreliable
      do.call(arrow::concat_arrays, chunks)
    }
  }
  array <- nanoarrow::as_nanoarrow_array(column)
  schema <- nanoarrow::infer_nanoarrow_schema(array)
  if (!isTRUE(nullable)) {
    schema <- na_set_nullable(schema, FALSE)
  }
  awkward_layout_from_nanoarrow(array, schema, where = name)
}

#' @rdname AwkwardArrowArray
#' @export
# nolint start: object_name_linter
is_AwkwardArrowArray <- function(x) {
  # nolint end: object_name_linter
  inherits(x, "AwkwardArrowArray")
}

#' @rdname AwkwardArrowArray
#' @export
length.AwkwardArrowArray <- function(x) {
  awk_field(x, "length")
}

#' @export
dim.AwkwardArrowArray <- function(x) {
  NULL
}

#' @rdname AwkwardArrowArray
#' @export
`[.AwkwardArrowArray` <- function(x, i, ...) {
  if (...length() > 0) {
    cli_abort("Awkward arrays are indexed along their first axis only: use {.code x[i]}")
  }
  idx <- awkward_resolve_index(if (missing(i)) NULL else i, length(x))
  if (length(idx) == length(x) && all(idx == seq_along(idx) - 1)) {
    return(x)
  }
  awkward_from_layout(
    awkward_take(awkward_layout(x), idx),
    source = awk_field(x, "source")
  )
}

#' @export
`[[.AwkwardArrowArray` <- function(x, i, ...) {
  cli_abort("Use {.code x[i]} to select elements or {.fn as_arrow_table} to access fields")
}

#' @export
`$.AwkwardArrowArray` <- function(x, name) {
  cli_abort(
    "Use {.fn as_arrow_table} to access fields of an {.cls AwkwardArrowArray}"
  )
}

#' @export
format.AwkwardArrowArray <- function(x, ...) {
  type <- awkward_type_string(awk_field(x, "form"))
  paste0(
    "<AwkwardArrowArray> ",
    format(length(x), scientific = FALSE),
    " * ",
    type
  )
}

#' @export
print.AwkwardArrowArray <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
str.AwkwardArrowArray <- function(object, ...) {
  cat(" ", format(object), "\n", sep = "")
}

#' @rdname AwkwardArrowArray
#' @exportS3Method nanoarrow::as_nanoarrow_array
as_nanoarrow_array.AwkwardArrowArray <- function(x, ..., schema = NULL) {
  if (!is.null(schema)) {
    cli_abort("{.arg schema} is not supported for {.cls AwkwardArrowArray}")
  }
  state <- awk_field(x, "state")
  if (is.null(state$nanoarrow)) {
    state$nanoarrow <- awkward_layout_to_nanoarrow(awkward_layout(x))
  }
  state$nanoarrow
}

#' @rdname AwkwardArrowArray
#' @exportS3Method arrow::as_arrow_array
as_arrow_array.AwkwardArrowArray <- function(x, ..., type = NULL) {
  awkward_check_requires("Converting Awkward Arrays to Arrow", "arrow")
  if (!is.null(type)) {
    cli_abort("{.arg type} is not supported for {.cls AwkwardArrowArray}")
  }
  arrow::as_arrow_array(as_nanoarrow_array.AwkwardArrowArray(x))
}

#' @rdname AwkwardArrowArray
#' @exportS3Method arrow::as_arrow_table
as_arrow_table.AwkwardArrowArray <- function(x, ..., schema = NULL, name = "value") {
  awkward_check_requires("Converting Awkward Arrays to Arrow", "arrow")
  if (!is.null(schema)) {
    cli_abort("{.arg schema} is not supported for {.cls AwkwardArrowArray}")
  }
  array <- as_nanoarrow_array.AwkwardArrowArray(x)
  schema <- nanoarrow::infer_nanoarrow_schema(array)
  if (identical(schema$format, "+s") && !na_flag_nullable(schema)) {
    return(arrow::as_arrow_table(arrow::as_record_batch(array)))
  }
  column <- arrow::as_arrow_array(array)
  columns <- stats::setNames(list(column), name)
  arrow::Table$create(
    !!!columns,
    schema = arrow::schema(arrow::field(
      name,
      column$type,
      nullable = na_flag_nullable(schema)
    ))
  )
}

#' @rdname AwkwardArrowArray
#' @export
as.data.frame.AwkwardArrowArray <- function(x, ..., unnest = FALSE) {
  layout <- awkward_layout(x)
  if (!isTRUE(unnest)) {
    out <- nanoarrow::convert_array(as_nanoarrow_array.AwkwardArrowArray(x))
    if (!is.data.frame(out)) {
      out <- data.frame(value = I(out))
    }
    return(out)
  }
  if (!layout$class %in% c("ListOffsetArray", "ListArray") || !is.null(string_kind(layout))) {
    cli_abort("{.arg unnest} requires a top-level list (e.g. cells with chains)")
  }
  layout <- awkward_take(layout, seq_len(layout$length) - 1)
  lens <- diff(layout$offsets_data)
  total <- sum(lens)
  if (total > 1e7) {
    cli_warn("Unnesting produces {format(total, big.mark = ',')} rows")
  }
  inner <- awkward_from_layout(layout[["content"]])
  values <- as.data.frame.AwkwardArrowArray(inner)
  cbind(.row = rep(seq_len(layout$length), lens), values)
}

#' @rdname AwkwardArrowArray
#' @export
# nolint start: object_name_linter
validate_AwkwardArrowArray <- function(x) {
  # nolint end: object_name_linter
  awkward_layout(x)
  awkward_layout_to_nanoarrow(awkward_layout(x))
  invisible(x)
}
