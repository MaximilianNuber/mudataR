# Awkward Array forms, buffers and layouts
#
# AnnData's `awkward-array` encoding (0.1.0) stores an Awkward Array as a
# group with
#   * attribute `form`:   JSON describing the layout (awkward v2 form),
#   * attribute `length`: the top-level length,
#   * one dataset per buffer, named `<form_key>-<role>` (e.g. `node0-offsets`).
#
# This file implements the format-level pieces, independent of any store:
#   * `awkward_parse_form()`     JSON -> validated, normalised form AST
#   * `awkward_buffer_specs()`   form -> required buffers and their dtypes
#   * `awkward_decode()`         form + buffers + length -> layout tree
#   * `awkward_take()`           vectorised row selection on a layout
#   * `awkward_encode()`         layout -> form + buffers (canonical keys)
#
# A *layout* is the form AST with data attached (`$data`, `$offsets`,
# `$index`, `$mask`, ...) and a `$length`. All operations work on whole
# buffers; there are no per-element R loops.
#
# Index-like buffers (offsets, starts, stops, index) are held as R doubles,
# which represent every integer up to 2^53 exactly. NumpyArray data keep a
# type that is exact for their primitive (see `awkward_normalize_data()`).

# Supported node classes. Anything else (notably UnionArray) is rejected with
# a classed condition naming the node.
.awkward_supported_classes <- c(
  "NumpyArray",
  "EmptyArray",
  "ListOffsetArray",
  "ListArray",
  "RegularArray",
  "RecordArray",
  "IndexedArray",
  "IndexedOptionArray",
  "ByteMaskedArray",
  "BitMaskedArray",
  "UnmaskedArray"
)

.awkward_index_types <- c(
  i8 = "int8",
  u8 = "uint8",
  i32 = "int32",
  u32 = "uint32",
  i64 = "int64"
)

.awkward_primitives <- c(
  "bool",
  "int8",
  "uint8",
  "int16",
  "uint16",
  "int32",
  "uint32",
  "int64",
  "uint64",
  "float32",
  "float64"
)

#' Signal an Awkward Array error
#'
#' @param message Message (cli format)
#' @param path Element path in the store, if known
#' @param node Path of the offending node inside the form, if known
#' @param class Additional condition class
#' @noRd
awkward_abort <- function(
  message,
  ...,
  path = NULL,
  node = NULL,
  class = NULL,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  bullets <- c(message, ...)
  if (!is.null(path)) {
    bullets <- c(bullets, "i" = paste0("Element: {.path ", path, "}"))
  }
  if (!is.null(node)) {
    bullets <- c(bullets, "i" = paste0("Form node: {.code ", node, "}"))
  }
  cli_abort(
    bullets,
    class = c(class, "mudatar_awkward_error", "mudatar_error"),
    element_path = path,
    form_node = node,
    call = call,
    .envir = .envir
  )
}

#' Signal an unsupported Awkward form
#'
#' The condition carries `store_kind`, `element_path`, `encoding_type`,
#' `encoding_version` and `form_node` fields so callers can report or handle it
#' programmatically. It is never downgraded to a warning by the element
#' readers.
#' @noRd
awkward_abort_unsupported <- function(
  what,
  node,
  path = NULL,
  call = rlang::caller_env()
) {
  supported <- .awkward_supported_classes # nolint: object_usage_linter. used in cli glue
  cli_abort(
    c(
      paste("Unsupported Awkward Array form:", what),
      "i" = "Form node: {.code {node}}",
      if (!is.null(path)) c("i" = "Element: {.path {path}}"),
      "i" = "mudataR supports the node classes {.val {supported}}.",
      ">" = "Please report the file at {.url https://github.com/MaximilianNuber/mudataR/issues} if you need this form."
    ),
    class = c(
      "mudatar_unsupported_encoding",
      "mudatar_awkward_error"
    ),
    encoding_type = "awkward-array",
    encoding_version = "0.1.0",
    element_path = path,
    form_node = node,
    call = call
  )
}

# Form parsing -----------------------------------------------------------------

#' Parse an Awkward form
#'
#' @param form A JSON string (as stored in the `form` attribute) or an already
#'   parsed list
#' @param path Element path, for error messages
#'
#' @return A normalised form AST: nested lists with a `class` field
#' @noRd
awkward_parse_form <- function(form, path = NULL) {
  if (is.character(form)) {
    if (length(form) != 1L) {
      awkward_abort("The {.field form} attribute must be a single string", path = path)
    }
    form <- tryCatch(
      jsonlite::fromJSON(form, simplifyVector = FALSE),
      error = function(e) {
        awkward_abort("The {.field form} attribute is not valid JSON", path = path)
      }
    )
  }
  awkward_normalize_form(form, "form", path)
}

awkward_normalize_form <- function(node, where, path) {
  # awkward permits a bare primitive name as shorthand for a NumpyArray
  if (is.character(node) && length(node) == 1L) {
    node <- list(class = "NumpyArray", primitive = node)
  }
  if (!is.list(node) || is.null(node[["class"]])) {
    awkward_abort("Malformed form node", path = path, node = where)
  }

  cls <- node[["class"]]
  # awkward v1 class names carry the index width as a suffix
  # (e.g. ListOffsetArray64); v2 forms store it in a separate field
  v1 <- regmatches(cls, regexec("^(.*Array)(64|32|U32|8|U8)$", cls))[[1]]
  if (length(v1) == 3L) {
    cls <- v1[[2]]
    width <- c("64" = "i64", "32" = "i32", "U32" = "u32", "8" = "i8", "U8" = "u8")[[v1[[3]]]]
    for (field in c("offsets", "starts", "index")) {
      node[[field]] <- node[[field]] %||% width
    }
  }

  if (!cls %in% .awkward_supported_classes) {
    awkward_abort_unsupported(paste0("node class {.cls ", cls, "}"), where, path)
  }

  params <- node[["parameters"]] %||% list()
  out <- list(
    class = cls,
    parameters = params,
    form_key = node[["form_key"]]
  )

  check_index <- function(field) {
    value <- node[[field]] %||% "i64"
    if (!value %in% names(.awkward_index_types)) {
      awkward_abort_unsupported(
        paste0("index type {.val ", value, "} for {.field ", field, "}"),
        where,
        path
      )
    }
    value
  }
  sub_node <- function(child, label) {
    awkward_normalize_form(child, paste0(where, ".", label), path)
  }

  switch(
    cls,
    NumpyArray = {
      primitive <- node[["primitive"]]
      if (is.null(primitive) && !is.null(node[["format"]])) {
        awkward_abort_unsupported("NumpyArray without a {.field primitive}", where, path)
      }
      if (!primitive %in% .awkward_primitives) {
        awkward_abort_unsupported(
          paste0("NumpyArray primitive {.val ", primitive, "}"),
          where,
          path
        )
      }
      out$primitive <- primitive
      inner_shape <- unlist(node[["inner_shape"]] %||% list())
      if (length(inner_shape) > 0) {
        # A NumpyArray with inner dimensions is equivalent to nested
        # RegularArrays around a flat NumpyArray holding the same buffer
        out$inner_shape <- list()
        leaf <- out
        for (size in rev(inner_shape)) {
          leaf <- list(
            class = "RegularArray",
            parameters = list(),
            form_key = NULL,
            size = as.numeric(size),
            content = leaf
          )
        }
        leaf$parameters <- params
        return(leaf)
      }
      out$inner_shape <- list()
    },
    EmptyArray = NULL,
    ListOffsetArray = {
      out$offsets <- check_index("offsets")
      out[["content"]] <- sub_node(node[["content"]], "content")
    },
    ListArray = {
      out$starts <- check_index("starts")
      out$stops <- node[["stops"]] %||% out$starts
      out[["content"]] <- sub_node(node[["content"]], "content")
    },
    RegularArray = {
      out$size <- as.numeric(node[["size"]])
      out[["content"]] <- sub_node(node[["content"]], "content")
    },
    RecordArray = {
      contents <- node[["contents"]]
      fields <- node[["fields"]]
      if (is.list(contents) && !is.null(names(contents)) && is.null(fields)) {
        # awkward v1 records store contents as an object keyed by field
        fields <- names(contents)
        contents <- unname(contents)
      }
      out$fields <- if (is.null(fields)) NULL else as.character(unlist(fields))
      out$contents <- lapply(seq_along(contents), function(i) {
        sub_node(contents[[i]], paste0("contents[", i - 1L, "]"))
      })
      if (!is.null(out$fields) && length(out$fields) != length(out$contents)) {
        awkward_abort("RecordArray has mismatched fields and contents", path = path, node = where)
      }
    },
    IndexedArray = ,
    IndexedOptionArray = {
      out$index <- check_index("index")
      out[["content"]] <- sub_node(node[["content"]], "content")
    },
    ByteMaskedArray = {
      out$mask <- node[["mask"]] %||% "i8"
      out$valid_when <- isTRUE(node[["valid_when"]])
      out[["content"]] <- sub_node(node[["content"]], "content")
    },
    BitMaskedArray = {
      out$mask <- node[["mask"]] %||% "u8"
      out$valid_when <- isTRUE(node[["valid_when"]])
      out$lsb_order <- isTRUE(node[["lsb_order"]])
      out[["content"]] <- sub_node(node[["content"]], "content")
    },
    UnmaskedArray = {
      out[["content"]] <- sub_node(node[["content"]], "content")
    }
  )
  out
}

#' Buffers required by a form
#'
#' @return A data frame with columns `key` (dataset name), `role`, `dtype`
#'   (numpy-style type name) and `node` (form node class)
#' @noRd
awkward_buffer_specs <- function(form) {
  specs <- list()
  add <- function(node, role, dtype) {
    if (is.null(node$form_key)) {
      return()
    }
    specs[[length(specs) + 1L]] <<- data.frame(
      key = paste0(node$form_key, "-", role),
      role = role,
      dtype = dtype,
      node = node$class
    )
  }
  walk <- function(node) {
    switch(
      node$class,
      NumpyArray = add(node, "data", node$primitive),
      ListOffsetArray = add(node, "offsets", .awkward_index_types[[node$offsets]]),
      ListArray = {
        add(node, "starts", .awkward_index_types[[node$starts]])
        add(node, "stops", .awkward_index_types[[node$stops]])
      },
      IndexedArray = ,
      IndexedOptionArray = add(node, "index", .awkward_index_types[[node$index]]),
      ByteMaskedArray = add(node, "mask", "int8"),
      BitMaskedArray = add(node, "mask", "uint8")
    )
    if (!is.null(node[["content"]])) {
      walk(node[["content"]])
    }
    for (child in node$contents) {
      walk(child)
    }
  }
  walk(form)
  if (length(specs) == 0) {
    return(data.frame(key = character(), role = character(), dtype = character(), node = character()))
  }
  do.call(rbind, specs)
}

# Data normalisation ----------------------------------------------------------------

#' Normalise a buffer read from a store to its in-memory representation
#'
#' Stores return different R types for the same dtype (rhdf5 returns
#' `integer64` for int64 and a factor for HDF5 booleans; Rarr returns integers
#' for int64). The form is the authority on the dtype, so buffers are coerced
#' here:
#'   * index dtypes (`int8`..`int64` as index) -> double
#'   * bool -> logical
#'   * (u)int8/16, int32 -> integer
#'   * uint32, float32, float64 -> double
#'   * int64/uint64 -> double when every value is exactly representable,
#'     otherwise `bit64::integer64` (never silently rounded)
#'   * uint8 carrying characters/bytes -> raw
#'
#' @noRd
awkward_normalize_data <- function(x, dtype, raw_bytes = FALSE, key = NULL, path = NULL) {
  if (is.factor(x)) {
    # HDF5 enum booleans (h5py convention)
    x <- as.character(x) == "TRUE"
  }
  x <- awkward_flatten(x)
  if (raw_bytes) {
    return(if (is.raw(x)) x else as.raw(as.integer(x)))
  }
  if (is.raw(x)) {
    x <- as.integer(x)
  }
  switch(
    dtype,
    bool = as.logical(x),
    int8 = ,
    uint8 = ,
    int16 = ,
    uint16 = ,
    int32 = as.integer(x),
    uint32 = ,
    float32 = ,
    float64 = as.double(x),
    int64 = ,
    uint64 = {
      if (inherits(x, "integer64")) {
        big <- !is.na(x) & abs(x) > bit64::as.integer64(2^53)
        if (any(big)) x else as.double(x)
      } else {
        as.double(x)
      }
    },
    awkward_abort(
      "Unsupported buffer dtype {.val {dtype}}",
      path = path,
      node = key
    )
  )
}

#' Drop dimensions without losing the `integer64` class
#'
#' `as.vector()` strips the class of `bit64::integer64` vectors and would
#' reinterpret their bits as doubles.
#' @noRd
awkward_flatten <- function(x) {
  if (inherits(x, "integer64")) {
    dim(x) <- NULL
    return(x)
  }
  as.vector(x)
}

as_awkward_index <- function(x, key, path) {
  x <- awkward_flatten(x)
  if (inherits(x, "integer64")) {
    if (any(!is.na(x) & abs(x) > bit64::as.integer64(2^53))) {
      awkward_abort(
        "Index buffer {.val {key}} has values beyond 2^53, which cannot be indexed in R",
        path = path
      )
    }
    x <- as.double(x)
  }
  x <- as.double(if (is.raw(x)) as.integer(x) else x)
  if (anyNA(x)) {
    awkward_abort("Index buffer {.val {key}} contains missing values", path = path)
  }
  x
}

is_string_like <- function(node) {
  identical(node$class, "ListOffsetArray") || identical(node$class, "ListArray")
}

string_kind <- function(node) {
  kind <- node$parameters[["__array__"]]
  if (is.null(kind) || !kind %in% c("string", "bytestring")) {
    return(NULL)
  }
  kind
}

# Decoding -------------------------------------------------------------------

#' Decode buffers into a layout
#'
#' @param form A parsed form (see `awkward_parse_form()`)
#' @param buffers A named list of buffers keyed by `<form_key>-<role>`
#' @param length The top-level length
#' @param path Element path, for error messages
#'
#' @return A layout: the form AST with buffers attached and `$length` set on
#'   every node
#' @noRd
awkward_decode <- function(form, buffers, length, path = NULL) {
  length <- as.numeric(length)
  if (length(length) != 1L || is.na(length) || length < 0) {
    awkward_abort("Invalid Awkward Array {.field length}", path = path)
  }
  awkward_decode_node(form, length, buffers, "form", path, raw_bytes = FALSE)
}

awkward_decode_node <- function(node, length, buffers, where, path, raw_bytes) {
  get_buffer <- function(role, n) {
    key <- paste0(node$form_key, "-", role)
    if (is.null(node$form_key) || is.null(buffers[[key]])) {
      awkward_abort(
        "Missing buffer {.val {key}}",
        path = path,
        node = where
      )
    }
    buf <- buffers[[key]]
    if (base::length(buf) < n) {
      awkward_abort(
        "Buffer {.val {key}} has {base::length(buf)} elements but at least {n} are required",
        path = path,
        node = where
      )
    }
    if (base::length(buf) > n) {
      buf <- buf[seq_len(n)]
    }
    list(key = key, value = buf)
  }

  node$length <- length

  switch(
    node$class,
    NumpyArray = {
      b <- get_buffer("data", length)
      node$data <- awkward_normalize_data(
        b$value,
        node$primitive,
        raw_bytes = raw_bytes && node$primitive == "uint8",
        key = b$key,
        path = path
      )
    },
    EmptyArray = {
      if (length != 0) {
        awkward_abort("EmptyArray must have length 0, got {length}", path = path, node = where)
      }
    },
    ListOffsetArray = {
      b <- get_buffer("offsets", length + 1)
      offsets <- as_awkward_index(b$value, b$key, path)
      if (any(offsets < 0) || is.unsorted(offsets)) {
        awkward_abort(
          "Offsets {.val {b$key}} must be non-negative and non-decreasing",
          path = path,
          node = where
        )
      }
      node$offsets_data <- offsets
      content_length <- if (length == 0) 0 else offsets[[length + 1]]
      node[["content"]] <- awkward_decode_node(
        node[["content"]],
        content_length,
        buffers,
        paste0(where, ".content"),
        path,
        raw_bytes = !is.null(string_kind(node))
      )
    },
    ListArray = {
      starts <- get_buffer("starts", length)
      stops <- get_buffer("stops", length)
      node$starts_data <- as_awkward_index(starts$value, starts$key, path)
      node$stops_data <- as_awkward_index(stops$value, stops$key, path)
      if (any(node$stops_data < node$starts_data) || any(node$starts_data < 0)) {
        awkward_abort("ListArray stops must not precede starts", path = path, node = where)
      }
      nonempty <- node$stops_data > node$starts_data
      content_length <- if (any(nonempty)) max(node$stops_data[nonempty]) else 0
      node[["content"]] <- awkward_decode_node(
        node[["content"]],
        content_length,
        buffers,
        paste0(where, ".content"),
        path,
        raw_bytes = !is.null(string_kind(node))
      )
    },
    RegularArray = {
      if (is.na(node$size) || node$size < 0) {
        awkward_abort("RegularArray size must be non-negative", path = path, node = where)
      }
      node[["content"]] <- awkward_decode_node(
        node[["content"]],
        length * node$size,
        buffers,
        paste0(where, ".content"),
        path,
        raw_bytes = FALSE
      )
    },
    RecordArray = {
      node$contents <- lapply(seq_along(node$contents), function(i) {
        awkward_decode_node(
          node$contents[[i]],
          length,
          buffers,
          paste0(where, ".contents[", i - 1L, "]"),
          path,
          raw_bytes = FALSE
        )
      })
    },
    IndexedArray = ,
    IndexedOptionArray = {
      b <- get_buffer("index", length)
      index <- as_awkward_index(b$value, b$key, path)
      option <- node$class == "IndexedOptionArray"
      if (!option && any(index < 0)) {
        awkward_abort("IndexedArray index must be non-negative", path = path, node = where)
      }
      valid <- index >= 0
      content_length <- if (any(valid)) max(index[valid]) + 1 else 0
      node$index_data <- index
      node[["content"]] <- awkward_decode_node(
        node[["content"]],
        content_length,
        buffers,
        paste0(where, ".content"),
        path,
        raw_bytes = raw_bytes
      )
    },
    ByteMaskedArray = {
      b <- get_buffer("mask", length)
      node$mask_data <- as.integer(as.vector(if (is.raw(b$value)) as.integer(b$value) else b$value)) != 0L
      node[["content"]] <- awkward_decode_node(
        node[["content"]],
        length,
        buffers,
        paste0(where, ".content"),
        path,
        raw_bytes = raw_bytes
      )
    },
    BitMaskedArray = {
      b <- get_buffer("mask", ceiling(length / 8))
      bytes <- awkward_normalize_data(b$value, "uint8", raw_bytes = TRUE)
      bits <- matrix(as.logical(rawToBits(bytes)), nrow = 8L)
      if (!node$lsb_order) {
        bits <- bits[8:1, , drop = FALSE]
      }
      # Normalise to a ByteMaskedArray: same semantics, simpler to operate on
      node$class <- "ByteMaskedArray"
      node$mask <- "i8"
      node$lsb_order <- NULL
      node$mask_data <- as.vector(bits)[seq_len(length)]
      node[["content"]] <- awkward_decode_node(
        node[["content"]],
        length,
        buffers,
        paste0(where, ".content"),
        path,
        raw_bytes = raw_bytes
      )
    },
    UnmaskedArray = {
      node[["content"]] <- awkward_decode_node(
        node[["content"]],
        length,
        buffers,
        paste0(where, ".content"),
        path,
        raw_bytes = raw_bytes
      )
    }
  )
  node
}

# Take -------------------------------------------------------------------------

#' Select elements of a layout
#'
#' @param node A layout
#' @param idx Zero-based double indices. `NA` produces a placeholder element
#'   (zero-length list, zero value, null option); placeholders are only valid
#'   underneath a null parent and are used internally by the Arrow converter.
#'
#' @return A layout of length `length(idx)`. List nodes come back as
#'   `ListOffsetArray` with offsets starting at zero; option nodes compose
#'   their index instead of materialising content.
#' @noRd
awkward_take <- function(node, idx) {
  n <- length(idx)
  if (n == node$length && !anyNA(idx) && all(idx == seq_len(n) - 1)) {
    return(node)
  }
  out <- node
  out$length <- n
  switch(
    node$class,
    NumpyArray = {
      out$data <- node$data[idx + 1]
    },
    NullArray = NULL,
    EmptyArray = {
      if (n > 0) {
        if (!all(is.na(idx))) {
          cli_abort("Cannot take non-placeholder elements from an EmptyArray")
        }
        out <- list(class = "NullArray", parameters = list(), length = n)
      }
    },
    ListOffsetArray = ,
    ListArray = {
      if (node$class == "ListOffsetArray") {
        starts <- node$offsets_data[idx + 1]
        stops <- node$offsets_data[idx + 2]
      } else {
        starts <- node$starts_data[idx + 1]
        stops <- node$stops_data[idx + 1]
      }
      starts[is.na(starts)] <- 0
      lens <- stops - starts
      lens[is.na(lens)] <- 0
      out$class <- "ListOffsetArray"
      out$offsets <- if (identical(node$class, "ListOffsetArray")) node$offsets else node$starts
      out$starts <- NULL
      out$stops <- NULL
      out$starts_data <- NULL
      out$stops_data <- NULL
      out$offsets_data <- c(0, cumsum(lens))
      inner <- if (sum(lens) > 0) {
        check_int_range(starts + lens)
        sequence(as.integer(lens), from = as.integer(starts))
      } else {
        numeric(0)
      }
      out[["content"]] <- awkward_take(node[["content"]], as.double(inner))
    },
    RegularArray = {
      size <- node$size
      inner <- if (size > 0 && n > 0) {
        rep(idx * size, each = size) + rep(seq_len(size) - 1, times = n)
      } else {
        numeric(0)
      }
      out[["content"]] <- awkward_take(node[["content"]], inner)
    },
    RecordArray = {
      out$contents <- lapply(node$contents, awkward_take, idx = idx)
    },
    IndexedArray = {
      out$index_data <- node$index_data[idx + 1]
    },
    IndexedOptionArray = {
      index <- node$index_data[idx + 1]
      index[is.na(index)] <- -1
      out$index_data <- index
    },
    ByteMaskedArray = {
      mask <- node$mask_data[idx + 1]
      mask[is.na(mask)] <- !node$valid_when
      out$mask_data <- mask
      out[["content"]] <- awkward_take(node[["content"]], idx)
    },
    UnmaskedArray = {
      out[["content"]] <- awkward_take(node[["content"]], idx)
    }
  )
  out
}

check_int_range <- function(x) {
  if (length(x) && max(x) > .Machine$integer.max) {
    cli_abort(
      "Awkward Array content with more than {.Machine$integer.max} elements is not supported yet"
    )
  }
}

#' Resolve a row selector to zero-based indices
#' @noRd
awkward_resolve_index <- function(i, n) {
  if (missing(i) || is.null(i)) {
    return(seq_len(n) - 1)
  }
  if (anyNA(i)) {
    cli_abort("Index must not contain {.val NA}")
  }
  if (is.logical(i)) {
    if (length(i) != n) {
      cli_abort("Logical index must have length {n}, not {length(i)}")
    }
    if (anyNA(i)) {
      cli_abort("Logical index must not contain {.val NA}")
    }
    return(which(i) - 1)
  }
  if (!is.numeric(i)) {
    cli_abort("Awkward arrays can only be indexed by position or logical vectors")
  }
  if (anyNA(i)) {
    cli_abort("Index must not contain {.val NA}")
  }
  i <- as.double(i)
  if (any(i < 0)) {
    if (any(i > 0)) {
      cli_abort("Cannot mix positive and negative indices")
    }
    keep <- rep(TRUE, n)
    keep[-i] <- FALSE
    return(which(keep) - 1)
  }
  i <- i[i != 0]
  if (any(i > n)) {
    cli_abort("Index out of bounds: the array has length {n}")
  }
  i - 1
}

# Encoding ---------------------------------------------------------------------

#' Encode a layout into a form and buffers
#'
#' Form keys are reassigned as `node<k>` in pre-order, as awkward does.
#'
#' @return A list with `form` (AST), `form_json` (string), `length` and
#'   `buffers` (named list)
#' @noRd
awkward_encode <- function(layout) {
  counter <- 0L
  buffers <- list()
  emit <- function(key, value) {
    buffers[[key]] <<- value
  }
  enc <- function(node) {
    key <- paste0("node", counter)
    counter <<- counter + 1L
    form <- list(class = node$class)
    switch(
      node$class,
      NumpyArray = {
        form$primitive <- node$primitive
        form$inner_shape <- list()
        emit(paste0(key, "-data"), node$data)
      },
      EmptyArray = NULL,
      NullArray = cli_abort("Internal error: NullArray cannot be encoded"),
      ListOffsetArray = {
        form$offsets <- node$offsets
        emit(paste0(key, "-offsets"), node$offsets_data)
        form[["content"]] <- enc(node[["content"]])
      },
      ListArray = {
        form$starts <- node$starts
        form$stops <- node$stops
        emit(paste0(key, "-starts"), node$starts_data)
        emit(paste0(key, "-stops"), node$stops_data)
        form[["content"]] <- enc(node[["content"]])
      },
      RegularArray = {
        form$size <- node$size
        form[["content"]] <- enc(node[["content"]])
      },
      RecordArray = {
        form$fields <- if (is.null(node$fields)) NULL else as.list(node$fields)
        form$contents <- lapply(node$contents, enc)
      },
      IndexedArray = ,
      IndexedOptionArray = {
        form$index <- node$index
        emit(paste0(key, "-index"), node$index_data)
        form[["content"]] <- enc(node[["content"]])
      },
      ByteMaskedArray = {
        form$mask <- "i8"
        form$valid_when <- node$valid_when
        emit(paste0(key, "-mask"), as.integer(node$mask_data))
        form[["content"]] <- enc(node[["content"]])
      },
      UnmaskedArray = {
        form[["content"]] <- enc(node[["content"]])
      }
    )
    form$parameters <- node$parameters %||% list()
    form$form_key <- key
    form
  }
  form <- enc(layout)
  parsed <- awkward_normalize_form(form, "form", NULL)
  list(
    form = parsed,
    form_json = awkward_form_to_json(form),
    length = layout$length,
    buffers = buffers
  )
}

#' Serialise a (raw, un-normalised) form list to JSON
#' @noRd
awkward_form_to_json <- function(form) {
  fix <- function(node) {
    if (!is.list(node)) {
      return(node)
    }
    if (!is.null(node$class)) {
      node$parameters <- if (length(node$parameters %||% list()) == 0) {
        structure(list(), names = character(0))
      } else {
        node$parameters
      }
      if (!is.null(node[["content"]])) node[["content"]] <- fix(node[["content"]])
      if (!is.null(node$contents)) node$contents <- lapply(node$contents, fix)
      if ("fields" %in% names(node) && is.null(node$fields)) {
        node["fields"] <- list(NULL)
      }
    }
    node
  }
  as.character(jsonlite::toJSON(fix(form), auto_unbox = TRUE, null = "null", digits = NA))
}

# Type strings -----------------------------------------------------------------

#' Human-readable type string, in the style of awkward's `ak.type`
#' @noRd
awkward_type_string <- function(node) {
  kind <- if (is_string_like(node)) string_kind(node) else NULL
  if (!is.null(kind)) {
    return(kind)
  }
  switch(
    node$class,
    NumpyArray = node$primitive,
    EmptyArray = "unknown",
    NullArray = "?unknown",
    ListOffsetArray = ,
    ListArray = paste("var *", awkward_type_string(node[["content"]])),
    RegularArray = paste(node$size, "*", awkward_type_string(node[["content"]])),
    RecordArray = {
      inner <- vapply(node$contents, awkward_type_string, character(1))
      if (is.null(node$fields)) {
        paste0("(", paste(inner, collapse = ", "), ")")
      } else {
        paste0("{", paste0(node$fields, ": ", inner, collapse = ", "), "}")
      }
    },
    IndexedArray = awkward_type_string(node[["content"]]),
    IndexedOptionArray = ,
    ByteMaskedArray = ,
    BitMaskedArray = ,
    UnmaskedArray = {
      inner <- awkward_type_string(node[["content"]])
      if (startsWith(inner, "?")) inner else paste0("?", inner)
    }
  )
}
