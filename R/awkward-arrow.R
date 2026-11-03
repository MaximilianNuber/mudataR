# Conversion between Awkward layouts and Arrow (via nanoarrow)
#
# Layout -> Arrow mapping:
# nolint start: commented_code_linter
#   NumpyArray                -> primitive (bool, int*, uint*, float, double)
#   ListOffsetArray/ListArray -> large_list / list (by offset width)
#     with __array__="string" -> large_string / string
#     with __array__="bytestring" -> large_binary / binary
#   RegularArray              -> fixed_size_list
#   RecordArray               -> struct (tuples get fields "0", "1", ...)
#   option nodes              -> validity bitmap on the content's Arrow array
#   option of EmptyArray      -> null type
# nolint end: commented_code_linter
#
# Arrow fields are marked nullable exactly when the Awkward node is an option
# type, so the Arrow schema preserves the distinction between `string` and
# `?string`. The reverse conversion uses that flag to decide whether a node
# is written as an option.

.na_nullable_flag <- 2L

na_set_nullable <- function(schema, nullable) {
  nanoarrow::nanoarrow_schema_modify(
    schema,
    list(flags = if (nullable) .na_nullable_flag else 0L)
  )
}

na_set_name <- function(schema, name) {
  nanoarrow::nanoarrow_schema_modify(schema, list(name = name))
}

# list-like schema constructors drop the child name; Arrow convention is "item"
na_set_child_item <- function(schema, child_schema) {
  nanoarrow::nanoarrow_schema_modify(schema, list(children = list(item = child_schema)))
}

na_validity_buffer <- function(valid) {
  if (is.null(valid)) {
    return(NULL)
  }
  nanoarrow::as_nanoarrow_array(valid)$buffers[[2]]
}

na_index_buffer <- function(x, large = TRUE) {
  nanoarrow::as_nanoarrow_array(
    x,
    schema = if (large) nanoarrow::na_int64() else nanoarrow::na_int32()
  )$buffers[[2]]
}

na_primitive_schema <- function(primitive) {
  switch(
    primitive,
    bool = nanoarrow::na_bool(),
    int8 = nanoarrow::na_int8(),
    uint8 = nanoarrow::na_uint8(),
    int16 = nanoarrow::na_int16(),
    uint16 = nanoarrow::na_uint16(),
    int32 = nanoarrow::na_int32(),
    uint32 = nanoarrow::na_uint32(),
    int64 = nanoarrow::na_int64(),
    uint64 = nanoarrow::na_uint64(),
    float32 = nanoarrow::na_float(),
    float64 = nanoarrow::na_double()
  )
}

na_primitive_data_buffer <- function(x, primitive) {
  if (primitive == "bool") {
    x[is.na(x)] <- FALSE
    return(nanoarrow::as_nanoarrow_array(x)$buffers[[2]])
  }
  if (primitive == "float64") {
    # writeBin keeps NaN payloads and signed zeros bit-exact
    return(nanoarrow::as_nanoarrow_buffer(writeBin(as.double(x), raw(), size = 8L, endian = "little")))
  }
  if (anyNA(x)) {
    # placeholders below a null parent
    x[is.na(x)] <- 0
  }
  nanoarrow::as_nanoarrow_array(x, schema = na_primitive_schema(primitive))$buffers[[2]]
}

na_build <- function(schema, length, valid, buffers, children = NULL) {
  null_count <- if (is.null(valid)) 0 else sum(!valid)
  args <- list(
    length = length,
    null_count = null_count,
    buffers = c(list(na_validity_buffer(valid)), buffers)
  )
  if (!is.null(children)) {
    args$children <- children
  }
  nanoarrow::nanoarrow_array_modify(
    nanoarrow::nanoarrow_array_init(schema),
    args,
    validate = TRUE
  )
}

#' Convert a layout to a nanoarrow array
#'
#' @param node A layout
#' @param valid `NULL` (non-option) or a logical validity vector of length
#'   `node$length` contributed by enclosing option nodes
#' @noRd
awkward_layout_to_nanoarrow <- function(node, valid = NULL) {
  n <- node$length
  nullable <- !is.null(valid)
  switch(
    node$class,
    NumpyArray = {
      schema <- na_set_nullable(na_primitive_schema(node$primitive), nullable)
      data <- if (is.raw(node$data)) as.integer(node$data) else node$data
      na_build(schema, n, valid, list(na_primitive_data_buffer(data, node$primitive)))
    },
    EmptyArray = ,
    NullArray = {
      nanoarrow::nanoarrow_array_modify(
        nanoarrow::nanoarrow_array_init(nanoarrow::na_na()),
        list(length = n, null_count = n)
      )
    },
    ListArray = awkward_layout_to_nanoarrow(
      awkward_take(node, seq_len(n) - 1),
      valid
    ),
    ListOffsetArray = {
      large <- node$offsets != "i32"
      offsets <- node$offsets_data
      if (!large && max(offsets, 0) > .Machine$integer.max) {
        large <- TRUE
      }
      kind <- string_kind(node)
      content <- node[["content"]]
      if (!is.null(kind) && content$class == "NumpyArray" && is.raw(content$data)) {
        schema <- if (kind == "string") {
          if (large) nanoarrow::na_large_string() else nanoarrow::na_string()
        } else {
          if (large) nanoarrow::na_large_binary() else nanoarrow::na_binary()
        }
        schema <- na_set_nullable(schema, nullable)
        return(na_build(
          schema,
          n,
          valid,
          list(na_index_buffer(offsets, large), nanoarrow::as_nanoarrow_buffer(content$data))
        ))
      }
      child <- awkward_layout_to_nanoarrow(content)
      child_schema <- na_set_name(nanoarrow::infer_nanoarrow_schema(child), "item")
      schema <- if (large) {
        nanoarrow::na_large_list(child_schema)
      } else {
        nanoarrow::na_list(child_schema)
      }
      schema <- na_set_child_item(na_set_nullable(schema, nullable), child_schema)
      na_build(schema, n, valid, list(na_index_buffer(offsets, large)), list(item = child))
    },
    RegularArray = {
      child <- awkward_layout_to_nanoarrow(node[["content"]])
      child_schema <- na_set_name(nanoarrow::infer_nanoarrow_schema(child), "item")
      schema <- na_set_child_item(
        na_set_nullable(
          nanoarrow::na_fixed_size_list(child_schema, as.integer(node$size)),
          nullable
        ),
        child_schema
      )
      na_build(schema, n, valid, list(), list(item = child))
    },
    RecordArray = {
      fields <- node$fields %||% as.character(seq_along(node$contents) - 1L)
      children <- lapply(node$contents, awkward_layout_to_nanoarrow)
      names(children) <- fields
      child_schemas <- lapply(children, nanoarrow::infer_nanoarrow_schema)
      schema <- na_set_nullable(nanoarrow::na_struct(child_schemas), nullable)
      na_build(schema, n, valid, list(), children)
    },
    IndexedArray = awkward_layout_to_nanoarrow(
      awkward_take(node[["content"]], node$index_data),
      valid
    ),
    IndexedOptionArray = {
      ok <- node$index_data >= 0
      idx <- node$index_data
      idx[!ok] <- NA
      awkward_layout_to_nanoarrow(
        awkward_take(node[["content"]], idx),
        if (is.null(valid)) ok else valid & ok
      )
    },
    ByteMaskedArray = {
      ok <- node$mask_data == node$valid_when
      awkward_layout_to_nanoarrow(node[["content"]], if (is.null(valid)) ok else valid & ok)
    },
    UnmaskedArray = {
      awkward_layout_to_nanoarrow(node[["content"]], valid %||% rep(TRUE, n))
    }
  )
}

# Arrow -> layout --------------------------------------------------------------

na_buffer_raw <- function(buffer) {
  if (is.null(buffer)) {
    return(NULL)
  }
  as.raw(buffer)
}

na_validity <- function(array, start, n) {
  if (n == 0) {
    return(logical(0))
  }
  bytes <- na_buffer_raw(array$buffers[[1]])
  if (array$null_count == 0 || length(bytes) == 0) {
    return(rep(TRUE, n))
  }
  pos <- array$offset + start + seq_len(n)
  as.logical(rawToBits(bytes))[pos]
}

na_offsets <- function(array, start, n, large) {
  buffer_array <- nanoarrow::nanoarrow_array_modify(
    nanoarrow::nanoarrow_array_init(if (large) nanoarrow::na_int64() else nanoarrow::na_int32()),
    list(
      length = array$offset + start + n + 1,
      null_count = 0,
      buffers = list(NULL, array$buffers[[2]])
    )
  )
  offsets <- as.double(nanoarrow::convert_array(buffer_array, double()))
  offsets[array$offset + start + seq_len(n + 1)]
}

na_flag_nullable <- function(schema) {
  bitwAnd(as.integer(schema$flags), .na_nullable_flag) != 0L
}

.na_primitive_formats <- c(
  b = "bool",
  c = "int8",
  C = "uint8",
  s = "int16",
  S = "uint16",
  i = "int32",
  I = "uint32",
  l = "int64",
  L = "uint64",
  f = "float32",
  g = "float64"
)

#' Convert a nanoarrow array (logical slice) to a layout
#'
#' @param array A nanoarrow array (with its schema)
#' @param schema Its schema
#' @param start,n Logical slice of the array (before `array$offset`)
#' @param where Path for error messages
#' @noRd
awkward_layout_from_nanoarrow <- function(
  array,
  schema = nanoarrow::infer_nanoarrow_schema(array),
  start = 0,
  n = array$length - start,
  where = "value"
) {
  fmt <- schema$format
  nullable <- na_flag_nullable(schema) || array$null_count != 0
  valid <- if (nullable) na_validity(array, start, n) else NULL

  node <- if (fmt %in% names(.na_primitive_formats)) {
    primitive <- .na_primitive_formats[[fmt]]
    sliced <- nanoarrow::nanoarrow_array_modify(
      array,
      list(offset = array$offset + start, length = n),
      validate = FALSE
    )
    data <- nanoarrow::convert_array(sliced)
    if (primitive %in% c("int64", "uint64")) {
      data <- as.double(data)
    }
    if (primitive == "float32") {
      data <- as.double(data)
    }
    list(
      class = "NumpyArray",
      primitive = primitive,
      inner_shape = list(),
      parameters = list(),
      length = n,
      data = data
    )
  } else if (fmt %in% c("u", "U", "z", "Z")) {
    large <- fmt %in% c("U", "Z")
    offsets <- na_offsets(array, start, n, large)
    bytes <- na_buffer_raw(array$buffers[[3]]) %||% raw(0)
    first <- offsets[[1]]
    last <- offsets[[n + 1]]
    list(
      class = "ListOffsetArray",
      offsets = if (large) "i64" else "i32",
      parameters = list(`__array__` = if (fmt %in% c("u", "U")) "string" else "bytestring"),
      length = n,
      offsets_data = offsets - first,
      content = list(
        class = "NumpyArray",
        primitive = "uint8",
        inner_shape = list(),
        parameters = list(`__array__` = if (fmt %in% c("u", "U")) "char" else "byte"),
        length = last - first,
        data = if (last > first) bytes[(first + 1):last] else raw(0)
      )
    )
  } else if (fmt %in% c("+l", "+L")) {
    large <- fmt == "+L"
    offsets <- na_offsets(array, start, n, large)
    first <- offsets[[1]]
    last <- offsets[[n + 1]]
    list(
      class = "ListOffsetArray",
      offsets = if (large) "i64" else "i32",
      parameters = list(),
      length = n,
      offsets_data = offsets - first,
      content = awkward_layout_from_nanoarrow(
        array$children[[1]],
        schema$children[[1]],
        start = first,
        n = last - first,
        where = paste0(where, "[]")
      )
    )
  } else if (startsWith(fmt, "+w:")) {
    size <- as.numeric(sub("^\\+w:", "", fmt))
    list(
      class = "RegularArray",
      size = size,
      parameters = list(),
      length = n,
      content = awkward_layout_from_nanoarrow(
        array$children[[1]],
        schema$children[[1]],
        start = (array$offset + start) * size,
        n = n * size,
        where = paste0(where, "[]")
      )
    )
  } else if (fmt == "+s") {
    fields <- names(schema$children)
    list(
      class = "RecordArray",
      fields = fields,
      parameters = list(),
      length = n,
      contents = lapply(seq_along(fields), function(i) {
        awkward_layout_from_nanoarrow(
          array$children[[i]],
          schema$children[[i]],
          start = array$offset + start,
          n = n,
          where = paste0(where, "$", fields[[i]])
        )
      })
    )
  } else if (fmt == "n") {
    # all-null: option of EmptyArray, as awkward writes `?unknown`
    return(list(
      class = "IndexedOptionArray",
      index = "i64",
      parameters = list(),
      length = n,
      index_data = rep(-1, n),
      content = list(class = "EmptyArray", parameters = list(), length = 0)
    ))
  } else {
    cli_abort(
      c(
        "Cannot encode Arrow type {.val {fmt}} at {.field {where}} as an Awkward Array",
        "i" = "Supported: booleans, integers, floats, (large) strings/binary, (large/fixed-size) lists and structs."
      ),
      class = c("mudatar_unsupported_encoding", "mudatar_awkward_error")
    )
  }

  if (is.null(valid)) {
    return(node)
  }
  # Option node: IndexedOptionArray over the compacted valid content, which
  # is how awkward itself represents missing values from Python objects
  index <- rep(-1, n)
  index[valid] <- seq_len(sum(valid)) - 1
  list(
    class = "IndexedOptionArray",
    index = "i64",
    parameters = list(),
    length = n,
    index_data = index,
    content = awkward_take(node, which(valid) - 1)
  )
}
