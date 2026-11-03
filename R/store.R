# Store-level helpers for H5MU files and MuData Zarr stores
#
# Element *encodings* (dataframes, arrays, mappings, ...) are anndataR's,
# reached through its public API (see anndatar-compat.R); Awkward Arrays are
# read and written by mudataR (awkward-store.R). This file handles the
# store primitives and what the MuData format adds on top: root/`mod`
# attributes, the modality groups, the typed maps and the H5MU user block.

store_kind <- function(path) {
  if (dir.exists(path)) "zarr" else "h5mu"
}

store_check_requires <- function(kind, call = rlang::caller_env()) {
  pkgs <- if (kind == "h5mu") "rhdf5" else "Rarr"
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    mudatar_abort(
      "Package{?s} {.pkg {missing}} {?is/are} required for {if (kind == 'h5mu') 'H5MU files' else 'Zarr stores'}",
      call = call
    )
  }
}

h5_open <- function(path, readonly = TRUE) {
  rhdf5::H5Fopen(path, flags = if (readonly) "H5F_ACC_RDONLY" else "H5F_ACC_RDWR")
}

store_read_attrs <- function(kind, path, name = "/") {
  if (kind == "h5mu") {
    rhdf5::h5readAttributes(path, name)
  } else {
    Rarr::read_zarr_attributes(if (name == "/") path else file.path(path, sub("^/", "", name)))
  }
}

#' Write attributes; character vectors wrapped in `I()` are always arrays
#' @noRd
store_write_attrs <- function(kind, path, name, attrs, zarr_format = NULL) {
  if (kind == "h5mu") {
    fid <- h5_open(path, readonly = FALSE)
    on.exit(rhdf5::H5Fclose(fid), add = TRUE)
    obj <- if (name == "/") fid else rhdf5::H5Oopen(fid, name)
    if (name != "/") on.exit(rhdf5::H5Oclose(obj), add = TRUE)
    for (a in names(attrs)) {
      value <- attrs[[a]]
      if (rhdf5::H5Aexists(obj, a)) rhdf5::H5Adelete(obj, a)
      is_array <- inherits(value, "AsIs") || length(value) != 1L
      value <- unclass(value)
      if (is.character(value)) {
        rhdf5::h5writeAttribute(
          enc2utf8(value),
          h5obj = obj,
          name = a,
          asScalar = !is_array,
          variableLengthString = TRUE,
          encoding = "UTF-8"
        )
      } else {
        rhdf5::h5writeAttribute(value, h5obj = obj, name = a, asScalar = !is_array)
      }
    }
  } else {
    zpath <- if (name == "/") path else file.path(path, sub("^/", "", name))
    attrs <- lapply(attrs, function(v) if (inherits(v, "AsIs") || length(v) != 1L) as.list(unclass(v)) else v)
    Rarr::write_zarr_attributes(zpath, new.zattrs = attrs, zarr_version = zarr_format)
  }
  invisible(NULL)
}

store_list_groups <- function(kind, path, name) {
  kids <- store_list_children(kind, path, name, info = TRUE)
  sort(kids$name[kids$group])
}

#' Direct children (groups and arrays/datasets) of a group, sorted by name
#' @noRd
store_list_children <- function(kind, path, name, info = FALSE) {
  name <- sub("^/", "", name)
  empty <- data.frame(name = character(), group = logical())
  if (kind == "h5mu") {
    fid <- h5_open(path)
    on.exit(rhdf5::H5Fclose(fid), add = TRUE)
    if (nzchar(name) && !rhdf5::H5Lexists(fid, name)) {
      return(if (info) empty else character())
    }
    gid <- if (nzchar(name)) rhdf5::H5Gopen(fid, name) else fid
    if (nzchar(name)) on.exit(rhdf5::H5Gclose(gid), add = TRUE, after = FALSE)
    ls <- rhdf5::h5ls(gid, recursive = FALSE, datasetinfo = FALSE)
    out <- data.frame(name = ls$name, group = ls$otype == "H5I_GROUP")
  } else {
    dir <- if (nzchar(name)) file.path(path, name) else path
    if (!dir.exists(dir)) {
      return(if (info) empty else character())
    }
    subs <- list.dirs(dir, recursive = FALSE, full.names = FALSE)
    kind_of <- vapply(subs, function(s) {
      files <- list.files(file.path(dir, s), all.files = TRUE)
      if (".zarray" %in% files) {
        "array"
      } else if (any(c(".zgroup", ".zattrs") %in% files)) {
        "group"
      } else if ("zarr.json" %in% files) {
        meta <- jsonlite::read_json(file.path(dir, s, "zarr.json"))
        if (identical(meta$node_type, "array")) "array" else "group"
      } else {
        NA_character_
      }
    }, character(1))
    keep <- !is.na(kind_of)
    out <- data.frame(name = subs[keep], group = unname(kind_of[keep] == "group"))
  }
  out <- out[order(out$name), , drop = FALSE]
  if (info) out else out$name
}

store_exists <- function(kind, path, name) {
  name <- sub("^/", "", name)
  if (kind == "h5mu") {
    fid <- h5_open(path)
    on.exit(rhdf5::H5Fclose(fid), add = TRUE)
    # check each level: H5Lexists fails on a missing intermediate group
    parts <- strsplit(name, "/", fixed = TRUE)[[1]]
    for (i in seq_along(parts)) {
      if (!rhdf5::H5Lexists(fid, paste(parts[seq_len(i)], collapse = "/"))) {
        return(FALSE)
      }
    }
    TRUE
  } else {
    dir.exists(file.path(path, name))
  }
}

#' Delete an element if it exists
#' @noRd
store_delete <- function(kind, path, name) {
  name <- sub("^/", "", name)
  if (!store_exists(kind, path, name)) {
    return(invisible(FALSE))
  }
  if (kind == "h5mu") {
    fid <- h5_open(path, readonly = FALSE)
    on.exit(rhdf5::H5Fclose(fid), add = TRUE)
    rhdf5::H5Ldelete(fid, name)
  } else {
    unlink(file.path(path, name), recursive = TRUE)
  }
  invisible(TRUE)
}

#' Create an empty `dict` group (replacing an existing element)
#' @noRd
store_create_dict <- function(kind, path, name, zarr_format = NULL) {
  name <- sub("^/", "", name)
  store_delete(kind, path, name)
  if (kind == "h5mu") {
    rhdf5::h5createGroup(path, name)
  } else {
    zarr_format <- zarr_format %||% store_zarr_format(path)
    Rarr::write_zarr_group(path, name, zarr_version = zarr_format)
  }
  store_write_attrs(
    kind,
    path,
    name,
    list(`encoding-type` = "dict", `encoding-version` = "0.1.0"),
    zarr_format = zarr_format
  )
}

store_zarr_format <- function(path) {
  files <- list.files(path, all.files = TRUE)
  if ("zarr.json" %in% files) 3L else 2L
}

#' Create an empty container with MuData root attributes
#' @noRd
store_create <- function(kind, path, axis, zarr_format = 3L) {
  if (kind == "h5mu") {
    # Reserve the 512-byte user block that marks a file as H5MU; the marker
    # itself is written by store_finalize() once the HDF5 file is closed.
    fcpl <- rhdf5::H5Pcreate("H5P_FILE_CREATE")
    on.exit(rhdf5::H5Pclose(fcpl), add = TRUE)
    rhdf5::H5Pset_userblock(fcpl, .h5mu_userblock_size)
    fid <- rhdf5::H5Fcreate(path, flags = "H5F_ACC_EXCL", fcpl = fcpl)
    rhdf5::H5Fclose(fid)
  } else {
    Rarr::write_zarr_group(path, "", zarr_version = zarr_format)
  }
  store_write_attrs(
    kind,
    path,
    "/",
    list(
      `encoding-type` = .mudata_encoding_type,
      `encoding-version` = .mudata_encoding_version,
      axis = as.integer(axis),
      encoder = "mudataR",
      `encoder-version` = mudatar_version()
    ),
    zarr_format = zarr_format
  )
}

h5mu_marker <- function() {
  paste0(
    "MuData (format-version=",
    .mudata_format_version,
    ";creator=mudataR;creator-version=",
    mudatar_version(),
    ")"
  )
}

#' Write the H5MU marker into the reserved user block
#'
#' The user block is space HDF5 reserves at the start of the file and never
#' touches; Python mudata writes the same marker the same way.
#' @noRd
store_finalize <- function(kind, path) {
  if (kind != "h5mu") {
    return(invisible(NULL))
  }
  fid <- h5_open(path)
  fcpl <- rhdf5::H5Fget_create_plist(fid)
  size <- rhdf5::H5Pget_userblock(fcpl)
  rhdf5::H5Pclose(fcpl)
  rhdf5::H5Fclose(fid)
  if (size < .h5mu_userblock_size) {
    mudatar_abort("Internal error: H5MU user block was not reserved")
  }
  marker <- charToRaw(enc2utf8(h5mu_marker()))
  con <- file(path, open = "r+b")
  on.exit(close(con), add = TRUE)
  seek(con, 0, rw = "write")
  writeBin(c(marker, raw(.h5mu_userblock_size - length(marker))), con)
  invisible(NULL)
}

#' Read the H5MU user block marker (or `NULL`)
#' @noRd
h5mu_read_marker <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  bytes <- readBin(con, "raw", n = .h5mu_userblock_size)
  if (length(bytes) < 6 || rawToChar(bytes[1:6]) != "MuData") {
    return(NULL)
  }
  rawToChar(bytes[seq_len(max(which(bytes != 0)))])
}

#' Write maps as unsigned 32-bit arrays, as Python mudata does
#' @noRd
store_write_maps <- function(kind, path, attr, maps, zarr_format = NULL) {
  name <- paste0(attr, "map")
  store_create_dict(kind, path, name, zarr_format)
  for (m in names(maps)) {
    values <- as.integer(maps[[m]])
    dataset <- paste0(name, "/", m)
    if (kind == "h5mu") {
      rhdf5::h5createDataset(
        path,
        dataset,
        dims = length(values),
        storage.mode = "integer",
        H5type = "H5T_STD_U32LE",
        chunk = NULL,
        level = 0
      )
      if (length(values)) rhdf5::h5write(values, path, dataset)
    } else {
      Rarr::write_zarr_array(
        array(values, dim = length(values)),
        zarr_array_path = file.path(path, dataset),
        chunk_dim = max(length(values), 1L),
        data_type = "<u4",
        zarr_version = zarr_format
      )
    }
    store_write_attrs(
      kind,
      path,
      dataset,
      list(`encoding-type` = "array", `encoding-version` = "0.2.0"),
      zarr_format = zarr_format
    )
  }
}

store_read_maps <- function(kind, path, attr) {
  name <- paste0(attr, "map")
  keys <- store_list_children(kind, path, name)
  maps <- lapply(keys, function(k) {
    value <- if (kind == "h5mu") {
      rhdf5::h5read(path, paste0(name, "/", k))
    } else {
      Rarr::read_zarr_array(file.path(path, name, k))
    }
    as.integer(as.vector(value))
  })
  names(maps) <- keys
  maps
}

#' Ordered modality names: `mod-order` when it lists exactly existing
#' modalities (RFC), otherwise the stored group order
#' @noRd
store_mod_names <- function(kind, path) {
  groups <- store_list_groups(kind, path, "/mod")
  if (!length(groups)) {
    return(character())
  }
  attrs <- store_read_attrs(kind, path, "/mod")
  order <- as.vector(unlist(attrs[["mod-order"]]))
  if (!is.null(order) && all(order %in% groups)) {
    # as Python mudata: the listed order (and only the listed modalities)
    order
  } else {
    groups
  }
}
