# Reading and writing `awkward-array` (0.1.0) elements in H5MU files and Zarr
# stores, directly with rhdf5 / Rarr.
#
# Released anndataR has no reader or writer for this encoding, so mudataR
# handles Awkward elements itself and keeps them away from anndataR (see
# anndatar-compat.R). Both stores share the logical layout: a group with
# `encoding-type`, `encoding-version`, `form` (JSON string) and `length`
# attributes, plus one dataset per buffer named `<form_key>-<role>`. Only the
# buffers required by the form are read; dtypes come from the form.
#
# The codec itself (awkward-form.R, awkward-arrow.R, AwkwardArrowArray.R) is
# adapted from the `awkward-array` support proposed for anndataR (MIT).

.awkward_encoding_versions <- "0.1.0"

# numpy dtype strings understood by Rarr
.awkward_zarr_dtypes <- c(
  bool = "|b1",
  int8 = "|i1",
  uint8 = "|u1",
  int16 = "<i2",
  uint16 = "<u2",
  int32 = "<i4",
  uint32 = "<u4",
  int64 = "<i8",
  uint64 = "<u8",
  float32 = "<f4",
  float64 = "<f8"
)

.awkward_hdf5_dtypes <- c(
  int8 = "H5T_STD_I8LE",
  uint8 = "H5T_STD_U8LE",
  int16 = "H5T_STD_I16LE",
  uint16 = "H5T_STD_U16LE",
  int32 = "H5T_STD_I32LE",
  uint32 = "H5T_STD_U32LE",
  int64 = "H5T_STD_I64LE",
  uint64 = "H5T_STD_U64LE",
  float32 = "H5T_IEEE_F32LE",
  float64 = "H5T_IEEE_F64LE"
)

awkward_check_requires <- function(what, pkgs, call = rlang::caller_env()) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    mudatar_abort(
      "{what} requires the package{?s} {.pkg {missing}}",
      call = call
    )
  }
  invisible(TRUE)
}

check_awkward_version <- function(version, path, kind) {
  if (!identical(version, "0.1.0")) {
    mudatar_abort(
      c(
        "Unsupported {.val awkward-array} encoding version {.val {version}}",
        "i" = "Element: {.path {path}} ({kind})",
        "i" = "Supported versions: {.val {(.awkward_encoding_versions)}}."
      ),
      class = c("mudatar_unsupported_encoding", "mudatar_awkward_error")
    )
  }
}

awkward_form_json <- function(form) {
  if (is.character(form)) {
    return(form)
  }
  # some writers store the form as a JSON object rather than a string
  as.character(jsonlite::toJSON(form, auto_unbox = TRUE, null = "null", digits = NA))
}

awkward_length_attribute <- function(n) {
  if (n <= .Machine$integer.max) as.integer(n) else as.double(n)
}

#' Prepare buffers for writing
#'
#' Validates every buffer against the form before anything is written, so a
#' failure never leaves a partially written element behind.
#' @noRd
awkward_prepare_buffers <- function(value, path, kind) {
  form <- awk_field(value, "form")
  specs <- awkward_buffer_specs(form)
  buffers <- awkward_buffers(value)
  missing <- setdiff(specs$key, names(buffers))
  if (length(missing)) {
    awkward_abort("Missing Awkward buffers {.val {missing}}", path = path)
  }
  # decoding validates offsets, indices and lengths against the form
  awkward_layout(value)

  out <- lapply(seq_len(nrow(specs)), function(i) {
    key <- specs$key[[i]]
    dtype <- specs$dtype[[i]]
    x <- buffers[[key]]
    if (is.factor(x)) {
      x <- as.character(x) == "TRUE"
    }
    x <- awkward_flatten(x)
    if (is.raw(x)) {
      x <- as.integer(x)
    }
    if (dtype == "bool") {
      x <- as.logical(x)
    } else if (dtype %in% c("int64", "uint64", "uint32")) {
      if (inherits(x, "integer64")) {
        if (any(!is.na(x) & abs(x) > bit64::as.integer64(2^53))) {
          awkward_abort(
            "Buffer {.val {key}} has 64-bit integers beyond 2^53, which cannot be written exactly yet",
            path = path
          )
        }
        x <- as.double(x)
      }
      if (kind == "zarr") {
        # Rarr writes 64-bit integer arrays from R integers only
        if (length(x) && max(abs(x)) > .Machine$integer.max) {
          awkward_abort(
            "Buffer {.val {key}} has values beyond the 32-bit range, which cannot be written to Zarr yet",
            path = path
          )
        }
        x <- as.integer(x)
      } else {
        x <- as.double(x)
      }
    } else if (dtype %in% c("float32", "float64")) {
      x <- as.double(x)
    } else {
      x <- as.integer(x)
    }
    list(key = key, dtype = dtype, value = x)
  })
  names(out) <- specs$key
  out
}

# Detection --------------------------------------------------------------------

#' Encoding type of an element (or `NULL`)
#' @noRd
store_encoding_type <- function(kind, path, name) {
  attrs <- tryCatch(store_read_attrs(kind, path, name), error = function(e) NULL)
  attrs[["encoding-type"]]
}

#' Names of the direct children of a group that are Awkward Arrays
#' @noRd
store_awkward_children <- function(kind, path, name) {
  kids <- store_list_children(kind, path, name)
  kids[vapply(kids, function(k) {
    identical(store_encoding_type(kind, path, paste0(name, "/", k)), "awkward-array")
  }, logical(1))]
}

# Reading ----------------------------------------------------------------------

#' Read an `awkward-array` element as an [AwkwardArrowArray]
#'
#' @param lazy Whether buffers are read on first use (HDF5 only)
#' @noRd
store_read_awkward <- function(kind, path, name, lazy = FALSE) {
  awkward_check_requires("Reading Awkward Arrays", c("jsonlite", "nanoarrow", "bit64"))
  name <- sub("^/", "", name)
  label <- paste0(path, "/", name)
  attrs <- store_read_attrs(kind, path, name)
  check_awkward_version(attrs[["encoding-version"]], label, kind)
  form_json <- awkward_form_json(attrs[["form"]])
  form <- awkward_parse_form(form_json, path = label)
  keys <- awkward_buffer_specs(form)$key

  load <- function() {
    present <- store_list_children(kind, path, name)
    missing <- setdiff(keys, present)
    if (length(missing)) {
      awkward_abort("Missing Awkward buffers {.val {missing}}", path = label)
    }
    buffers <- lapply(keys, function(key) {
      if (kind == "h5mu") {
        rhdf5::h5read(path, paste0(name, "/", key), bit64conversion = "bit64", native = FALSE)
      } else {
        as.vector(Rarr::read_zarr_array(file.path(path, name, key)))
      }
    })
    names(buffers) <- keys
    buffers
  }

  lazy <- isTRUE(lazy) && kind == "h5mu"
  new_awkward_arrow_array(
    form = form,
    form_json = form_json,
    length = as.double(attrs[["length"]]),
    buffers = if (lazy) NULL else load(),
    loader = if (lazy) load else NULL,
    source = list(store_kind = kind, path = label, encoding_version = "0.1.0")
  )
}

# Writing ----------------------------------------------------------------------

#' Write an [AwkwardArrowArray] as an `awkward-array` element
#'
#' Replaces an existing element of the same name.
#' @noRd
store_write_awkward <- function(kind, path, name, value, zarr_format = NULL, compression = NULL) {
  name <- sub("^/", "", name)
  label <- paste0(path, "/", name)
  # validate before anything on disk changes
  buffers <- awkward_prepare_buffers(value, label, kind)
  store_delete(kind, path, name)

  if (kind == "h5mu") {
    fid <- h5_open(path, readonly = FALSE)
    on.exit(rhdf5::H5Fclose(fid), add = TRUE)
    rhdf5::H5Gclose(rhdf5::H5Gcreate(fid, name))
    for (buf in buffers) {
      dataset <- paste0(name, "/", buf$key)
      if (buf$dtype == "bool") {
        h5_write_boolean(fid, dataset, buf$value)
      } else {
        rhdf5::h5createDataset(
          fid,
          dataset,
          dims = length(buf$value),
          storage.mode = storage.mode(buf$value),
          H5type = .awkward_hdf5_dtypes[[buf$dtype]],
          chunk = NULL,
          level = 0,
          native = FALSE
        )
        if (length(buf$value)) rhdf5::h5write(buf$value, fid, dataset, native = FALSE)
      }
    }
    rhdf5::H5Fclose(fid)
    on.exit()
  } else {
    zarr_format <- zarr_format %||% store_zarr_format(path)
    Rarr::write_zarr_group(path, name, zarr_version = zarr_format)
    for (buf in buffers) {
      n <- length(buf$value)
      Rarr::write_zarr_array(
        array(buf$value, dim = n),
        zarr_array_path = file.path(path, name, buf$key),
        chunk_dim = max(n, 1L),
        data_type = .awkward_zarr_dtypes[[buf$dtype]],
        compressor = zarr_compressor(compression),
        zarr_version = zarr_format
      )
    }
  }
  store_write_attrs(
    kind,
    path,
    name,
    list(
      `encoding-type` = "awkward-array",
      `encoding-version` = "0.1.0",
      form = awk_field(value, "form_json"),
      length = awkward_length_attribute(length(value))
    ),
    zarr_format = zarr_format
  )
  invisible(NULL)
}

#' Write a boolean dataset the way h5py does (an 8-bit FALSE/TRUE enum)
#' @noRd
h5_write_boolean <- function(h5loc, name, value) {
  space <- rhdf5::H5Screate_simple(dims = length(value), NULL, native = FALSE)
  on.exit(rhdf5::H5Sclose(space), add = TRUE)
  tid <- rhdf5::H5Tenum_create(dtype_id = "H5T_NATIVE_SCHAR")
  rhdf5::H5Tenum_insert(tid, name = "FALSE", value = 0L)
  rhdf5::H5Tenum_insert(tid, name = "TRUE", value = 1L)
  dcpl <- rhdf5::H5Pcreate("H5P_DATASET_CREATE")
  on.exit(rhdf5::H5Pclose(dcpl), add = TRUE)
  rhdf5::H5Pset_fill_time(dcpl, "H5D_FILL_TIME_ALLOC")
  rhdf5::H5Pset_obj_track_times(dcpl, FALSE)
  did <- rhdf5::H5Dcreate(h5loc = h5loc, name = name, dtype_id = tid, h5space = space, dcpl = dcpl)
  on.exit(rhdf5::H5Dclose(did), add = TRUE)
  if (length(value)) rhdf5::H5Dwrite(did, as.raw(as.integer(value)), h5type = tid)
  invisible(NULL)
}

zarr_compressor <- function(compression) {
  switch(
    compression %||% "none",
    "zstd" = Rarr::use_zstd(),
    "gzip" = Rarr::use_gzip(),
    "blosc" = Rarr::use_blosc(),
    "lzma" = Rarr::use_lzma(),
    "bz2" = Rarr::use_bz2(),
    "zlib" = Rarr::use_zlib(),
    "lz4" = Rarr::use_lz4(),
    NULL
  )
}
