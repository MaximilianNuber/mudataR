#' Read a MuData H5MU file
#'
#' Read an H5MU file natively (no Python), either fully into memory or as a
#' lazily-read [HDF5MuData].
#'
#' @param path Path to the `.h5mu` file
#' @param as `"InMemoryMuData"` (default) to read everything into memory, or
#'   `"HDF5MuData"` for an object that reads on access
#' @param mode For `as = "HDF5MuData"`: `"r"` (default) or `"r+"`
#' @param backed For `as = "HDF5MuData"`: whether arrays are returned as
#'   `DelayedArray` objects and Awkward Array buffers (see [get_awkward()])
#'   are loaded lazily
#' @param modalities Optional character vector of modalities to expose; the
#'   others are not read (and not modified)
#'
#' @return An [InMemoryMuData] or [HDF5MuData]
#' @family MuData I/O
#' @export
read_h5mu <- function(
  path,
  as = c("InMemoryMuData", "HDF5MuData"),
  mode = c("r", "r+"),
  backed = FALSE,
  modalities = NULL
) {
  as <- match.arg(as)
  mode <- match.arg(mode)
  if (as == "InMemoryMuData" && (mode != "r" || isTRUE(backed))) {
    mudatar_abort("{.arg mode} and {.arg backed} only apply to {.code as = \"HDF5MuData\"}")
  }
  store_check_requires("h5mu")
  if (file.exists(path) && !dir.exists(path) && is.null(h5mu_read_marker(path))) {
    cli_warn(
      c(
        "{.file {path}} has no MuData user block marker.",
        "i" = "It was not written by mudata, muon or mudataR; reading may not work correctly."
      ),
      class = "mudatar_missing_marker"
    )
  }
  # nolint start: object_name_linter
  h5 <- HDF5MuData$new(path, mode = mode, backed = backed, modalities = modalities)
  # nolint end: object_name_linter
  if (as == "HDF5MuData") {
    return(h5)
  }
  on.exit(h5$close(), add = TRUE)
  h5$as_InMemoryMuData()
}

#' Read a MuData Zarr store
#'
#' Read a MuData Zarr store (v2 or v3) natively, either into memory or as a
#' [ZarrMuData] that reads on access.
#'
#' @inheritParams read_h5mu
#' @param path Path to the Zarr store
#' @param as `"InMemoryMuData"` (default) or `"ZarrMuData"`
#' @param mode For `as = "ZarrMuData"`: `"r"` (default) or `"r+"`
#'
#' @return An [InMemoryMuData] or [ZarrMuData]
#' @family MuData I/O
#' @export
read_mudata_zarr <- function(
  path,
  as = c("InMemoryMuData", "ZarrMuData"),
  mode = c("r", "r+"),
  modalities = NULL
) {
  as <- match.arg(as)
  mode <- match.arg(mode)
  if (as == "InMemoryMuData" && mode != "r") {
    mudatar_abort("{.arg mode} only applies to {.code as = \"ZarrMuData\"}")
  }
  store_check_requires("zarr")
  # nolint start: object_name_linter
  z <- ZarrMuData$new(path, mode = mode, modalities = modalities)
  # nolint end: object_name_linter
  if (as == "ZarrMuData") {
    return(z)
  }
  on.exit(z$close(), add = TRUE)
  z$as_InMemoryMuData()
}

#' Write a MuData object to an H5MU file
#'
#' Writes natively (no Python) a file readable by Python `mudata`, including
#' the 512-byte MuData user block, root attributes, `mod-order`, unsigned
#' 32-bit maps and each modality as an AnnData group under `/mod/<name>`.
#'
#' The object is validated and the file is written to a temporary sibling
#' path that is renamed over `path` only on success, so a failed write never
#' leaves a truncated file behind.
#'
#' @param x A MuData object ([InMemoryMuData], [HDF5MuData], [ZarrMuData])
#' @param path Output path
#' @param mode `"w-"`/`"x"` (default; fail if `path` exists) or `"w"`
#'   (replace)
#' @param compression Compression for datasets, passed to anndataR
#'   (`"none"`, `"gzip"` or `"lzf"`)
#'
#' @return `path`, invisibly
#' @family MuData I/O
#' @export
write_h5mu <- function(x, path, mode = c("w-", "w", "x"), compression = c("none", "gzip", "lzf")) {
  mode <- match.arg(mode)
  compression <- match.arg(compression)
  write_mudata_store(x, path, kind = "h5mu", mode = mode, compression = compression)
}

#' Write a MuData object to a Zarr store
#'
#' @inheritParams write_h5mu
#' @param compression Zarr codec passed to anndataR (default `"zstd"`)
#' @param zarr_format Zarr format, `2` or `3` (default)
#'
#' @return `path`, invisibly
#' @family MuData I/O
#' @export
write_mudata_zarr <- function(x, path, mode = c("w-", "w", "x"), compression = "zstd", zarr_format = 3L) {
  mode <- match.arg(mode)
  zarr_format <- as.integer(zarr_format)
  if (!zarr_format %in% c(2L, 3L)) {
    mudatar_abort("{.arg zarr_format} must be 2 or 3")
  }
  write_mudata_store(x, path, kind = "zarr", mode = mode, compression = compression, zarr_format = zarr_format)
}

write_mudata_store <- function(x, path, kind, mode, compression = NULL, zarr_format = 3L, call = rlang::caller_env()) {
  if (!inherits(x, "AbstractMuData")) {
    mudatar_abort("{.arg x} must be a MuData object, not {.cls {class(x)}}", call = call)
  }
  store_check_requires(kind, call = call)
  exists <- file.exists(path)
  if (exists && mode %in% c("w-", "x")) {
    mudatar_abort(
      c("{.file {path}} already exists.", "i" = "Use {.code mode = \"w\"} to replace it."),
      call = call
    )
  }
  if (!dir.exists(dirname(path))) {
    mudatar_abort("Directory {.file {dirname(path)}} does not exist", call = call)
  }

  # Validate the complete object before anything is written
  mods <- as.list(x$mod)
  mod_names <- names(mods)
  obs_names <- x$obs_names
  var_names <- x$var_names
  for (attr in c("obs", "var")) {
    names <- if (attr == "obs") obs_names else var_names
    if (anyDuplicated(names)) {
      mudatar_abort_unsupported("Writing duplicated global {attr}_names", call = call)
    }
    maps <- validate_maps(
      .subset2(x, paste0(attr, "map")),
      attr,
      names,
      vapply(mods, function(a) length(a[[paste0(attr, "_names")]]), integer(1)),
      call = call
    )
    check_maps_match_names(maps, attr, names, lapply(mods, function(a) a[[paste0(attr, "_names")]]), call = call)
  }
  globals <- list(
    obs = x$obs,
    var = x$var,
    obsm = x$obsm,
    varm = x$varm,
    obsp = x$obsp,
    varp = x$varp,
    uns = x$uns
  )
  rownames(globals$obs) <- obs_names
  rownames(globals$var) <- var_names
  maps <- list(obs = x$obsmap, var = x$varmap)
  axis <- x$axis

  tmp <- tempfile(
    pattern = paste0(".", basename(path), "-"),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  done <- FALSE
  on.exit(if (!done) unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

  store_create(kind, tmp, axis, zarr_format)
  for (slot in names(globals)) {
    store_write_elem(kind, tmp, slot, globals[[slot]], compression = compression, zarr_format = if (kind == "zarr") zarr_format)
  }
  for (attr in c("obs", "var")) {
    store_write_maps(kind, tmp, attr, maps[[attr]], zarr_format)
  }
  if (kind == "h5mu") {
    rhdf5::h5createGroup(tmp, "mod")
  } else {
    Rarr::write_zarr_group(tmp, "mod", zarr_version = zarr_format)
  }
  store_write_attrs(kind, tmp, "/mod", list(`mod-order` = I(mod_names)), zarr_format)
  for (name in mod_names) {
    store_write_anndata(kind, tmp, name, mods[[name]], compression = compression, zarr_format = if (kind == "zarr") zarr_format)
  }
  store_finalize(kind, tmp)

  if (exists) {
    if (kind == "h5mu" && !dir.exists(path)) {
      # rename(2) replaces an existing file atomically
      ok <- file.rename(tmp, path)
    } else {
      # directories cannot be replaced atomically: remove, then move
      unlink(path, recursive = TRUE, force = TRUE)
      ok <- file.rename(tmp, path)
    }
  } else {
    ok <- file.rename(tmp, path)
  }
  if (!ok) {
    mudatar_abort("Could not move the written store to {.file {path}}", call = call)
  }
  done <- TRUE
  invisible(path)
}
