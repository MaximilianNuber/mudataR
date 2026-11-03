# Working with released anndataR
#
# mudataR needs three things released anndataR does not offer:
#
#   1. reading and writing single elements (MuData's global slots),
#   2. AnnData objects that are not the root of their file (`/mod/<name>`),
#   3. Awkward Arrays (scirpy).
#
# All three are solved here, on top of anndataR's public API only:
#
# * HDF5 reads go through a *shim*: a small temporary .h5ad whose entries are
#   HDF5 external links into the H5MU. anndataR reads the shim like any H5AD.
#   Writes never go through a shim (anndataR would replace the links inside
#   the shim and the H5MU would never change); they go to a temporary H5AD
#   whose groups are copied into the H5MU with H5Ocopy. Shim-backed
#   modalities are therefore opened read-only.
# * Zarr modalities are directories, so anndataR opens `<store>/mod/<name>`
#   directly (also for writing). Shims made of symlinks are used for single
#   elements and for modalities that contain Awkward Arrays.
# * Awkward Arrays are never shown to anndataR: shims leave them out, mudataR
#   reads and writes them itself (awkward-store.R) and puts them into
#   in-memory AnnData objects as AwkwardFrame (aligned rows anndataR subsets
#   with its standard data-frame rules). Backed modalities expose them
#   through get_awkward().

.mudatar_env <- new.env(parent = emptyenv())

# Slots whose direct children may be Awkward Arrays
.awkward_slots <- c("obsm", "varm", "uns")

# anndataR helpers ---------------------------------------------------------------

#' Run anndataR, turning its "could not read/write element" warnings into
#' errors (they would otherwise silently drop data) and muffling the generic
#' read/write-mode caution that does not apply to mudataR's own files
#' @noRd
anndatar_strict <- function(expr, call = rlang::caller_env()) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("non-empty file is opened in read/write mode", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
      if (grepl("Could not write element|Error reading element|Error writing element", msg)) {
        mudatar_abort(
          c("anndataR could not process an element", "x" = msg),
          class = "mudatar_anndatar_error",
          call = call
        )
      }
    }
  )
}

anndatar_read <- function(kind, path, as = "InMemoryAnnData", mode = "r", backed = FALSE) {
  anndatar_strict(
    if (kind == "h5mu") {
      anndataR::read_h5ad(path, as = as, mode = mode, backed = backed)
    } else {
      anndataR::read_zarr(path, as = as, mode = mode)
    }
  )
}

anndatar_write <- function(kind, adata, path, compression = NULL, zarr_format = NULL) {
  args <- list(adata, path, mode = "w")
  if (!is.null(compression)) args$compression <- compression
  if (kind == "zarr") args$zarr_format <- zarr_format
  anndatar_strict(do.call(if (kind == "h5mu") anndataR::write_h5ad else anndataR::write_zarr, args))
}

# Shims --------------------------------------------------------------------------

#' A fresh, empty AnnData store to put links into
#' @noRd
shim_template <- function(kind, zarr_format = 3L) {
  key <- paste0("template_", kind, if (kind == "zarr") zarr_format)
  tpl <- .mudatar_env[[key]]
  if (is.null(tpl) || !file.exists(tpl)) {
    tpl <- tempfile("mudatar-template-", fileext = if (kind == "h5mu") ".h5ad" else ".zarr")
    anndatar_write(kind, anndataR::AnnData(), tpl, zarr_format = zarr_format)
    .mudatar_env[[key]] <- tpl
  }
  out <- tempfile("mudatar-shim-", fileext = if (kind == "h5mu") ".h5ad" else ".zarr")
  if (kind == "h5mu") {
    file.copy(tpl, out)
  } else {
    dir.create(out)
    file.copy(list.files(tpl, all.files = TRUE, full.names = TRUE, no.. = TRUE), out, recursive = TRUE)
  }
  out
}

#' Make `shim:<dest>` show `path:<name>`, leaving out `exclude` children
#'
#' With no exclusions this is a single link; otherwise `dest` is a real group
#' with the source's attributes and one link per remaining child.
#' @noRd
shim_link <- function(kind, shim, dest, path, name, exclude = character()) {
  name <- sub("^/", "", name)
  if (kind == "h5mu") {
    fid <- h5_open(shim, readonly = FALSE)
    on.exit(rhdf5::H5Fclose(fid), add = TRUE)
    if (rhdf5::H5Lexists(fid, dest)) rhdf5::H5Ldelete(fid, dest)
    if (!length(exclude)) {
      rhdf5::H5Lcreate_external(path, name, fid, dest)
      return(invisible(NULL))
    }
    gid <- rhdf5::H5Gcreate(fid, dest)
    on.exit(rhdf5::H5Gclose(gid), add = TRUE, after = FALSE)
    for (k in setdiff(store_list_children(kind, path, name), exclude)) {
      rhdf5::H5Lcreate_external(path, paste0(name, "/", k), gid, k)
    }
    rhdf5::H5Gclose(gid)
    rhdf5::H5Fclose(fid)
    on.exit()
    attrs <- store_read_attrs(kind, path, name)
    if (length(attrs)) store_write_attrs(kind, shim, dest, attrs)
  } else {
    target <- file.path(shim, dest)
    unlink(target, recursive = TRUE)
    src <- file.path(path, name)
    if (!length(exclude)) {
      file.symlink(normalizePath(src), target)
      return(invisible(NULL))
    }
    dir.create(target)
    meta <- intersect(list.files(src, all.files = TRUE, no.. = TRUE), c(".zgroup", ".zattrs", "zarr.json"))
    file.copy(file.path(src, meta), target)
    for (k in setdiff(store_list_children(kind, path, name), exclude)) {
      file.symlink(normalizePath(file.path(src, k)), file.path(target, k))
    }
  }
  invisible(NULL)
}

#' Shim exposing the AnnData at `path:<root>` without its Awkward elements
#'
#' @return `list(path = <shim>, awkward = list(obsm = <keys>, ...))`
#' @noRd
shim_anndata <- function(kind, path, root, zarr_format = NULL) {
  root <- sub("^/", "", root)
  awkward <- awkward_keys(kind, path, root)
  if (kind == "h5mu") {
    shim <- tempfile("mudatar-shim-", fileext = ".h5ad")
    fid <- rhdf5::H5Fcreate(shim)
    rhdf5::H5Fclose(fid)
  } else {
    shim <- tempfile("mudatar-shim-", fileext = ".zarr")
    dir.create(shim)
    src <- file.path(path, root)
    meta <- intersect(list.files(src, all.files = TRUE, no.. = TRUE), c(".zgroup", ".zattrs", "zarr.json"))
    file.copy(file.path(src, meta), shim)
  }
  for (k in store_list_children(kind, path, root)) {
    shim_link(kind, shim, k, path, paste0(root, "/", k), exclude = awkward[[k]] %||% character())
  }
  if (kind == "h5mu") {
    attrs <- store_read_attrs(kind, path, if (nzchar(root)) root else "/")
    if (length(attrs)) store_write_attrs(kind, shim, "/", attrs)
  }
  list(path = shim, awkward = awkward)
}

#' Awkward children per slot of the AnnData at `path:<root>`
#' @noRd
awkward_keys <- function(kind, path, root) {
  root <- sub("^/", "", root)
  kids <- store_list_children(kind, path, root)
  out <- lapply(intersect(.awkward_slots, kids), function(slot) {
    store_awkward_children(kind, path, paste0(root, "/", slot))
  })
  names(out) <- intersect(.awkward_slots, kids)
  Filter(length, out)
}

# Single elements -----------------------------------------------------------------

#' Read one element of a store with anndataR's decoders
#'
#' The element is linked into an empty AnnData as `uns/value`; Awkward
#' children of a mapping are read by mudataR.
#' @noRd
store_read_elem <- function(kind, path, name, lazy = FALSE) {
  name <- sub("^/", "", name)
  if (!store_exists(kind, path, name)) {
    return(NULL)
  }
  type <- store_encoding_type(kind, path, name)
  if (identical(type, "awkward-array")) {
    return(store_read_awkward(kind, path, name, lazy = lazy))
  }
  awk <- if (identical(type, "dict")) store_awkward_children(kind, path, name) else character()
  shim <- shim_template(kind, if (kind == "zarr") store_zarr_format(path) else 3L)
  on.exit(unlink(shim, recursive = TRUE), add = TRUE)
  shim_link(kind, shim, "uns/value", path, name, exclude = awk)
  value <- anndatar_read(kind, shim)$uns$value
  if (length(awk)) {
    value <- value %||% list()
    for (k in awk) value[[k]] <- store_read_awkward(kind, path, paste0(name, "/", k), lazy = lazy)
    value <- value[intersect(store_list_children(kind, path, name), names(value))]
  }
  value
}

#' Write one element of a store with anndataR's encoders
#'
#' The value is written as `uns/value` of a temporary AnnData and moved into
#' place; Awkward values (also as children of a mapping) are written by
#' mudataR. Replaces an existing element.
#' @noRd
store_write_elem <- function(kind, path, name, value, compression = NULL, zarr_format = NULL) {
  name <- sub("^/", "", name)
  if (is_awkward_value(value)) {
    return(store_write_awkward(kind, path, name, as_AwkwardArrowArray(value), zarr_format, compression))
  }
  awk <- list()
  if (is.list(value) && !is.data.frame(value)) {
    is_awk <- vapply(value, is_awkward_value, logical(1))
    awk <- value[is_awk]
    value <- value[!is_awk]
  }
  # validate Awkward values before anything on disk changes
  for (k in names(awk)) awkward_prepare_buffers(as_AwkwardArrowArray(awk[[k]]), paste0(path, "/", name, "/", k), kind)

  if (kind == "h5mu") {
    tmp <- tempfile("mudatar-elem-", fileext = ".h5ad")
    on.exit(unlink(tmp), add = TRUE)
    anndatar_write(kind, anndataR::AnnData(uns = list(value = value)), tmp, compression = compression)
    store_delete(kind, path, name)
    src <- h5_open(tmp)
    dst <- h5_open(path, readonly = FALSE)
    rhdf5::H5Ocopy(src, "uns/value", dst, name)
    rhdf5::H5Fclose(src)
    rhdf5::H5Fclose(dst)
  } else {
    zarr_format <- zarr_format %||% store_zarr_format(path)
    # next to the store, so the result can be renamed into place
    tmp <- tempfile(".mudatar-elem-", tmpdir = dirname(path), fileext = ".zarr")
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    anndatar_write(
      kind,
      anndataR::AnnData(uns = list(value = value)),
      tmp,
      compression = compression,
      zarr_format = zarr_format
    )
    store_delete(kind, path, name)
    if (!file.rename(file.path(tmp, "uns", "value"), file.path(path, name))) {
      mudatar_abort("Could not move element {.val {name}} into {.file {path}}")
    }
  }
  for (k in names(awk)) {
    store_write_awkward(kind, path, paste0(name, "/", k), as_AwkwardArrowArray(awk[[k]]), zarr_format, compression)
  }
  invisible(NULL)
}

# Modalities ---------------------------------------------------------------------

#' Open the modality at `/mod/<name>` as an anndataR object
#'
#' HDF5 modalities are read through a shim and therefore read-only. Zarr
#' modalities are opened in place (honouring `mode`) unless they contain
#' Awkward Arrays, which need a read-only shim. Every returned object is
#' registered so mudataR can find its Awkward elements again.
#' @noRd
store_open_anndata <- function(kind, path, name, mode = "r", backed = FALSE) {
  root <- paste0("mod/", name)
  awkward <- awkward_keys(kind, path, root)
  if (kind == "zarr" && !length(awkward)) {
    adata <- anndatar_read(kind, file.path(path, root), as = "ZarrAnnData", mode = mode)
  } else {
    shim <- shim_anndata(kind, path, root)
    adata <- anndatar_read(
      kind,
      shim$path,
      as = if (kind == "h5mu") "HDF5AnnData" else "ZarrAnnData",
      mode = "r",
      backed = backed
    )
  }
  register_modality(adata, kind, path, name, lazy = backed)
  adata
}

#' Read the modality at `/mod/<name>` (or the AnnData at `root`) into an
#' `InMemoryAnnData`, with its Awkward elements as AwkwardFrame /
#' AwkwardArrowArray
#' @noRd
store_read_anndata <- function(kind, path, name, root = paste0("mod/", name)) {
  shim <- shim_anndata(kind, path, root)
  on.exit(unlink(shim$path, recursive = TRUE), add = TRUE)
  adata <- anndatar_read(kind, shim$path)
  if (length(shim$awkward)) {
    add_awkward(adata, read_awkward_elements(kind, path, root, shim$awkward))
  }
  adata
}

read_awkward_elements <- function(kind, path, root, keys, lazy = FALSE) {
  out <- lapply(names(keys), function(slot) {
    vals <- lapply(keys[[slot]], function(k) {
      store_read_awkward(kind, path, paste0(root, "/", slot, "/", k), lazy = lazy)
    })
    names(vals) <- keys[[slot]]
    vals
  })
  names(out) <- names(keys)
  out
}

#' Write an AnnData to `/mod/<name>`, replacing an existing modality
#' @noRd
store_write_anndata <- function(kind, path, name, adata, compression = NULL, zarr_format = NULL) {
  parts <- split_awkward(adata)
  root <- paste0("mod/", name)
  for (slot in names(parts$awkward)) {
    for (k in names(parts$awkward[[slot]])) {
      awkward_prepare_buffers(parts$awkward[[slot]][[k]], paste0(path, "/", root, "/", slot, "/", k), kind)
    }
  }
  if (kind == "h5mu") {
    tmp <- tempfile("mudatar-mod-", fileext = ".h5ad")
    on.exit(unlink(tmp), add = TRUE)
    anndatar_write(kind, parts$adata, tmp, compression = compression)
    store_delete(kind, path, root)
    src <- h5_open(tmp)
    dst <- h5_open(path, readonly = FALSE)
    rhdf5::H5Ocopy(src, "/", dst, root)
    rhdf5::H5Fclose(src)
    rhdf5::H5Fclose(dst)
  } else {
    zarr_format <- zarr_format %||% store_zarr_format(path)
    store_delete(kind, path, root)
    anndatar_write(kind, parts$adata, file.path(path, root), compression = compression, zarr_format = zarr_format)
  }
  for (slot in names(parts$awkward)) {
    for (k in names(parts$awkward[[slot]])) {
      store_write_awkward(
        kind,
        path,
        paste0(root, "/", slot, "/", k),
        parts$awkward[[slot]][[k]],
        zarr_format = zarr_format,
        compression = compression
      )
    }
  }
  invisible(NULL)
}

# In-memory AnnData with Awkward elements ------------------------------------------

#' Split an AnnData into an `InMemoryAnnData` without Awkward values (safe
#' for anndataR's writers) and the Awkward values per slot
#' @noRd
split_awkward <- function(adata) {
  adata <- as_in_memory_anndata(adata)
  slots <- lapply(stats::setNames(nm = .awkward_slots), function(s) adata[[s]] %||% list())
  awkward <- lapply(slots, function(v) lapply(Filter(is_awkward_value, v), as_AwkwardArrowArray))
  awkward <- Filter(length, awkward)
  if (!length(awkward)) {
    return(list(adata = adata, awkward = list()))
  }
  plain <- lapply(slots, function(v) Filter(Negate(is_awkward_value), v))
  stripped <- anndataR::AnnData(
    X = adata$X,
    obs = adata$obs,
    var = adata$var,
    layers = adata$layers,
    obsm = plain$obsm,
    varm = plain$varm,
    obsp = adata$obsp,
    varp = adata$varp,
    uns = plain$uns
  )
  list(adata = stripped, awkward = awkward)
}

#' Add Awkward values to an `InMemoryAnnData` (aligned ones as AwkwardFrame)
#' @noRd
add_awkward <- function(adata, awkward) {
  for (slot in names(awkward)) {
    vals <- awkward[[slot]]
    if (!length(vals)) next
    if (slot %in% c("obsm", "varm")) vals <- lapply(vals, as_AwkwardFrame)
    current <- adata[[slot]] %||% list()
    current[names(vals)] <- vals
    adata[[slot]] <- current[sort(names(current))]
  }
  invisible(adata)
}

#' Any anndataR object as an `InMemoryAnnData` that keeps its Awkward
#' elements
#'
#' `InMemoryAnnData` objects are returned as they are. Modalities opened by
#' mudataR are re-read from their store; other objects are copied slot by
#' slot so that Awkward values (e.g. from a `ReticulateAnnData`) become
#' AwkwardFrame instead of failing anndataR's validation.
#' @noRd
as_in_memory_anndata <- function(adata) {
  if (inherits(adata, "InMemoryAnnData")) {
    return(adata)
  }
  info <- modality_info(adata)
  if (!is.null(info)) {
    return(store_read_anndata(info$kind, info$path, info$name))
  }
  if (inherits(adata, "AnnDataView")) {
    return(adata$as_InMemoryAnnData())
  }
  if (inherits(adata, "ReticulateAnnData")) {
    return(py_anndata_to_r(adata$py_anndata()))
  }
  slots <- lapply(stats::setNames(nm = .awkward_slots), function(s) adata[[s]] %||% list())
  if (!any(vapply(slots, function(v) any(vapply(v, is_awkward_value, logical(1))), logical(1)))) {
    return(adata$as_InMemoryAnnData())
  }
  out <- anndataR::AnnData(
    X = adata$X,
    obs = adata$obs,
    var = adata$var,
    layers = adata$layers,
    obsm = Filter(Negate(is_awkward_value), slots$obsm),
    varm = Filter(Negate(is_awkward_value), slots$varm),
    obsp = adata$obsp,
    varp = adata$varp,
    uns = Filter(Negate(is_awkward_value), slots$uns)
  )
  add_awkward(out, lapply(slots, function(v) lapply(Filter(is_awkward_value, v), as_AwkwardArrowArray)))
}

#' A Python anndata.AnnData as an `InMemoryAnnData`: written by Python to a
#' temporary H5AD (Awkward Arrays included) and read natively
#' @noRd
py_anndata_to_r <- function(py_adata) {
  tmp <- tempfile("mudatar-py-", fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  py_adata$write_h5ad(tmp)
  store_read_anndata("h5mu", tmp, root = "")
}

# Registry of modalities opened from a store ---------------------------------------

#' Remember where a modality object came from (weakly: the registry never
#' keeps an object alive)
#' @noRd
register_modality <- function(adata, kind, path, name, lazy = FALSE) {
  reg <- Filter(function(w) !is.null(rlang::wref_key(w)), .mudatar_env$modalities %||% list())
  info <- list(kind = kind, path = path, name = name, lazy = lazy)
  .mudatar_env$modalities <- c(reg, list(rlang::new_weakref(adata, value = info)))
  invisible(adata)
}

modality_info <- function(adata) {
  if (!is.environment(adata)) {
    return(NULL)
  }
  for (w in .mudatar_env$modalities %||% list()) {
    if (identical(rlang::wref_key(w), adata)) {
      return(rlang::wref_value(w))
    }
  }
  NULL
}

#' Awkward Arrays of an anndataR object
#'
#' Returns the Awkward Arrays (e.g. scirpy's `obsm[["airr"]]`) of a modality
#' as [AwkwardArrowArray] objects. This is the way to reach them for backed
#' modalities (`HDF5MuData`/`ZarrMuData`), which anndataR opens without
#' their Awkward elements. For in-memory AnnData objects the
#' [AwkwardFrame]s in `obsm`/`varm` (and Awkward values in `uns`) are
#' returned as arrays, subset to the object's current rows.
#'
#' @param adata An anndataR AnnData object, typically a MuData modality
#' @param slot `"obsm"` (default), `"varm"` or `"uns"`
#'
#' @return A named list of [AwkwardArrowArray] objects (possibly empty)
#' @examples
#' path <- system.file("extdata", "gex_airr.h5mu", package = "mudataR")
#' if (nzchar(path) && requireNamespace("rhdf5", quietly = TRUE) &&
#'     requireNamespace("nanoarrow", quietly = TRUE)) {
#'   mdata <- read_h5mu(path, as = "HDF5MuData")
#'   get_awkward(mdata[["airr"]])$airr
#' }
#' @export
get_awkward <- function(adata, slot = c("obsm", "varm", "uns")) {
  slot <- match.arg(slot)
  if (!inherits(adata, "AbstractAnnData")) {
    mudatar_abort("{.arg adata} must be an anndataR AnnData object, not {.cls {class(adata)}}")
  }
  info <- modality_info(adata)
  if (!is.null(info)) {
    root <- paste0("mod/", info$name)
    keys <- awkward_keys(info$kind, info$path, root)[slot]
    keys <- Filter(length, stats::setNames(keys, slot))
    if (!length(keys)) {
      return(list())
    }
    return(read_awkward_elements(info$kind, info$path, root, keys, lazy = info$lazy)[[slot]])
  }
  lapply(Filter(is_awkward_value, adata[[slot]] %||% list()), as_AwkwardArrowArray)
}
