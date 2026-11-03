# Shared implementation of store-backed MuData objects. HDF5MuData and
# ZarrMuData only differ in the store kind passed to the helpers in store.R.
#
# Nothing is read at construction beyond the root attributes and the list of
# modalities. Global slots are read on access with anndataR's decoders
# (store_read_elem()); modalities are opened on access as anndataR objects
# (store_open_anndata()): read-only HDF5AnnData objects reading the H5MU
# through a shim, and ZarrAnnData objects opened at `<store>/mod/<name>`.
#' Shared base class of HDF5MuData and ZarrMuData
#' @noRd
BackedMuData <- R6::R6Class(
  "BackedMuData", # nolint
  inherit = AbstractMuData,
  cloneable = FALSE,
  private = list(
    .kind = NULL,
    .path = NULL,
    .mode = "r",
    .backed = FALSE,
    .is_open = TRUE,
    .selected = NULL,
    .zarr_format = NULL,

    .check_open = function() {
      if (!isTRUE(private$.is_open)) {
        mudatar_abort_closed(private$.path)
      }
    },
    .check_writeable = function() {
      private$.check_open()
      if (private$.mode == "r") {
        mudatar_abort_readonly(paste0(class(self)[[1]], " opened with mode = \"r\""))
      }
      if (!is.null(private$.selected)) {
        mudatar_abort_unsupported("Modifying a MuData opened with a subset of {.arg modalities}")
      }
    },
    .read = function(name, backed = FALSE) {
      private$.check_open()
      store_read_elem(private$.kind, private$.path, name, lazy = backed)
    },
    .write = function(value, name) {
      private$.check_writeable()
      store_write_elem(private$.kind, private$.path, name, value, zarr_format = private$.zarr_format)
    },
    .read_mapping = function(name) {
      value <- private$.read(name, backed = private$.backed) %||% list()
      lapply(value, function(x) {
        if (is.null(dim(x)) || inherits(x, "AwkwardArrowArray")) x else drop_1d(x)
      })
    },
    .names = function(attr) {
      private$.check_open()
      df <- private$.read(attr)
      rownames(df) %||% character()
    },
    .backing_label = function() {
      private$.path
    },
    .replace_attr = function(attr, names, df, maps, m, p) {
      private$.check_writeable()
      rownames(df) <- if (anyDuplicated(names)) NULL else names
      if (anyDuplicated(names)) {
        mudatar_abort_unsupported("Writing duplicated global {attr}_names to a backed MuData")
      }
      # write everything from validated in-memory values
      private$.write(df, attr)
      store_write_maps(private$.kind, private$.path, attr, maps, private$.zarr_format)
      private$.write(m, paste0(attr, "m"))
      private$.write(p, paste0(attr, "p"))
      invisible(NULL)
    }
  ),
  active = list(
    #' @field obs See [AbstractMuData]
    obs = function(value) {
      if (missing(value)) {
        return(private$.read("obs"))
      }
      private$.check_writeable()
      names <- self$obs_names
      df <- private$.prepare_df(value, "obs", names)
      rownames(df) <- names
      private$.write(df, "obs")
    },
    #' @field var See [AbstractMuData]
    var = function(value) {
      if (missing(value)) {
        return(private$.read("var"))
      }
      private$.check_writeable()
      names <- self$var_names
      df <- private$.prepare_df(value, "var", names)
      rownames(df) <- names
      private$.write(df, "var")
    },
    #' @field obsm See [AbstractMuData]
    obsm = function(value) {
      if (missing(value)) return(private$.read_mapping("obsm"))
      private$.write(check_aligned_mapping(value, "obsm", self$n_obs()), "obsm")
    },
    #' @field varm See [AbstractMuData]
    varm = function(value) {
      if (missing(value)) return(private$.read_mapping("varm"))
      private$.write(check_aligned_mapping(value, "varm", self$n_vars()), "varm")
    },
    #' @field obsp See [AbstractMuData]
    obsp = function(value) {
      if (missing(value)) return(private$.read_mapping("obsp"))
      private$.write(check_aligned_mapping(value, "obsp", self$n_obs(), pairwise = TRUE), "obsp")
    },
    #' @field varp See [AbstractMuData]
    varp = function(value) {
      if (missing(value)) return(private$.read_mapping("varp"))
      private$.write(check_aligned_mapping(value, "varp", self$n_vars(), pairwise = TRUE), "varp")
    },
    #' @field uns See [AbstractMuData]
    uns = function(value) {
      if (missing(value)) return(private$.read("uns", backed = private$.backed) %||% list())
      private$.write(value %||% list(), "uns")
    },
    #' @field obsmap See [AbstractMuData]
    obsmap = function(value) {
      if (missing(value)) {
        private$.check_open()
        maps <- store_read_maps(private$.kind, private$.path, "obs")
        # stored mappings are listed alphabetically; present them in modality order
        return(maps[intersect(self$mod_names(), names(maps))])
      }
      private$.check_writeable()
      sizes <- vapply(as.list(self$mod), function(a) a$n_obs(), integer(1))
      store_write_maps(private$.kind, private$.path, "obs", validate_maps(value, "obs", self$obs_names, sizes), private$.zarr_format)
    },
    #' @field varmap See [AbstractMuData]
    varmap = function(value) {
      if (missing(value)) {
        private$.check_open()
        maps <- store_read_maps(private$.kind, private$.path, "var")
        # stored mappings are listed alphabetically; present them in modality order
        return(maps[intersect(self$mod_names(), names(maps))])
      }
      private$.check_writeable()
      sizes <- vapply(as.list(self$mod), function(a) a$n_vars(), integer(1))
      store_write_maps(private$.kind, private$.path, "var", validate_maps(value, "var", self$var_names, sizes), private$.zarr_format)
    },
    #' @field obs_names See [AbstractMuData] (read-only for backed objects)
    obs_names = function(value) {
      if (!missing(value)) mudatar_abort_unsupported("Renaming observations of a backed MuData")
      private$.names("obs")
    },
    #' @field var_names See [AbstractMuData] (read-only for backed objects)
    var_names = function(value) {
      if (!missing(value)) mudatar_abort_unsupported("Renaming variables of a backed MuData")
      private$.names("var")
    },
    #' @field axis See [AbstractMuData]
    axis = function(value) {
      if (!missing(value)) mudatar_abort_readonly("axis")
      private$.check_open()
      axis <- store_read_attrs(private$.kind, private$.path, "/")[["axis"]]
      # files written before axis was recorded are axis-0 MuData
      if (is.null(axis)) 0L else as.integer(axis)
    },
    #' @field filename See [AbstractMuData]
    filename = function(value) {
      if (!missing(value)) {
        mudatar_abort_unsupported(
          "Re-binding {.field filename}",
          "i" = "Use {.fn write_h5mu} and {.fn read_h5mu} instead."
        )
      }
      private$.path
    },
    #' @field isbacked See [AbstractMuData]
    isbacked = function(value) {
      if (!missing(value)) mudatar_abort_readonly("field")
      TRUE
    }
  ),
  public = list(
    #' @description Open an existing MuData store
    #' @param path Path to the store
    #' @param mode `"r"` (read-only, default) or `"r+"` (read/write)
    #' @param backed Whether arrays are returned as `DelayedArray` and Awkward
    #'   buffers (see [get_awkward()]) loaded lazily (HDF5 only)
    #' @param modalities Optional subset of modality names to expose. Other
    #'   modalities stay untouched on disk.
    #' @param kind Store kind (internal)
    initialize = function(path, mode = c("r", "r+"), backed = FALSE, modalities = NULL, kind) {
      mode <- match.arg(mode)
      store_check_requires(kind)
      validate_mudata_store(kind, path)
      private$.kind <- kind
      private$.path <- normalizePath(path, mustWork = TRUE)
      private$.mode <- mode
      private$.backed <- isTRUE(backed)
      if (kind == "zarr") {
        private$.zarr_format <- store_zarr_format(path)
      }
      if (!is.null(modalities)) {
        available <- store_mod_names(kind, private$.path)
        missing <- setdiff(modalities, available)
        if (length(missing)) {
          mudatar_abort("Modalities {.val {missing}} are not in {.file {path}}; available: {.val {available}}")
        }
        private$.selected <- available[available %in% modalities]
      }
      invisible(self)
    },

    #' @description See [AbstractMuData]
    mod_names = function() {
      private$.check_open()
      private$.selected %||% store_mod_names(private$.kind, private$.path)
    },

    #' @description Get a modality as an anndataR object reading this store.
    #'   H5MU modalities are read-only `HDF5AnnData` objects (modify them with
    #'   `mdata[[name]] <- adata`); Zarr modalities are `ZarrAnnData` objects
    #'   opened with this object's `mode`, read-only if they hold Awkward
    #'   Arrays. Awkward Arrays are available with [get_awkward()].
    #' @param name Modality name
    get_mod = function(name) {
      private$.check_open()
      name <- resolve_mod_name(self$mod_names(), name)
      store_open_anndata(
        private$.kind,
        private$.path,
        name,
        mode = private$.mode,
        backed = private$.backed
      )
    },

    #' @description See [AbstractMuData]
    #' @param name Modality name
    #' @param value An anndataR AnnData object or `NULL`
    set_mod = function(name, value) {
      private$.check_writeable()
      validate_mod_name(name)
      if (is.null(value)) {
        mudatar_abort_unsupported("Removing a modality from a backed MuData")
      }
      check_modalities(stats::setNames(list(value), name))
      # read everything (including Awkward Arrays) before the store changes:
      # the value may be a modality of this very store
      value <- as_in_memory_anndata(value)
      order <- store_mod_names(private$.kind, private$.path)
      store_write_anndata(private$.kind, private$.path, name, value, zarr_format = private$.zarr_format)
      store_write_attrs(private$.kind, private$.path, "/mod", list(`mod-order` = I(unique(c(order, name)))), private$.zarr_format)
      self$update()
      invisible(self)
    },

    #' @description Re-open a closed object
    open = function() {
      if (!file.exists(private$.path)) {
        mudatar_abort("Cannot re-open {.file {private$.path}}: it no longer exists")
      }
      private$.is_open <- TRUE
      invisible(self)
    },

    #' @description Close the object. Its global slots and modalities can no
    #'   longer be accessed through it; modality objects obtained earlier stay
    #'   usable (they hold their own handles) and closing one never affects
    #'   this object or other modalities.
    close = function() {
      private$.is_open <- FALSE
      invisible(self)
    },

    #' @description Whether the object is open
    is_open = function() {
      isTRUE(private$.is_open)
    }
  )
)

#' @title HDF5MuData
#'
#' @description
#' A MuData object backed by an H5MU file. Construction reads only the root
#' attributes; global slots are read on access and modalities are returned as
#' read-only `anndataR` `HDF5AnnData` objects reading `/mod/<name>` of the
#' same file (replace a modality with `mdata[[name]] <- adata`). Their Awkward
#' Arrays are available with [get_awkward()]. Open one with [read_h5mu()].
#'
#' @return An `HDF5MuData` object
#' @family MuData classes
#' @export
HDF5MuData <- R6::R6Class(
  "HDF5MuData", # nolint
  inherit = BackedMuData,
  cloneable = FALSE,
  public = list(
    #' @description Open an H5MU file, see [read_h5mu()]
    #' @param path Path to the `.h5mu` file
    #' @param mode `"r"` or `"r+"`
    #' @param backed See [read_h5mu()]
    #' @param modalities See [read_h5mu()]
    initialize = function(path, mode = c("r", "r+"), backed = FALSE, modalities = NULL) {
      super$initialize(path, mode = mode, backed = backed, modalities = modalities, kind = "h5mu")
    }
  )
)

#' @title ZarrMuData
#'
#' @description
#' A MuData object backed by a MuData Zarr store (v2 or v3). Modalities are
#' returned as `anndataR` `ZarrAnnData` objects opened at `<store>/mod/<name>`
#' (read-only when they hold Awkward Arrays, which [get_awkward()] returns).
#' Open one with [read_mudata_zarr()].
#'
#' @return A `ZarrMuData` object
#' @family MuData classes
#' @export
ZarrMuData <- R6::R6Class(
  "ZarrMuData", # nolint
  inherit = BackedMuData,
  cloneable = FALSE,
  public = list(
    #' @description Open a MuData Zarr store, see [read_mudata_zarr()]
    #' @param path Path to the Zarr store
    #' @param mode `"r"` or `"r+"`
    #' @param modalities See [read_mudata_zarr()]
    initialize = function(path, mode = c("r", "r+"), modalities = NULL) {
      super$initialize(path, mode = mode, backed = FALSE, modalities = modalities, kind = "zarr")
    },

    #' @description The Zarr format (2 or 3) of the store
    zarr_format = function() {
      private$.zarr_format
    }
  )
)

validate_mudata_store <- function(kind, path, call = rlang::caller_env()) {
  if (kind == "h5mu") {
    if (!file.exists(path) || dir.exists(path)) {
      mudatar_abort("H5MU file {.file {path}} does not exist", call = call)
    }
    if (!rhdf5::H5Fis_hdf5(path)) {
      mudatar_abort("{.file {path}} is not an HDF5 file", call = call)
    }
  } else if (!dir.exists(path)) {
    mudatar_abort("Zarr store {.file {path}} does not exist", call = call)
  }
  attrs <- store_read_attrs(kind, path, "/")
  type <- attrs[["encoding-type"]]
  groups <- store_list_groups(kind, path, "/")
  if (!identical(type, .mudata_encoding_type) && !"mod" %in% groups) {
    mudatar_abort(
      c(
        "{.file {path}} is not a MuData store",
        "i" = "Expected {.code encoding-type = \"MuData\"} or a {.field mod} group; got {.val {type %||% 'none'}}."
      ),
      call = call
    )
  }
  version <- attrs[["encoding-version"]]
  if (!is.null(version) && !identical(version, .mudata_encoding_version)) {
    mudatar_abort_unsupported(
      "MuData encoding version {.val {version}}",
      "i" = "Supported: {.val {(.mudata_encoding_version)}}.",
      call = call
    )
  }
  invisible(TRUE)
}

drop_1d <- function(x) {
  if (length(dim(x)) == 1L) as.vector(x) else x
}

check_aligned_mapping <- function(value, attr, n, pairwise = FALSE, call = rlang::caller_env()) {
  value <- value %||% list()
  if (!is.list(value) || is.data.frame(value) || (length(value) && is.null(names(value)))) {
    mudatar_abort("{.field {attr}} must be a named list", call = call)
  }
  for (key in names(value)) {
    x <- value[[key]]
    len <- if (inherits(x, "AwkwardArrowArray")) length(x) else NROW(x)
    if (len != n || (pairwise && NCOL(x) != n)) {
      mudatar_abort("{.field {attr}[[{.val {key}}]]} must be aligned to {n} global names", call = call)
    }
  }
  value
}
