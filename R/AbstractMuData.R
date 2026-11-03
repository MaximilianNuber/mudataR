#' @title AbstractMuData
#'
#' @description
#' Abstract R6 class defining the interface shared by all MuData backends
#' ([InMemoryMuData], [HDF5MuData], [ZarrMuData]). It mirrors
#' `anndataR::AbstractAnnData`: slots are active bindings, operations are
#' methods, and standard R generics (`dim()`, `[[`, `print()`, ...) are
#' provided as S3 methods.
#'
#' Global observation and variable names are stored separately from the
#' `obs`/`var` data frames, as in anndataR and Python.
#'
#' @section Slots (active bindings):
#' * `mod`: a [ModalityMapping] of the modalities, each a genuine anndataR
#'   AnnData object of the matching backend
#' * `obs`, `var`: global annotation data frames
#' * `obsm`, `varm`, `obsp`, `varp`, `uns`: aligned, pairwise and unstructured
#'   global annotation
#' * `obsmap`, `varmap`: named lists of integer maps (0 = absent, otherwise a
#'   one-based position in the modality)
#' * `obs_names`, `var_names`: global names
#' * `axis`: `0` (shared observations), `1` (shared variables) or `-1` (both)
#' * `filename`, `isbacked`, `is_view`: read-only
#'
#' @return An `AbstractMuData` subclass object
#' @family MuData classes
#' @export
AbstractMuData <- R6::R6Class(
  "AbstractMuData", # nolint
  cloneable = FALSE,
  active = list(
    #' @field mod The modalities, see [ModalityMapping]
    mod = function(value) {
      if (!missing(value)) {
        # `mdata$mod[["x"]] <- adata` assigns the (already updated) mapping
        # back; anything else would bypass validation
        if (inherits(value, "ModalityMapping") && identical(mapping_owner(value), self)) {
          return(invisible(self))
        }
        mudatar_abort("Assign individual modalities with {.code mdata[[\"name\"]] <- adata}")
      }
      new_modality_mapping(self)
    },
    #' @field obs Global observation annotation
    obs = function(value) private$.abstract("obs"),
    #' @field var Global variable annotation
    var = function(value) private$.abstract("var"),
    #' @field obsm Global aligned observation annotation
    obsm = function(value) private$.abstract("obsm"),
    #' @field varm Global aligned variable annotation
    varm = function(value) private$.abstract("varm"),
    #' @field obsp Global pairwise observation annotation
    obsp = function(value) private$.abstract("obsp"),
    #' @field varp Global pairwise variable annotation
    varp = function(value) private$.abstract("varp"),
    #' @field uns Global unstructured annotation
    uns = function(value) private$.abstract("uns"),
    #' @field obsmap Observation maps, one integer vector per modality
    obsmap = function(value) private$.abstract("obsmap"),
    #' @field varmap Variable maps, one integer vector per modality
    varmap = function(value) private$.abstract("varmap"),
    #' @field obs_names Global observation names
    obs_names = function(value) private$.abstract("obs_names"),
    #' @field var_names Global variable names
    var_names = function(value) private$.abstract("var_names"),
    #' @field axis The shared axis (read-only)
    axis = function(value) private$.abstract("axis"),
    #' @field filename The backing file, or `NULL` (read-only)
    filename = function(value) {
      if (!missing(value)) {
        mudatar_abort_unsupported(
          "Re-binding {.field filename}",
          "i" = "Use {.fn write_h5mu} and {.fn read_h5mu} instead."
        )
      }
      NULL
    },
    #' @field isbacked Whether the object is backed by a file (read-only)
    isbacked = function(value) {
      if (!missing(value)) mudatar_abort_readonly("field")
      FALSE
    },
    #' @field is_view Whether the object is a view (read-only)
    is_view = function(value) {
      if (!missing(value)) mudatar_abort_readonly("field")
      FALSE
    }
  ),
  public = list(
    #' @description Names of the modalities, in order
    mod_names = function() private$.abstract("mod_names"),

    #' @description Get a modality (used by [ModalityMapping])
    #' @param name Modality name
    get_mod = function(name) private$.abstract("get_mod"),

    #' @description Add, replace (or with `value = NULL` remove) a modality
    #'   and rebuild the global indices, as Python's `mdata[name] = adata`
    #' @param name Modality name
    #' @param value An anndataR AnnData object or `NULL`
    set_mod = function(name, value) private$.abstract("set_mod"),

    #' @description Dimensions `c(n_obs, n_vars)`
    shape = function() {
      c(self$n_obs(), self$n_vars())
    },
    #' @description Number of global observations
    n_obs = function() {
      length(self$obs_names)
    },
    #' @description Number of global variables
    n_vars = function() {
      length(self$var_names)
    },
    #' @description Number of modalities
    n_mod = function() {
      length(self$mod_names())
    },
    #' @description Global `obs` column names
    obs_keys = function() {
      colnames(self$obs)
    },
    #' @description Global `var` column names
    var_keys = function() {
      colnames(self$var)
    },
    #' @description Global `obsm` keys
    obsm_keys = function() {
      names(self$obsm)
    },
    #' @description Global `varm` keys
    varm_keys = function() {
      names(self$varm)
    },
    #' @description Global `obsp` keys
    obsp_keys = function() {
      names(self$obsp)
    },
    #' @description Global `varp` keys
    varp_keys = function() {
      names(self$varp)
    },
    #' @description Global `uns` keys
    uns_keys = function() {
      names(self$uns)
    },

    #' @description Rebuild global observation names and `obsmap` from the
    #'   modalities (port of Python `MuData.update_obs()`). Global `obs`
    #'   columns, `obsm` and `obsp` are reindexed; modality columns are not
    #'   pulled.
    update_obs = function() {
      private$.update_attr("obs")
      invisible(self)
    },
    #' @description Rebuild global variable names and `varmap`
    update_var = function() {
      private$.update_attr("var")
      invisible(self)
    },
    #' @description Rebuild both observation and variable indices
    update = function() {
      self$update_obs()
      self$update_var()
      invisible(self)
    },

    #' @description Convert to an [InMemoryMuData], materialising every
    #'   modality as an `anndataR` `InMemoryAnnData`
    as_InMemoryMuData = function() {
      # nolint start: object_name_linter
      mods <- lapply(as.list(self$mod), as_in_memory_anndata)
      # nolint end: object_name_linter
      InMemoryMuData$new(
        mod = mods,
        obs = self$obs,
        var = self$var,
        obsm = self$obsm,
        varm = self$varm,
        obsp = self$obsp,
        varp = self$varp,
        uns = self$uns,
        obsmap = self$obsmap,
        varmap = self$varmap,
        obs_names = self$obs_names,
        var_names = self$var_names,
        axis = self$axis
      )
    },

    #' @description Convert to a `SingleCellExperiment`: one modality (`main`,
    #'   by default `"rna"`) becomes the main experiment, every other
    #'   modality an `altExp()` of the same name. Each modality is converted
    #'   with `anndataR::as_SingleCellExperiment()` using its entry in
    #'   `mod_args`.
    #'
    #'   Global `obs` columns are added to the main `colData` (prefixed with
    #'   `mudata:` on a name clash), global `obsm` matrices become
    #'   `reducedDims` of the main experiment, and global `uns`/`obsm` are kept
    #'   in `metadata(sce)$mudata`. Awkward Arrays (e.g. scirpy AIRR data) are
    #'   kept in `metadata(<experiment>)$awkward_obsm`.
    #' @param main Name of the modality used as the main experiment. Defaults
    #'   to `"rna"` if present, otherwise the first modality.
    #' @param mod_args Named list with one list of
    #'   `anndataR::as_SingleCellExperiment()` arguments per modality, e.g.
    #'   `list(rna = list(x_mapping = "counts"), prot = list(x_mapping = "X"))`.
    #'   Modalities not listed use anndataR's defaults.
    #' @param join Which cells to keep. `"intersect"` (default): cells present
    #'   in every modality. `"main"`: every cell of `main`; altExps are padded
    #'   with `NA` for cells they lack and get a logical `present` column.
    #'   Either way, cells keep the order they have in `main`.
    as_SingleCellExperiment = function(main = NULL, mod_args = list(), join = c("intersect", "main")) {
      # nolint start: object_name_linter
      mudata_to_sce(self, main = main, mod_args = mod_args, join = match.arg(join))
      # nolint end: object_name_linter
    },

    #' @description Convert to a `MultiAssayExperiment`, MuData's closest
    #'   Bioconductor counterpart: every modality becomes an experiment
    #'   (converted with `anndataR::as_SingleCellExperiment()` and its
    #'   `mod_args` entry), global `obs` becomes `colData`, and `obsmap`
    #'   becomes the `sampleMap`. Nothing is dropped or padded. Global
    #'   `uns`, `obsm` and `obsp` are kept in `metadata(mae)$mudata`.
    #' @param mod_args See `as_SingleCellExperiment()`
    as_MultiAssayExperiment = function(mod_args = list()) {
      # nolint start: object_name_linter
      mudata_to_mae(self, mod_args = mod_args)
      # nolint end: object_name_linter
    },

    #' @description Convert to a [ReticulateMuData] (a Python `mudata.MuData`
    #'   whose modalities are `anndataR::ReticulateAnnData`); needs Python
    #'   with `mudata`, see [as_ReticulateMuData()]
    as_ReticulateMuData = function() {
      # nolint start: object_name_linter
      mudata_to_reticulate(self)
      # nolint end: object_name_linter
    },

    #' @description An independent in-memory copy
    copy = function() {
      self$as_InMemoryMuData()
    },

    #' @description Write to an H5MU file, see [write_h5mu()]
    #' @param path Output path
    #' @param ... Passed to [write_h5mu()]
    write_h5mu = function(path, ...) {
      write_h5mu(self, path, ...)
    },

    #' @description Write to a MuData Zarr store, see [write_mudata_zarr()]
    #' @param path Output path
    #' @param ... Passed to [write_mudata_zarr()]
    write_zarr = function(path, ...) {
      write_mudata_zarr(self, path, ...)
    },

    #' @description Format a summary
    #' @param ... Unused
    format = function(...) {
      shape <- self$shape()
      out <- paste0(
        "MuData object with n_obs \u00d7 n_vars = ",
        shape[[1]],
        " \u00d7 ",
        shape[[2]],
        " (axis ",
        self$axis,
        ")"
      )
      backing <- private$.backing_label()
      if (!is.null(backing)) {
        out <- c(out, paste0("  backed at ", backing))
      }
      for (slot in c("obs", "var")) {
        keys <- .subset2(self, paste0(slot, "_keys"))()
        if (length(keys)) out <- c(out, paste0("  ", slot, ": ", paste0("'", keys, "'", collapse = ", ")))
      }
      for (slot in c("obsm", "varm", "obsp", "varp", "uns")) {
        keys <- .subset2(self, paste0(slot, "_keys"))()
        if (length(keys)) out <- c(out, paste0("  ", slot, ": ", paste0("'", keys, "'", collapse = ", ")))
      }
      nms <- self$mod_names()
      out <- c(out, paste0("  ", length(nms), " modalit", if (length(nms) == 1) "y" else "ies"))
      for (n in nms) {
        a <- self$get_mod(n)
        out <- c(out, paste0("    ", n, ": ", a$n_obs(), " x ", a$n_vars()))
        for (slot in c("obs", "var", "obsm", "varm", "layers")) {
          # some backends cannot list every slot (e.g. anndataR's
          # ReticulateAnnData layers with Python anndata >= 0.13)
          keys <- tryCatch(a[[paste0(slot, "_keys")]](), error = function(e) NULL)
          if (length(keys)) {
            out <- c(out, paste0("      ", slot, ": ", paste0("'", keys, "'", collapse = ", ")))
          }
        }
      }
      out
    },
    #' @description Print a summary
    #' @param ... Unused
    print = function(...) {
      cat(self$format(), sep = "\n")
      invisible(self)
    }
  ),
  private = list(
    .abstract = function(what) {
      mudatar_abort(
        "{.fn {what}} is not implemented for {.cls {class(self)[[1]]}}",
        class = "mudatar_not_implemented"
      )
    },

    .backing_label = function() NULL,

    # Atomically replace the global names, maps and aligned elements of one
    # attribute. Backends implement this; update logic stays here.
    .replace_attr = function(attr, names, df, maps, m, p) {
      private$.abstract("replace_attr")
    },

    .update_attr = function(attr) {
      names_field <- paste0(attr, "_names")
      mods <- as.list(self$mod)
      mod_names <- lapply(mods, function(a) a[[names_field]])
      prev_names <- .subset2(self, names_field) %||% character()
      prev_df <- .subset2(self, attr)
      res <- compute_global_attr(mod_names, attr, self$axis, prev_names)

      n <- length(res$names)
      df <- if (is.null(prev_df) || ncol(prev_df) == 0 || is.null(res$index_order)) {
        as.data.frame(matrix(nrow = n, ncol = 0))
      } else {
        reindex_rows(prev_df, res$index_order)
      }
      rownames(df) <- NULL

      m <- .subset2(self, paste0(attr, "m")) %||% list()
      m <- m[setdiff(names(m), names(mods))]
      m <- lapply(m, reindex_rows, index_order = res$index_order)
      for (mod in names(res$maps)) {
        m[[mod]] <- res$maps[[mod]] > 0L
      }
      p <- .subset2(self, paste0(attr, "p")) %||% list()
      p <- lapply(p, reindex_pairwise, index_order = res$index_order)

      private$.replace_attr(attr, res$names, df, res$maps, m, p)
    },

    # Strip/validate rownames of a global data frame
    .prepare_df = function(df, attr, names) {
      if (is.null(df)) {
        df <- as.data.frame(matrix(nrow = length(names), ncol = 0))
      }
      if (!is.data.frame(df)) {
        mudatar_abort("{.field {attr}} must be a data frame")
      }
      if (nrow(df) != length(names)) {
        mudatar_abort("{.field {attr}} must have {length(names)} rows, not {nrow(df)}")
      }
      if (.row_names_info(df) > 0 && !identical(rownames(df), names)) {
        mudatar_abort("The row names of {.field {attr}} do not match {.field {attr}_names}")
      }
      rownames(df) <- NULL
      df
    },

    .add_df_names = function(df, names) {
      if (!is.null(df) && !anyDuplicated(names) && nrow(df) == length(names)) {
        rownames(df) <- names
      }
      df
    }
  )
)
