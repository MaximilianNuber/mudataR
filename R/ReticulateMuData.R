#' @title ReticulateMuData
#'
#' @description
#' A MuData object backed by a Python `mudata.MuData` accessed through
#' \pkg{reticulate}. All data stay in Python; slots are converted on access.
#' Modalities are returned as `anndataR::ReticulateAnnData` objects (nested
#' MuData modalities as `ReticulateMuData`), so the shared [AbstractMuData]
#' interface works unchanged.
#'
#' Python is only needed for this backend: native H5MU/Zarr reading and
#' writing never load it. `update()`, `pull_obs()`/`pull_var()`,
#' `push_obs()`/`push_var()` and `copy()` call the Python implementation.
#'
#' Refactored from the `ReticulateMuData` prototype in MofaflexR.
#'
#' @section Python access:
#' * `$py_mudata()` returns the underlying Python object (`convert = FALSE`)
#' * `$py_call(method, ...)` calls any Python method; `AnnData` and `MuData`
#'   results are wrapped
#'
#' @return A `ReticulateMuData` object
#' @family MuData classes
#' @export
ReticulateMuData <- R6::R6Class(
  "ReticulateMuData", # nolint
  inherit = AbstractMuData,
  cloneable = FALSE,
  private = list(
    .py = NULL,

    .get = function(name) {
      reticulate::py_get_attr(private$.py, name)
    },
    .set = function(name, value) {
      reticulate::py_set_attr(private$.py, name, value)
      invisible(self)
    },
    .mapping = function(name) {
      py_mapping_to_list(private$.get(name))
    },
    .set_mapping = function(name, value) {
      value <- value %||% list()
      if (!is.list(value) || (length(value) && is.null(names(value)))) {
        mudatar_abort("{.field {name}} must be a named list")
      }
      mapping <- private$.get(name)
      for (key in setdiff(py_keys(mapping), names(value))) {
        mapping$`__delitem__`(key)
      }
      for (key in names(value)) {
        mapping$`__setitem__`(key, reticulate::r_to_py(value[[key]]))
      }
      invisible(self)
    },
    .set_df = function(name, value) {
      names <- if (name == "obs") self$obs_names else self$var_names
      df <- private$.prepare_df(value, name, names)
      rownames(df) <- names
      private$.set(name, reticulate::r_to_py(df))
    },
    .backing_label = function() {
      self$filename
    }
  ),
  active = list(
    #' @field obs See [AbstractMuData]
    obs = function(value) {
      if (missing(value)) return(reticulate::py_to_r(private$.get("obs")))
      private$.set_df("obs", value)
    },
    #' @field var See [AbstractMuData]
    var = function(value) {
      if (missing(value)) return(reticulate::py_to_r(private$.get("var")))
      private$.set_df("var", value)
    },
    #' @field obsm See [AbstractMuData]
    obsm = function(value) {
      if (missing(value)) return(private$.mapping("obsm"))
      private$.set_mapping("obsm", check_aligned_mapping(value, "obsm", self$n_obs()))
    },
    #' @field varm See [AbstractMuData]
    varm = function(value) {
      if (missing(value)) return(private$.mapping("varm"))
      private$.set_mapping("varm", check_aligned_mapping(value, "varm", self$n_vars()))
    },
    #' @field obsp See [AbstractMuData]
    obsp = function(value) {
      if (missing(value)) return(private$.mapping("obsp"))
      private$.set_mapping("obsp", check_aligned_mapping(value, "obsp", self$n_obs(), pairwise = TRUE))
    },
    #' @field varp See [AbstractMuData]
    varp = function(value) {
      if (missing(value)) return(private$.mapping("varp"))
      private$.set_mapping("varp", check_aligned_mapping(value, "varp", self$n_vars(), pairwise = TRUE))
    },
    #' @field uns See [AbstractMuData]
    uns = function(value) {
      if (missing(value)) return(reticulate::py_to_r(private$.get("uns")) %||% list())
      private$.set("uns", reticulate::r_to_py(value %||% list()))
    },
    #' @field obsmap See [AbstractMuData] (read-only: rebuilt by `update()`)
    obsmap = function(value) {
      if (!missing(value)) mudatar_readonly_map()
      lapply(private$.mapping("obsmap"), function(m) as.integer(as.vector(m)))
    },
    #' @field varmap See [AbstractMuData] (read-only: rebuilt by `update()`)
    varmap = function(value) {
      if (!missing(value)) mudatar_readonly_map()
      lapply(private$.mapping("varmap"), function(m) as.integer(as.vector(m)))
    },
    #' @field obs_names See [AbstractMuData]
    obs_names = function(value) {
      if (missing(value)) return(py_index_to_r(private$.get("obs_names")))
      private$.set("obs_names", reticulate::r_to_py(as.character(value)))
    },
    #' @field var_names See [AbstractMuData]
    var_names = function(value) {
      if (missing(value)) return(py_index_to_r(private$.get("var_names")))
      private$.set("var_names", reticulate::r_to_py(as.character(value)))
    },
    #' @field axis See [AbstractMuData]
    axis = function(value) {
      if (!missing(value)) mudatar_abort_readonly("axis")
      as.integer(reticulate::py_to_r(private$.get("axis")))
    },
    #' @field filename See [AbstractMuData]
    filename = function(value) {
      if (!missing(value)) {
        mudatar_abort_unsupported("Re-binding {.field filename}", "i" = "Use {.code $py_mudata()$filename} in Python.")
      }
      f <- private$.get("filename")
      if (inherits(f, "python.builtin.NoneType")) NULL else as.character(reticulate::py_str(f))
    },
    #' @field isbacked See [AbstractMuData]
    isbacked = function(value) {
      if (!missing(value)) mudatar_abort_readonly("field")
      isTRUE(reticulate::py_to_r(private$.get("isbacked")))
    },
    #' @field is_view See [AbstractMuData]
    is_view = function(value) {
      if (!missing(value)) mudatar_abort_readonly("field")
      isTRUE(reticulate::py_to_r(private$.get("is_view")))
    }
  ),
  public = list(
    #' @description Wrap a Python `mudata.MuData`
    #' @param py_mudata A Python `mudata.MuData` object, or a
    #'   `ReticulateMuData` (wrapping is idempotent)
    initialize = function(py_mudata) {
      check_reticulate_mudata()
      if (inherits(py_mudata, "ReticulateMuData")) {
        py_mudata <- py_mudata$py_mudata()
      }
      if (!is_py_mudata(py_mudata)) {
        mudatar_abort("{.arg py_mudata} must be a Python {.cls mudata.MuData}, not {.cls {class(py_mudata)}}")
      }
      private$.py <- py_mudata
      invisible(self)
    },

    #' @description The underlying Python `mudata.MuData` (`convert = FALSE`)
    py_mudata = function() {
      private$.py
    },

    #' @description Call a Python method of the MuData; `AnnData` and `MuData`
    #'   results are wrapped as `ReticulateAnnData` / `ReticulateMuData`
    #' @param method Method name
    #' @param ... Arguments passed to the method
    py_call = function(method, ...) {
      wrap_py_return(reticulate::py_get_attr(private$.py, method)(...))
    },

    #' @description See [AbstractMuData]
    mod_names = function() {
      py_keys(private$.get("mod"))
    },

    #' @description A modality as an `anndataR::ReticulateAnnData` (or a
    #'   `ReticulateMuData` for nested MuData)
    #' @param name Modality name
    get_mod = function(name) {
      name <- resolve_mod_name(self$mod_names(), name)
      wrap_py_return(private$.get("mod")$`__getitem__`(name))
    },

    #' @description Add, replace or (with `NULL`) remove a modality, then
    #'   rebuild the indices with Python `MuData.update()`. Non-reticulate
    #'   AnnData objects are converted with `$as_ReticulateAnnData()`.
    #' @param name Modality name
    #' @param value An anndataR AnnData object, a Python `anndata.AnnData`, or
    #'   `NULL`
    set_mod = function(name, value) {
      validate_mod_name(name)
      mod <- private$.get("mod")
      if (is.null(value)) {
        resolve_mod_name(self$mod_names(), name)
        mod$`__delitem__`(name)
      } else {
        mod$`__setitem__`(name, as_py_anndata(value))
      }
      self$update()
      invisible(self)
    },

    #' @description See [AbstractMuData]; runs Python `MuData.update_obs()`
    update_obs = function() {
      private$.py$update_obs()
      invisible(self)
    },
    #' @description See [AbstractMuData]; runs Python `MuData.update_var()`
    update_var = function() {
      private$.py$update_var()
      invisible(self)
    },
    #' @description See [AbstractMuData]; runs Python `MuData.update()`
    update = function() {
      private$.py$update()
      invisible(self)
    },

    #' @description Copy modality `obs` columns to the global `obs` (Python
    #'   `MuData.pull_obs()`; arguments as in Python)
    #' @param columns,mods Character vectors or `NULL`
    #' @param ... Further arguments of Python `pull_obs()` (`common`,
    #'   `join_common`, `nonunique`, `join_nonunique`, `unique`,
    #'   `prefix_unique`, `drop`, `only_drop`)
    pull_obs = function(columns = NULL, mods = NULL, ...) {
      private$.py$pull_obs(columns = py_list_or_null(columns), mods = py_list_or_null(mods), ...)
      invisible(self)
    },
    #' @description As `pull_obs()` for `var`
    #' @param columns,mods Character vectors or `NULL`
    #' @param ... Further arguments of Python `pull_var()`
    pull_var = function(columns = NULL, mods = NULL, ...) {
      private$.py$pull_var(columns = py_list_or_null(columns), mods = py_list_or_null(mods), ...)
      invisible(self)
    },
    #' @description Copy global `obs` columns to the modalities (Python
    #'   `MuData.push_obs()`)
    #' @param columns,mods Character vectors or `NULL`
    #' @param ... Further arguments of Python `push_obs()` (`common`,
    #'   `prefixed`, `drop`, `only_drop`)
    push_obs = function(columns = NULL, mods = NULL, ...) {
      private$.py$push_obs(columns = py_list_or_null(columns), mods = py_list_or_null(mods), ...)
      invisible(self)
    },
    #' @description As `push_obs()` for `var`
    #' @param columns,mods Character vectors or `NULL`
    #' @param ... Further arguments of Python `push_var()`
    push_var = function(columns = NULL, mods = NULL, ...) {
      private$.py$push_var(columns = py_list_or_null(columns), mods = py_list_or_null(mods), ...)
      invisible(self)
    },

    #' @description A deep copy, as a new `ReticulateMuData` (Python
    #'   `MuData.copy()`)
    copy = function() {
      ReticulateMuData$new(private$.py$copy())
    },

    #' @description Return `self` (already Python-backed)
    as_ReticulateMuData = function() {
      # nolint start: object_name_linter
      self
      # nolint end: object_name_linter
    }
  )
)

mudatar_readonly_map <- function(call = rlang::caller_env()) {
  mudatar_abort(
    "The maps of a {.cls ReticulateMuData} are maintained by Python; call {.code $update()} instead.",
    class = "mudatar_readonly_error",
    call = call
  )
}

check_reticulate_mudata <- function(call = rlang::caller_env()) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    mudatar_abort("Package {.pkg reticulate} is required for {.cls ReticulateMuData}", call = call)
  }
  if (!reticulate::py_module_available("mudata")) {
    mudatar_abort(
      c(
        "Python module {.pkg mudata} is not available to reticulate.",
        "i" = "{.cls ReticulateMuData} needs Python; native {.fn read_h5mu}/{.fn write_h5mu} do not."
      ),
      call = call
    )
  }
}

is_py_mudata <- function(x) {
  inherits(x, "python.builtin.object") &&
    reticulate::py_to_r(reticulate::import_builtins(convert = FALSE)$isinstance(
      x,
      reticulate::import("mudata", convert = FALSE)$MuData
    ))
}

is_py_anndata <- function(x) {
  inherits(x, "python.builtin.object") &&
    reticulate::py_to_r(reticulate::import_builtins(convert = FALSE)$isinstance(
      x,
      reticulate::import("anndata", convert = FALSE)$AnnData
    ))
}

# Python AnnData -> anndataR::ReticulateAnnData (public py_to_r() method of
# anndataR); Python MuData -> ReticulateMuData; anything else unchanged
wrap_py_return <- function(obj) {
  if (!inherits(obj, "python.builtin.object")) {
    return(obj)
  }
  if (is_py_mudata(obj)) {
    return(ReticulateMuData$new(obj))
  }
  if (is_py_anndata(obj)) {
    return(reticulate::py_to_r(obj))
  }
  obj
}

as_py_anndata <- function(value, call = rlang::caller_env()) {
  if (inherits(value, "ReticulateAnnData")) {
    return(value$py_anndata())
  }
  if (inherits(value, "AbstractAnnData")) {
    return(anndata_to_py(value))
  }
  if (is_py_anndata(value)) {
    return(value)
  }
  mudatar_abort("A modality must be an anndataR AnnData or a Python AnnData, not {.cls {class(value)}}", call = call)
}

# An R AnnData as a Python anndata.AnnData: written natively to a temporary
# H5AD (Awkward Arrays included) and read by Python anndata. Lossless, and
# independent of anndataR's object-by-object conversion (which does not
# handle Awkward Arrays and fails for zero-length axes with anndata 0.13).
anndata_to_py <- function(adata) {
  ad <- reticulate::import("anndata", convert = FALSE)
  parts <- split_awkward(adata)
  tmp <- tempfile("mudatar-py-", fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  anndatar_write("h5mu", parts$adata, tmp)
  for (slot in names(parts$awkward)) {
    for (k in names(parts$awkward[[slot]])) {
      store_write_awkward("h5mu", tmp, paste0(slot, "/", k), parts$awkward[[slot]][[k]])
    }
  }
  ad$read_h5ad(tmp)
}

py_keys <- function(mapping) {
  bi <- reticulate::import_builtins(convert = FALSE)
  as.character(reticulate::py_to_r(bi$list(mapping$keys())))
}

py_index_to_r <- function(index) {
  bi <- reticulate::import_builtins(convert = FALSE)
  as.character(reticulate::py_to_r(bi$list(index$astype("str"))))
}

py_mapping_to_list <- function(mapping) {
  keys <- py_keys(mapping)
  out <- lapply(keys, function(k) {
    value <- mapping$`__getitem__`(k)
    # keep values reticulate cannot convert (e.g. awkward arrays) as Python
    # objects rather than dropping them
    tryCatch(reticulate::py_to_r(value), error = function(e) value)
  })
  names(out) <- keys
  out
}

py_list_or_null <- function(x) {
  if (is.null(x)) NULL else reticulate::r_to_py(as.list(as.character(x)))
}

#' Convert to a ReticulateMuData
#'
#' Wrap a Python `mudata.MuData`, or build one from any MuData object: each
#' modality is converted with `anndataR`'s `$as_ReticulateAnnData()`, and the
#' global `obs`, `var`, `uns` and non-modality `obsm`/`varm` entries are
#' copied. Python recomputes the maps, which match mudataR's (see the
#' differential tests).
#'
#' @param x A MuData object, a Python `mudata.MuData` or a `ReticulateMuData`
#'
#' @return A [ReticulateMuData]
#' @family MuData classes
#' @export
as_ReticulateMuData <- function(x) {
  # nolint start: object_name_linter
  if (inherits(x, "AbstractMuData")) {
    return(x$as_ReticulateMuData())
  }
  ReticulateMuData$new(x)
  # nolint end: object_name_linter
}

#' @rdname is_mudata
#' @export
is_ReticulateMuData <- function(x) {
  # nolint start: object_name_linter
  inherits(x, "ReticulateMuData")
  # nolint end: object_name_linter
}

mudata_to_reticulate <- function(x, call = rlang::caller_env()) {
  check_reticulate_mudata(call = call)
  md <- reticulate::import("mudata", convert = FALSE)
  mod_names <- x$mod_names()
  py_mods <- reticulate::py_dict(
    keys = as.list(mod_names),
    values = lapply(mod_names, function(m) as_py_anndata(x[[m]], call = call)),
    convert = FALSE
  )
  obs <- x$obs
  rownames(obs) <- x$obs_names
  var <- x$var
  rownames(var) <- x$var_names
  # As Python's own read_h5mu(): pass the global slots and let
  # MuData.update() reconcile them with the modalities (it keeps the given
  # order and drops names no modality has)
  args <- list(
    data = py_mods,
    obs = reticulate::r_to_py(obs),
    var = reticulate::r_to_py(var),
    uns = reticulate::r_to_py(x$uns %||% list()),
    axis = as.integer(x$axis)
  )
  for (slot in c("obsm", "varm", "obsp", "varp")) {
    value <- if (slot %in% c("obsm", "varm")) {
      global_aligned(.subset2(x, slot), mod_names)
    } else {
      .subset2(x, slot)
    }
    if (length(value)) {
      args[[slot]] <- reticulate::py_dict(
        keys = as.list(names(value)),
        values = lapply(value, reticulate::r_to_py),
        convert = FALSE
      )
    }
  }
  ReticulateMuData$new(do.call(md$MuData, args))
}
