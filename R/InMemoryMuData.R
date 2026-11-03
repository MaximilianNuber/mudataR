#' @title InMemoryMuData
#'
#' @description
#' A MuData object held in memory. Every modality is an `anndataR`
#' `InMemoryAnnData`: other AnnData objects (backed, Python-backed or views)
#' are materialised when added. Use [MuData()] to construct one.
#'
#' @return An `InMemoryMuData` object
#' @family MuData classes
#' @export
InMemoryMuData <- R6::R6Class(
  "InMemoryMuData", # nolint
  inherit = AbstractMuData,
  cloneable = FALSE,
  private = list(
    .mod = list(),
    .obs = NULL,
    .var = NULL,
    .obsm = list(),
    .varm = list(),
    .obsp = list(),
    .varp = list(),
    .uns = list(),
    .obsmap = list(),
    .varmap = list(),
    .obs_names = character(),
    .var_names = character(),
    .axis = 0L,

    .replace_attr = function(attr, names, df, maps, m, p) {
      private[[paste0(".", attr, "_names")]] <- names
      private[[paste0(".", attr)]] <- df
      private[[paste0(".", attr, "map")]] <- maps
      private[[paste0(".", attr, "m")]] <- m
      private[[paste0(".", attr, "p")]] <- p
      invisible(NULL)
    },

    .validate_aligned = function(value, attr, pairwise = FALSE) {
      if (is.null(value)) {
        return(list())
      }
      if (!is.list(value) || is.data.frame(value) || (length(value) && is.null(names(value)))) {
        mudatar_abort("{.field {attr}} must be a named list")
      }
      n <- if (startsWith(attr, "obs")) length(private$.obs_names) else length(private$.var_names)
      for (key in names(value)) {
        x <- value[[key]]
        len <- if (inherits(x, "AwkwardArrowArray")) length(x) else NROW(x)
        if (len != n || (pairwise && NCOL(x) != n)) {
          mudatar_abort(
            "{.field {attr}[[{.val {key}}]]} must be aligned to the {n} global {substr(attr, 1, 3)} names"
          )
        }
        if (!inherits(x, "AwkwardArrowArray") && !is.null(dim(x)) && !is.data.frame(x)) {
          dimnames(x) <- NULL
          value[[key]] <- x
        }
      }
      value
    }
  ),
  active = list(
    #' @field obs See [AbstractMuData]
    obs = function(value) {
      if (missing(value)) {
        return(private$.add_df_names(private$.obs, private$.obs_names))
      }
      private$.obs <- private$.prepare_df(value, "obs", private$.obs_names)
    },
    #' @field var See [AbstractMuData]
    var = function(value) {
      if (missing(value)) {
        return(private$.add_df_names(private$.var, private$.var_names))
      }
      private$.var <- private$.prepare_df(value, "var", private$.var_names)
    },
    #' @field obsm See [AbstractMuData]
    obsm = function(value) {
      if (missing(value)) return(private$.obsm)
      private$.obsm <- private$.validate_aligned(value, "obsm")
    },
    #' @field varm See [AbstractMuData]
    varm = function(value) {
      if (missing(value)) return(private$.varm)
      private$.varm <- private$.validate_aligned(value, "varm")
    },
    #' @field obsp See [AbstractMuData]
    obsp = function(value) {
      if (missing(value)) return(private$.obsp)
      private$.obsp <- private$.validate_aligned(value, "obsp", pairwise = TRUE)
    },
    #' @field varp See [AbstractMuData]
    varp = function(value) {
      if (missing(value)) return(private$.varp)
      private$.varp <- private$.validate_aligned(value, "varp", pairwise = TRUE)
    },
    #' @field uns See [AbstractMuData]
    uns = function(value) {
      if (missing(value)) return(private$.uns)
      if (!is.null(value) && (!is.list(value) || (length(value) && is.null(names(value))))) {
        mudatar_abort("{.field uns} must be a named list")
      }
      private$.uns <- value %||% list()
    },
    #' @field obsmap See [AbstractMuData]
    obsmap = function(value) {
      if (missing(value)) return(private$.obsmap)
      private$.obsmap <- validate_maps(
        value,
        "obs",
        private$.obs_names,
        vapply(private$.mod, function(a) a$n_obs(), integer(1))
      )
    },
    #' @field varmap See [AbstractMuData]
    varmap = function(value) {
      if (missing(value)) return(private$.varmap)
      private$.varmap <- validate_maps(
        value,
        "var",
        private$.var_names,
        vapply(private$.mod, function(a) a$n_vars(), integer(1))
      )
    },
    #' @field obs_names See [AbstractMuData]. Assigning names of the same
    #'   length renames global observations (modalities are not renamed).
    obs_names = function(value) {
      if (missing(value)) return(private$.obs_names)
      value <- as.character(value)
      if (length(value) != length(private$.obs_names) || anyNA(value)) {
        mudatar_abort("{.field obs_names} must be {length(private$.obs_names)} non-missing names")
      }
      private$.obs_names <- value
    },
    #' @field var_names See [AbstractMuData]
    var_names = function(value) {
      if (missing(value)) return(private$.var_names)
      value <- as.character(value)
      if (length(value) != length(private$.var_names) || anyNA(value)) {
        mudatar_abort("{.field var_names} must be {length(private$.var_names)} non-missing names")
      }
      private$.var_names <- value
    },
    #' @field axis See [AbstractMuData]
    axis = function(value) {
      if (!missing(value)) {
        mudatar_abort_readonly("axis")
      }
      private$.axis
    }
  ),
  public = list(
    #' @description Create an `InMemoryMuData`. Prefer [MuData()].
    #' @param mod Named list of anndataR AnnData objects
    #' @param obs,var,obsm,varm,obsp,varp,uns Global slots
    #' @param obsmap,varmap Optional maps; when given they are validated
    #'   against `obs_names`/`var_names` and the modalities instead of being
    #'   recomputed
    #' @param obs_names,var_names Optional global names. Without maps they fix
    #'   the order of the global index (as Python does with a supplied `obs`)
    #' @param axis `0`, `1` or `-1`
    initialize = function(
      mod,
      obs = NULL,
      var = NULL,
      obsm = NULL,
      varm = NULL,
      obsp = NULL,
      varp = NULL,
      uns = NULL,
      obsmap = NULL,
      varmap = NULL,
      obs_names = NULL,
      var_names = NULL,
      axis = 0L
    ) {
      private$.axis <- check_axis(axis)
      private$.mod <- as_in_memory_modalities(check_modalities(mod))

      for (attr in c("obs", "var")) {
        df <- if (attr == "obs") obs else var
        names <- if (attr == "obs") obs_names else var_names
        maps <- if (attr == "obs") obsmap else varmap
        if (is.null(names) && is.data.frame(df) && .row_names_info(df) > 0) {
          names <- rownames(df)
        }
        mod_names <- lapply(private$.mod, function(a) a[[paste0(attr, "_names")]])
        if (!is.null(maps)) {
          if (is.null(names)) {
            mudatar_abort("{.arg {attr}map} requires {.arg {attr}_names}")
          }
          sizes <- lengths(mod_names)
          maps <- validate_maps(maps, attr, names, sizes)
          check_maps_match_names(maps, attr, names, mod_names)
          private[[paste0(".", attr, "_names")]] <- as.character(names)
          private[[paste0(".", attr, "map")]] <- maps
          private[[paste0(".", attr)]] <- private$.prepare_df(df, attr, names)
        } else {
          prev <- as.character(names %||% character())
          res <- compute_global_attr(mod_names, attr, private$.axis, prev)
          if (length(prev) && !is.null(df)) {
            df <- reindex_rows(private$.prepare_df(df, attr, prev), res$index_order)
          } else if (!is.null(df) && length(prev) == 0 && nrow(df) != length(res$names)) {
            mudatar_abort("{.arg {attr}} must have {length(res$names)} rows")
          }
          private[[paste0(".", attr, "_names")]] <- res$names
          private[[paste0(".", attr, "map")]] <- res$maps
          private[[paste0(".", attr)]] <- private$.prepare_df(df, attr, res$names)
        }
      }

      self$obsm <- obsm
      self$varm <- varm
      self$obsp <- obsp
      self$varp <- varp
      self$uns <- uns
      # Python records modality membership in obsm/varm
      for (m in names(private$.mod)) {
        if (is.null(private$.obsm[[m]])) private$.obsm[[m]] <- private$.obsmap[[m]] > 0L
        if (is.null(private$.varm[[m]])) private$.varm[[m]] <- private$.varmap[[m]] > 0L
      }
      invisible(self)
    },

    #' @description See [AbstractMuData]
    mod_names = function() {
      names(private$.mod)
    },

    #' @description See [AbstractMuData]
    #' @param name Modality name
    get_mod = function(name) {
      private$.mod[[resolve_mod_name(names(private$.mod), name)]]
    },

    #' @description See [AbstractMuData]
    #' @param name Modality name
    #' @param value An anndataR AnnData object or `NULL`
    set_mod = function(name, value) {
      validate_mod_name(name)
      mods <- private$.mod
      if (is.null(value)) {
        resolve_mod_name(names(mods), name)
        mods[[name]] <- NULL
      } else {
        mods[[name]] <- value
      }
      mods <- as_in_memory_modalities(check_modalities(mods, allow_empty = TRUE))
      old <- private$.mod
      private$.mod <- mods
      tryCatch(
        {
          for (attr in c("obs", "var")) {
            m <- private[[paste0(".", attr, "m")]]
            if (is.null(value)) m[[name]] <- NULL
            private[[paste0(".", attr, "m")]] <- m
            maps <- private[[paste0(".", attr, "map")]]
            maps[[name]] <- NULL
            private[[paste0(".", attr, "map")]] <- maps
          }
          self$update()
        },
        error = function(e) {
          private$.mod <- old
          stop(e)
        }
      )
      invisible(self)
    }
  )
)

#' Every modality of an InMemoryMuData is an anndataR InMemoryAnnData
#'
#' Backed (HDF5AnnData/ZarrAnnData), Python-backed (ReticulateAnnData) and
#' view (AnnDataView) modalities are materialised, so an in-memory MuData never
#' depends on a file handle or a Python session.
#' @noRd
as_in_memory_modalities <- function(mod) {
  # keeps Awkward Arrays (as AwkwardFrame), also of backed modalities
  lapply(mod, as_in_memory_anndata)
}

check_axis <- function(axis, call = rlang::caller_env()) {
  axis <- suppressWarnings(as.integer(axis))
  allowed <- .mudata_axes
  if (length(axis) != 1L || is.na(axis) || !axis %in% allowed) {
    mudatar_abort("{.arg axis} must be one of {.val {allowed}}", call = call)
  }
  axis
}

check_modalities <- function(mod, allow_empty = FALSE, call = rlang::caller_env()) {
  if (inherits(mod, "ModalityMapping")) {
    mod <- as.list(mod)
  }
  if (!is.list(mod) || (!allow_empty && length(mod) == 0)) {
    mudatar_abort("{.arg mod} must be a non-empty named list of AnnData objects", call = call)
  }
  nms <- names(mod)
  if (length(mod) && (is.null(nms) || anyNA(nms) || any(!nzchar(nms)))) {
    mudatar_abort("Every modality must be named", call = call)
  }
  if (anyDuplicated(nms)) {
    mudatar_abort("Modality names must be unique; duplicated: {.val {unique(nms[duplicated(nms)])}}", call = call)
  }
  for (n in nms) {
    validate_mod_name(n, call = call)
    if (inherits(mod[[n]], "AbstractMuData")) {
      mudatar_abort_unsupported("Nested MuData modalities", call = call)
    }
    if (!inherits(mod[[n]], "AbstractAnnData")) {
      mudatar_abort(
        "Modality {.val {n}} must be an {.pkg anndataR} AnnData object, not {.cls {class(mod[[n]])}}",
        call = call
      )
    }
  }
  mod
}

#' Create a MuData object
#'
#' Build an [InMemoryMuData] from named `anndataR` AnnData objects. Global
#' observation and variable names and the `obsmap`/`varmap` maps are computed
#' from the modalities as in Python `mudata.MuData()`.
#'
#' @param data A named list of `anndataR` AnnData objects
#' @param ... Global slots passed to [InMemoryMuData]: `obs`, `var`, `obsm`,
#'   `varm`, `obsp`, `varp`, `uns`, `obsmap`, `varmap`, `obs_names`,
#'   `var_names`
#' @param axis `0` for shared observations (default), `1` for shared
#'   variables, `-1` for both
#'
#' @return An [InMemoryMuData]
#' @export
#' @examples
#' rna <- anndataR::AnnData(
#'   obs = data.frame(row.names = c("c1", "c2", "c3")),
#'   var = data.frame(row.names = c("g1", "g2"))
#' )
#' prot <- anndataR::AnnData(
#'   obs = data.frame(row.names = c("c2", "c3", "c4")),
#'   var = data.frame(row.names = "p1")
#' )
#' mdata <- MuData(list(rna = rna, prot = prot))
#' mdata
#' mdata$obsmap
MuData <- function(data, ..., axis = 0L) {
  # nolint start: object_name_linter
  InMemoryMuData$new(mod = data, ..., axis = axis)
  # nolint end: object_name_linter
}
