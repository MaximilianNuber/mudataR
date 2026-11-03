# Observation/variable maps
#
# `obsmap[[mod]]` has one entry per global observation: 0 means the
# observation is absent from `mod`, a positive value is its one-based position
# in `mod` (the same convention as Python mudata and Muon.jl). Maps are stored
# as unsigned 32-bit integers on disk and as R integers in memory.

#' Validate a set of maps against global names and modalities
#' @noRd
validate_maps <- function(maps, attr, global_names, mod_sizes, call = rlang::caller_env()) {
  label <- paste0(attr, "map")
  if (!is.list(maps) || (length(maps) && is.null(names(maps)))) {
    mudatar_abort("{.field {label}} must be a named list", call = call)
  }
  missing <- setdiff(names(mod_sizes), names(maps))
  if (length(missing)) {
    mudatar_abort("{.field {label}} is missing modalities {.val {missing}}", call = call)
  }
  extra <- setdiff(names(maps), names(mod_sizes))
  if (length(extra)) {
    mudatar_abort("{.field {label}} has entries for unknown modalities {.val {extra}}", call = call)
  }
  n <- length(global_names)
  out <- lapply(names(mod_sizes), function(m) {
    map <- maps[[m]]
    if (is.matrix(map) && ncol(map) == 1L) {
      map <- map[, 1]
    }
    map <- as.vector(map)
    if (!is.numeric(map) || anyNA(map)) {
      mudatar_abort("{.field {label}[[{.val {m}}]]} must be a numeric vector without missing values", call = call)
    }
    if (length(map) != n) {
      mudatar_abort(
        "{.field {label}[[{.val {m}}]]} has length {length(map)} but there are {n} global {attr} names",
        call = call
      )
    }
    if (any(map < 0) || any(map != round(map))) {
      mudatar_abort("{.field {label}[[{.val {m}}]]} must contain non-negative integers", call = call)
    }
    if (any(map > mod_sizes[[m]])) {
      mudatar_abort(
        "{.field {label}[[{.val {m}}]]} points beyond the {mod_sizes[[m]]} {attr} of modality {.val {m}}",
        call = call
      )
    }
    present <- map[map > 0]
    if (anyDuplicated(present)) {
      mudatar_abort("{.field {label}[[{.val {m}}]]} maps several global {attr} to the same modality {attr}", call = call)
    }
    as.integer(map)
  })
  names(out) <- names(mod_sizes)
  out
}

#' Check maps are consistent with modality names
#' @noRd
check_maps_match_names <- function(maps, attr, global_names, mod_names, call = rlang::caller_env()) {
  for (m in names(maps)) {
    map <- maps[[m]]
    present <- map > 0
    # some backends return NULL rather than character(0) for empty names
    if (!identical(as.character(global_names[present]), as.character(mod_names[[m]][map[present]]))) {
      mudatar_abort(
        c(
          "{.field {attr}map[[{.val {m}}]]} does not match the {attr} names of modality {.val {m}}",
          "i" = "Call {.code $update()} to rebuild the maps."
        ),
        call = call
      )
    }
  }
  invisible(TRUE)
}

#' Occurrence-qualified keys, so duplicated names can be matched positionally
#' @noRd
occurrence_keys <- function(x) {
  if (!anyDuplicated(x)) {
    return(x)
  }
  paste0(x, "\r", stats::ave(seq_along(x), x, FUN = seq_along))
}

is_shared_attr <- function(attr, axis) {
  (attr == "obs" && axis %in% c(0L, -1L)) || (attr == "var" && axis %in% c(1L, -1L))
}

#' Rebuild global names and maps for one attribute
#'
#' Port of the main (no within-modality duplicates) case of
#' `mudata.MuData._update_attr()` (mudata 0.4.1):
#'
#' * on the shared axis the global names are the union of modality names in
#'   order of first appearance;
#' * on the other axis they are the concatenation of modality names
#'   (duplicates across modalities are kept, as in Python, with a warning);
#' * names that were already global keep their previous order, new names are
#'   appended;
#' * `index_order` gives, for each new global name, its previous position (or
#'   `NA`), used to reindex global columns and `obsm`/`obsp`.
#'
#' @param mod_names Named list of modality names along `attr`
#' @param prev_names Previous global names (possibly empty)
#' @return A list with `names`, `maps` and `index_order` (NULL if there were
#'   no previous names)
#' @noRd
compute_global_attr <- function(mod_names, attr, axis, prev_names = character(), call = rlang::caller_env()) {
  shared <- is_shared_attr(attr, axis)
  within_dup <- vapply(mod_names, function(n) anyDuplicated(n) > 0, logical(1))
  if (any(within_dup)) {
    mudatar_abort_unsupported(
      "Updating MuData with duplicated {attr}_names within a modality",
      "i" = "Modalities with duplicates: {.val {names(mod_names)[within_dup]}}.",
      "i" = "Make them unique first, e.g. with {.code make.unique()}.",
      call = call
    )
  }

  if (shared) {
    new_names <- unique(unlist(mod_names, use.names = FALSE))
    if (is.null(new_names)) new_names <- character()
    maps <- lapply(mod_names, function(n) {
      m <- match(new_names, n)
      m[is.na(m)] <- 0L
      as.integer(m)
    })
  } else {
    sizes <- lengths(mod_names)
    new_names <- unlist(mod_names, use.names = FALSE)
    if (is.null(new_names)) new_names <- character()
    ends <- cumsum(sizes)
    starts <- ends - sizes
    total <- sum(sizes)
    maps <- lapply(seq_along(mod_names), function(k) {
      m <- integer(total)
      if (sizes[[k]] > 0) m[starts[[k]] + seq_len(sizes[[k]])] <- seq_len(sizes[[k]])
      m
    })
    names(maps) <- names(mod_names)
    if (anyDuplicated(new_names)) {
      cli_warn(
        "{attr}_names are not unique. To make them unique, call {.code ${attr}_names_make_unique()}."
      )
    }
  }

  index_order <- NULL
  if (length(prev_names) > 0) {
    prev_keys <- occurrence_keys(prev_names)
    new_keys <- occurrence_keys(new_names)
    kept <- prev_keys[prev_keys %in% new_keys]
    added <- new_keys[!new_keys %in% prev_keys]
    ordered <- c(kept, added)
    perm <- match(ordered, new_keys)
    new_names <- new_names[perm]
    maps <- lapply(maps, function(m) m[perm])
    index_order <- match(ordered, prev_keys)

    can_update <- length(added) == 0 ||
      length(kept) == length(prev_keys) ||
      length(ordered) == length(prev_keys) ||
      (!shared && length(ordered) > length(prev_keys))
    if (!can_update) {
      mudatar_abort(
        c(
          "{attr}_names seem to have been renamed and filtered at the same time.",
          "i" = "There is no way to restore the order; re-create the MuData from its modalities: {.code MuData(as.list(mdata$mod))}."
        ),
        call = call
      )
    }
  }

  list(names = new_names, maps = maps, index_order = index_order)
}

#' Reindex an aligned element by `index_order` (NA rows become missing)
#' @noRd
reindex_rows <- function(x, index_order, call = rlang::caller_env()) {
  if (is.null(x) || is.null(index_order)) {
    return(x)
  }
  if (is_awkward_value(x)) {
    if (anyNA(index_order)) {
      mudatar_abort_unsupported("Adding observations to an Awkward Array element", call = call)
    }
    return(x[index_order])
  }
  if (is.data.frame(x)) {
    out <- x[index_order, , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  if (inherits(x, "sparseMatrix") && anyNA(index_order)) {
    mudatar_abort_unsupported("Filling new rows of a sparse aligned matrix", call = call)
  }
  if (is.null(dim(x))) {
    return(x[index_order])
  }
  x[index_order, , drop = FALSE]
}

reindex_pairwise <- function(x, index_order, call = rlang::caller_env()) {
  if (is.null(x) || is.null(index_order)) {
    return(x)
  }
  if (inherits(x, "sparseMatrix") && anyNA(index_order)) {
    mudatar_abort_unsupported("Filling new rows of a sparse pairwise matrix", call = call)
  }
  x[index_order, index_order, drop = FALSE]
}
