#' Modality mapping
#'
#' @description
#' The `$mod` field of a MuData object. It behaves like an ordered, named list
#' of AnnData objects (`names()`, `length()`, `[[`, `as.list()`, `lapply()`)
#' but assignment goes through the owning MuData object, so that modality
#' names are validated and the global maps are rebuilt.
#'
#' @param x A `ModalityMapping`
#' @param i,name Modality name (or position for `[[`)
#' @param value An `anndataR` AnnData object, or `NULL` to remove a modality
#' @param ... Unused
#'
#' @return `[[` returns the modality's AnnData object. Each backend returns
#'   its matching `anndataR` class: `InMemoryAnnData` for [InMemoryMuData],
#'   a rooted `HDF5AnnData` for [HDF5MuData] and a rooted `ZarrAnnData` for
#'   [ZarrMuData].
#' @name ModalityMapping
NULL

new_modality_mapping <- function(owner) {
  structure(list(owner = owner), class = "ModalityMapping")
}

mapping_owner <- function(x) {
  .subset2(x, "owner")
}

#' @rdname ModalityMapping
#' @export
names.ModalityMapping <- function(x) {
  mapping_owner(x)$mod_names()
}

#' @rdname ModalityMapping
#' @export
length.ModalityMapping <- function(x) {
  length(names(x))
}

#' @rdname ModalityMapping
#' @export
`[[.ModalityMapping` <- function(x, i, ...) {
  mapping_owner(x)$get_mod(resolve_mod_name(names(x), i))
}

#' @rdname ModalityMapping
#' @export
`$.ModalityMapping` <- function(x, name) {
  mapping_owner(x)$get_mod(resolve_mod_name(names(x), name, exact = TRUE))
}

#' @rdname ModalityMapping
#' @export
`[[<-.ModalityMapping` <- function(x, i, value) {
  if (is.numeric(i)) {
    i <- resolve_mod_name(names(x), i)
  }
  mapping_owner(x)$set_mod(i, value)
  x
}

#' @rdname ModalityMapping
#' @export
`$<-.ModalityMapping` <- function(x, name, value) {
  mapping_owner(x)$set_mod(name, value)
  x
}

#' @rdname ModalityMapping
#' @export
as.list.ModalityMapping <- function(x, ...) {
  nms <- names(x)
  out <- lapply(nms, function(n) mapping_owner(x)$get_mod(n))
  names(out) <- nms
  out
}

#' @export
print.ModalityMapping <- function(x, ...) {
  nms <- names(x)
  cat("<ModalityMapping> ", length(nms), " modalit", if (length(nms) == 1) "y" else "ies", "\n", sep = "")
  for (n in nms) {
    cat("  ", n, "\n", sep = "")
  }
  invisible(x)
}

resolve_mod_name <- function(nms, i, exact = TRUE, call = rlang::caller_env()) {
  if (is.numeric(i)) {
    if (length(i) != 1L || is.na(i) || i < 1 || i > length(nms)) {
      mudatar_abort("Modality index {.val {i}} is out of bounds", call = call)
    }
    return(nms[[i]])
  }
  if (!rlang::is_string(i)) {
    mudatar_abort("Modalities are selected by a single name or position", call = call)
  }
  if (!i %in% nms) {
    mudatar_abort(
      c("No modality named {.val {i}}", "i" = "Available: {.val {nms}}"),
      class = "mudatar_missing_modality",
      call = call
    )
  }
  i
}

validate_mod_name <- function(name, call = rlang::caller_env()) {
  if (!rlang::is_string(name) || is.na(name) || !nzchar(name)) {
    mudatar_abort("Modality names must be non-empty strings", call = call)
  }
  if (grepl("/", name, fixed = TRUE)) {
    mudatar_abort("Modality name {.val {name}} must not contain {.val /}", call = call)
  }
  name
}
