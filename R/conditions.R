#' Abort with a classed mudataR condition
#' @noRd
mudatar_abort <- function(message, ..., class = NULL, call = rlang::caller_env(), .envir = parent.frame()) {
  cli_abort(
    message,
    ...,
    class = c(class, "mudatar_error"),
    call = call,
    .envir = .envir
  )
}

#' Abort because a MuData feature is not implemented yet
#'
#' Used instead of silently approximating Python behaviour; each use is listed
#' in `docs/parity.md`.
#' @noRd
mudatar_abort_unsupported <- function(what, ..., call = rlang::caller_env(), .envir = parent.frame()) {
  mudatar_abort(
    c(paste(what, "is not supported by mudataR yet."), ...),
    class = "mudatar_unsupported",
    call = call,
    .envir = .envir
  )
}

#' Abort because the object (or its store) was closed
#' @noRd
mudatar_abort_closed <- function(path, call = rlang::caller_env()) {
  mudatar_abort(
    c(
      "The MuData store {.file {path}} has been closed.",
      "i" = "Call {.code $open()} to access it again."
    ),
    class = "mudatar_closed_error",
    call = call
  )
}

#' Abort because the object is read-only
#' @noRd
mudatar_abort_readonly <- function(what = "object", call = rlang::caller_env()) {
  mudatar_abort(
    "This {what} is read-only and cannot be modified.",
    class = "mudatar_readonly_error",
    call = call
  )
}
