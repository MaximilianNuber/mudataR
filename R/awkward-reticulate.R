#' Convert Awkward Arrays between R and Python
#'
#' `r_to_py()` turns an [AwkwardArrowArray] (or an [AwkwardFrame]) into a
#' Python `awkward.Array` with `ak.from_buffers()`, and `py_to_r()` turns an
#' `awkward.Array` into an [AwkwardArrowArray] with `ak.to_buffers()`. Both
#' directions keep the form and buffers, so ragged data such as scirpy AIRR
#' data are lossless (64-bit integers are exact up to 2^53).
#'
#' @param x An `AwkwardArrowArray`/`AwkwardFrame` or a Python `awkward.Array`
#' @param convert Passed to reticulate
#' @return The converted object
#' @name awkward-reticulate
#' @keywords internal
NULL

#' @rdname awkward-reticulate
#' @exportS3Method reticulate::r_to_py
r_to_py.AwkwardArrowArray <- function(x, convert = FALSE) {
  ak <- reticulate::import("awkward", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  buffers <- awkward_prepare_buffers(x, awk_field(x, "source")$path, "h5mu")
  values <- lapply(buffers, function(b) {
    # a 1-d R array always becomes a 1-d numpy array (never a scalar)
    np$asarray(
      reticulate::r_to_py(array(b$value, dim = length(b$value))),
      dtype = b$dtype
    )
  })
  container <- reticulate::py_dict(
    keys = as.list(names(buffers)),
    values = unname(values),
    convert = FALSE
  )
  ak$from_buffers(
    awk_field(x, "form_json"),
    as.integer(length(x)),
    container
  )
}

#' @rdname awkward-reticulate
#' @exportS3Method reticulate::r_to_py
r_to_py.AwkwardFrame <- function(x, convert = FALSE) {
  r_to_py.AwkwardArrowArray(as_AwkwardArrowArray(x), convert = convert)
}

#' @rdname awkward-reticulate
#' @exportS3Method reticulate::py_to_r
py_to_r.awkward.highlevel.Array <- function(x) {
  ak <- reticulate::import("awkward", convert = FALSE)
  bi <- reticulate::import_builtins(convert = FALSE)
  parts <- ak$to_buffers(x)
  form_json <- reticulate::py_to_r(parts[0]$to_json())
  length <- reticulate::py_to_r(parts[1])
  container <- parts[2]
  keys <- as.character(reticulate::py_to_r(bi$list(container$keys())))
  buffers <- lapply(keys, function(k) {
    reticulate::py_to_r(container$`__getitem__`(k)) |> as.vector()
  })
  names(buffers) <- keys
  form <- awkward_parse_form(form_json)
  new_awkward_arrow_array(
    form = form,
    form_json = form_json,
    length = length,
    buffers = buffers,
    source = list(store_kind = "python", path = NULL, encoding_version = "0.1.0")
  )
}
