# Conversion of MuData objects to Bioconductor containers
#
# Each modality is converted with anndataR's own as_SingleCellExperiment(), so
# every anndataR mapping argument (x_mapping, assays_mapping, colData_mapping,
# rowData_mapping, reducedDims_mapping, colPairs_mapping, rowPairs_mapping,
# metadata_mapping) can be given per modality through `mod_args`. mudataR only
# adds the MuData layer: aligning modalities through obsmap, the global
# annotation and Awkward Array data that Bioconductor containers cannot hold
# (kept as AwkwardArrowArray objects in `metadata(<experiment>)$awkward_obsm`).

.anndatar_sce_args <- c(
  "x_mapping",
  "assays_mapping",
  "colData_mapping",
  "rowData_mapping",
  "reducedDims_mapping",
  "colPairs_mapping",
  "rowPairs_mapping",
  "metadata_mapping"
)

check_conversion_requires <- function(pkgs, what, call = rlang::caller_env()) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    mudatar_abort(
      "Package{?s} {.pkg {missing}} {?is/are} required to convert to {.cls {what}}",
      call = call
    )
  }
}

check_obs_axis <- function(x, what, call = rlang::caller_env()) {
  if (!x$axis %in% c(0L, -1L)) {
    mudatar_abort_unsupported(
      "Converting a MuData with {.code axis = {x$axis}} (shared variables) to {.cls {what}}",
      "i" = "{.cls {what}} needs shared observations ({.code axis = 0} or {.code -1}).",
      call = call
    )
  }
  if (anyDuplicated(x$obs_names)) {
    mudatar_abort_unsupported("Converting a MuData with duplicated {.field obs_names}", call = call)
  }
}

check_mod_args <- function(mod_args, mod_names, call = rlang::caller_env()) {
  if (is.null(mod_args)) {
    return(list())
  }
  if (!is.list(mod_args) || (length(mod_args) && is.null(names(mod_args)))) {
    mudatar_abort(
      "{.arg mod_args} must be a named list with one list of anndataR arguments per modality",
      call = call
    )
  }
  unknown <- setdiff(names(mod_args), mod_names)
  if (length(unknown)) {
    mudatar_abort(
      c("{.arg mod_args} names unknown modalities {.val {unknown}}", "i" = "Modalities: {.val {mod_names}}"),
      call = call
    )
  }
  for (m in names(mod_args)) {
    args <- mod_args[[m]]
    if (!is.list(args) || (length(args) && is.null(names(args)))) {
      mudatar_abort("{.code mod_args${m}} must be a named list of arguments", call = call)
    }
    bad <- setdiff(names(args), .anndatar_sce_args)
    if (length(bad)) {
      mudatar_abort(
        c(
          "{.code mod_args${m}} has unknown arguments {.val {bad}}",
          "i" = "Allowed (see {.help anndataR::as_SingleCellExperiment}): {.val {(.anndatar_sce_args)}}"
        ),
        call = call
      )
    }
  }
  mod_args
}

#' Convert one modality with anndataR, keeping Awkward obsm data
#'
#' anndataR does not map Awkward Arrays (ragged, e.g. scirpy AIRR data) to
#' reducedDims; they are kept, aligned to the columns, in
#' `metadata(sce)$awkward_obsm`.
#' @noRd
modality_to_sce <- function(adata, args) {
  if (inherits(adata, "ReticulateAnnData")) {
    # read natively: anndataR's conversion of Python-backed objects does not
    # handle Awkward Arrays (nor `layers` with Python anndata >= 0.13)
    adata <- as_in_memory_anndata(adata)
  }
  sce <- do.call(adata$as_SingleCellExperiment, args)
  awkward <- get_awkward(adata, "obsm")
  if (length(awkward)) {
    meta <- S4Vectors::metadata(sce)
    meta$awkward_obsm <- awkward
    S4Vectors::metadata(sce) <- meta
  }
  sce
}

# Global obsm entries other than the per-modality presence masks
global_aligned <- function(m, mod_names) {
  m <- m %||% list()
  m[setdiff(names(m), mod_names)]
}

#' Subset or pad the columns of a converted modality to `idx`
#'
#' `idx` gives, for each output column, the column in `sce` (0 = absent).
#' Absent columns are filled with `NA` (assays, colData, matrix
#' reducedDims); colPairs cannot be padded and are dropped with a warning.
#' @noRd
align_sce_columns <- function(sce, idx, cells, mod) {
  present <- idx > 0
  if (all(present)) {
    out <- sce[, idx]
    colnames(out) <- cells
    out <- align_awkward(out, idx)
    return(out)
  }

  assays <- lapply(SummarizedExperiment::assays(sce), pad_columns, idx = idx)
  cd <- as.data.frame(SummarizedExperiment::colData(sce), optional = TRUE)
  cd <- cd[ifelse(present, idx, NA), , drop = FALSE]
  rownames(cd) <- cells
  rd <- lapply(SingleCellExperiment::reducedDims(sce), function(r) {
    if (!is.matrix(r)) {
      mudatar_abort_unsupported(
        "Padding the non-matrix reducedDim of modality {.val {mod}} with {.code join = \"main\"}"
      )
    }
    out <- r[ifelse(present, idx, NA), , drop = FALSE]
    rownames(out) <- cells
    out
  })
  if (length(SingleCellExperiment::colPairNames(sce))) {
    cli_warn("colPairs of modality {.val {mod}} cannot be padded to the main cells and are dropped")
  }
  meta <- S4Vectors::metadata(sce)
  out <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(cd, check.names = FALSE),
    rowData = SummarizedExperiment::rowData(sce),
    reducedDims = rd,
    rowPairs = SingleCellExperiment::rowPairs(sce),
    metadata = meta
  )
  colnames(out) <- cells
  out$present <- present
  align_awkward(out, idx)
}

# Awkward Arrays cannot hold padding rows: keep the cells that are present and
# record which output columns they belong to
align_awkward <- function(sce, idx) {
  meta <- S4Vectors::metadata(sce)
  if (is.null(meta$awkward_obsm)) {
    return(sce)
  }
  present <- idx > 0
  meta$awkward_obsm <- lapply(meta$awkward_obsm, function(a) a[idx[present]])
  if (!all(present)) {
    meta$awkward_obsm_columns <- colnames(sce)[present]
  }
  S4Vectors::metadata(sce) <- meta
  sce
}

pad_columns <- function(m, idx) {
  present <- idx > 0
  if (inherits(m, "DelayedArray")) {
    m <- if (DelayedArray::is_sparse(m)) methods::as(m, "CsparseMatrix") else as.matrix(m)
  }
  if (is.matrix(m)) {
    out <- m[, ifelse(present, idx, NA), drop = FALSE]
    colnames(out) <- NULL
    return(out)
  }
  if (inherits(m, "sparseMatrix")) {
    out <- Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = c(nrow(m), length(idx)),
      dimnames = list(rownames(m), NULL)
    )
    out[, present] <- m[, idx[present], drop = FALSE]
    out[, !present] <- NA
    return(out)
  }
  mudatar_abort_unsupported("Padding assays of class {.cls {class(m)}}")
}

mudata_to_sce <- function(x, main, mod_args, join, call = rlang::caller_env()) {
  check_conversion_requires(
    c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors"),
    "SingleCellExperiment",
    call = call
  )
  check_obs_axis(x, "SingleCellExperiment", call = call)
  mod_names <- x$mod_names()
  if (is.null(main)) {
    main <- if ("rna" %in% mod_names) "rna" else mod_names[[1]]
  }
  main <- resolve_mod_name(mod_names, main, call = call)
  mod_args <- check_mod_args(mod_args, mod_names, call = call)

  maps <- x$obsmap
  global <- x$obs_names
  keep <- if (join == "intersect") {
    Reduce(`&`, lapply(maps[mod_names], function(m) m > 0L))
  } else {
    maps[[main]] > 0L
  }
  if (!any(keep)) {
    mudatar_abort(
      c(
        "No observations are shared by {.val {if (join == 'intersect') mod_names else main}}",
        "i" = "Use {.code join = \"main\"} to keep every cell of {.val {main}}."
      ),
      call = call
    )
  }
  if (join == "intersect" && sum(keep) < length(global)) {
    cli_inform(
      "Keeping the {sum(keep)} of {length(global)} observations present in every modality ({.code join = \"intersect\"})."
    )
  }
  # the main experiment keeps the main modality's own cell order
  sel <- which(keep)
  sel <- sel[order(maps[[main]][sel])]
  cells <- global[sel]

  sces <- lapply(mod_names, function(m) modality_to_sce(x[[m]], mod_args[[m]] %||% list()))
  names(sces) <- mod_names

  sce <- align_sce_columns(sces[[main]], maps[[main]][sel], cells, main)
  sce$present <- NULL

  # global obs columns join the main colData; clashes are prefixed
  gobs <- x$obs[sel, , drop = FALSE]
  for (col in colnames(gobs)) {
    target <- if (col %in% colnames(SummarizedExperiment::colData(sce))) paste0("mudata:", col) else col
    SummarizedExperiment::colData(sce)[[target]] <- gobs[[col]]
  }
  # global aligned matrices become reducedDims of the main experiment
  for (key in names(global_aligned(x$obsm, mod_names))) {
    value <- x$obsm[[key]]
    if (is.matrix(value) && is.numeric(value) && !key %in% SingleCellExperiment::reducedDimNames(sce)) {
      value <- value[sel, , drop = FALSE]
      rownames(value) <- cells
      SingleCellExperiment::reducedDim(sce, key) <- value
    }
  }

  for (m in setdiff(mod_names, main)) {
    SingleCellExperiment::altExp(sce, m) <- align_sce_columns(sces[[m]], maps[[m]][sel], cells, m)
  }

  meta <- S4Vectors::metadata(sce)
  meta$mudata <- list(
    main = main,
    join = join,
    axis = x$axis,
    uns = x$uns,
    obsm = global_aligned(x$obsm, mod_names) |> lapply(reindex_rows, index_order = sel)
  )
  S4Vectors::metadata(sce) <- meta
  sce
}

mudata_to_mae <- function(x, mod_args, call = rlang::caller_env()) {
  check_conversion_requires(
    c("MultiAssayExperiment", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors"),
    "MultiAssayExperiment",
    call = call
  )
  check_obs_axis(x, "MultiAssayExperiment", call = call)
  mod_names <- x$mod_names()
  mod_args <- check_mod_args(mod_args, mod_names, call = call)
  global <- x$obs_names
  maps <- x$obsmap

  experiments <- lapply(mod_names, function(m) modality_to_sce(x[[m]], mod_args[[m]] %||% list()))
  names(experiments) <- mod_names

  # obsmap is exactly MultiAssayExperiment's sampleMap: one row per
  # (modality, observation) pair linking the global name to the column name
  sample_map <- do.call(rbind, lapply(mod_names, function(m) {
    map <- maps[[m]]
    present <- which(map > 0L)
    data.frame(
      assay = rep(m, length(present)),
      primary = global[present],
      colname = colnames(experiments[[m]])[map[present]]
    )
  }))
  sample_map <- S4Vectors::DataFrame(sample_map)
  sample_map$assay <- factor(sample_map$assay, levels = mod_names)

  col_data <- x$obs
  rownames(col_data) <- global
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(experiments),
    colData = S4Vectors::DataFrame(col_data, check.names = FALSE),
    sampleMap = sample_map,
    metadata = list(
      mudata = list(
        axis = x$axis,
        uns = x$uns,
        obsm = global_aligned(x$obsm, mod_names),
        obsp = x$obsp
      )
    )
  )
}
