skip_if_no_hdf5()

h5mu <- fixture("gex_airr.h5mu")

example_mudata <- function() {
  rna <- anndataR::AnnData(
    X = matrix(as.numeric(1:6), 3, 2),
    obs = data.frame(row.names = c("a", "b", "c"), ct = factor(c("T", "B", "T"))),
    var = data.frame(row.names = c("g1", "g2"))
  )
  prot <- make_adata(c("b", "c", "d"), "p1")
  m <- MuData(
    list(rna = rna, prot = prot),
    uns = list(note = "x", n = 3L)
  )
  m$obsm <- c(m$obsm, list(emb = matrix(as.numeric(1:8), 4, 2)))
  m
}

test_that("R writes H5MU files with the MuData layout and user block", {
  path <- withr::local_tempfile(fileext = ".h5mu")
  m <- example_mudata()
  write_h5mu(m, path)

  marker <- mudataR:::h5mu_read_marker(path)
  expect_match(marker, "^MuData \\(format-version=0\\.1\\.0;creator=mudataR;creator-version=")
  fid <- rhdf5::H5Fopen(path, "H5F_ACC_RDONLY")
  fcpl <- rhdf5::H5Fget_create_plist(fid)
  expect_equal(rhdf5::H5Pget_userblock(fcpl), 512)
  rhdf5::H5Pclose(fcpl)
  rhdf5::H5Fclose(fid)

  attrs <- rhdf5::h5readAttributes(path, "/")
  expect_equal(attrs[["encoding-type"]], "MuData")
  expect_equal(attrs[["encoding-version"]], "0.1.0")
  expect_equal(as.integer(attrs[["axis"]]), 0L)
  expect_equal(attrs[["encoder"]], "mudataR")
  expect_equal(as.vector(rhdf5::h5readAttributes(path, "/mod")[["mod-order"]]), c("rna", "prot"))
  expect_equal(rhdf5::h5readAttributes(path, "/mod/rna")[["encoding-type"]], "anndata")

  # maps are unsigned 32-bit, as written by Python mudata
  contents <- rhdf5::h5ls(path, all = TRUE)
  map_types <- contents$dtype[contents$group %in% c("/obsmap", "/varmap")]
  expect_equal(unique(map_types), "H5T_STD_U32LE")
})

test_that("R-written H5MU files round trip in R", {
  path <- withr::local_tempfile(fileext = ".h5mu")
  m <- example_mudata()
  write_h5mu(m, path)
  back <- read_h5mu(path)
  expect_equal(back$obs_names, m$obs_names)
  expect_equal(back$var_names, m$var_names)
  expect_equal(back$obsmap, m$obsmap)
  expect_equal(back$varmap, m$varmap)
  expect_equal(back$mod_names(), m$mod_names())
  expect_equal(back$uns$note, "x")
  expect_equal(unname(back$obsm$emb), unname(m$obsm$emb), ignore_attr = TRUE)
  expect_equal(unname(as.matrix(back[["rna"]]$X)), unname(m[["rna"]]$X), ignore_attr = TRUE)
  expect_equal(back[["rna"]]$obs$ct, m[["rna"]]$obs$ct, ignore_attr = TRUE)
})

test_that("Python-written H5MU files round trip through R", {
  path <- withr::local_tempfile(fileext = ".h5mu")
  src <- read_h5mu(h5mu, as = "HDF5MuData")
  write_h5mu(src, path)
  back <- read_h5mu(path, as = "HDF5MuData")
  expect_equal(back$obs_names, src$obs_names)
  expect_equal(back$obsmap, src$obsmap)
  expect_equal(back$mod_names(), src$mod_names())
})

test_that("writes never clobber an existing file unless asked", {
  path <- withr::local_tempfile(fileext = ".h5mu")
  m <- example_mudata()
  write_h5mu(m, path)
  before <- tools::md5sum(path)
  expect_error(write_h5mu(m, path), "already exists")
  expect_equal(tools::md5sum(path), before)
  write_h5mu(m, path, mode = "w")
  expect_true(file.exists(path))
  expect_length(list.files(dirname(path), pattern = paste0("^\\.", basename(path))), 0)
})

test_that("an invalid object fails before the destination is touched", {
  path <- withr::local_tempfile(fileext = ".h5mu")
  m <- example_mudata()
  write_h5mu(m, path)
  before <- tools::md5sum(path)
  # rename observations through the modality's reference semantics, which
  # bypasses update(): the maps are now inconsistent with the modalities
  rna <- m[["rna"]]
  rna$obs_names <- c("x", "y", "z")
  expect_error(write_h5mu(m, path, mode = "w"), "does not match")
  expect_equal(tools::md5sum(path), before)
  expect_length(list.files(dirname(path), pattern = paste0("^\\.", basename(path))), 0)
})

test_that("Zarr v2 and v3 MuData stores round trip in R", {
  skip_if_no_zarr()
  m <- example_mudata()
  for (zarr_format in c(2L, 3L)) {
    path <- withr::local_tempfile(fileext = ".zarr")
    write_mudata_zarr(m, path, zarr_format = zarr_format)
    z <- read_mudata_zarr(path, as = "ZarrMuData")
    expect_equal(z$zarr_format(), zarr_format)
    expect_equal(z$mod_names(), c("rna", "prot"))
    expect_equal(z$obsmap, m$obsmap)
    attrs <- Rarr::read_zarr_attributes(file.path(path, "mod"))
    expect_equal(unlist(attrs[["mod-order"]]), c("rna", "prot"))
    # a single modality still records mod-order as an array
    one <- withr::local_tempfile(fileext = ".zarr")
    write_mudata_zarr(MuData(list(rna = m[["rna"]])), one, zarr_format = zarr_format)
    expect_type(Rarr::read_zarr_attributes(file.path(one, "mod"))[["mod-order"]], "list")
  }
})
