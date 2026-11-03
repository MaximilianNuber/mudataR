skip_if_no_hdf5()
skip_if_no_awkward()

expected <- scirpy_expected()
h5mu <- fixture("gex_airr.h5mu")

py_fixture <- function() {
  skip_if_no_reticulate_mudata()
  md <- reticulate::import("mudata", convert = FALSE)
  quietly(md$read_h5mu(h5mu))
}

test_that("a Python MuData is wrapped with ReticulateAnnData modalities", {
  py <- py_fixture()
  r <- quietly(ReticulateMuData$new(py))
  expect_true(is_ReticulateMuData(r))
  expect_true(is_mudata(r))
  expect_identical(as_ReticulateMuData(r), r)
  expect_identical(r$py_mudata(), py)
  expect_equal(r$mod_names(), c("gex", "airr"))
  expect_equal(dim(r), c(9L, 5L))
  expect_equal(r$axis, 0L)
  expect_equal(r$obs_names, unlist(expected$mudata$obs_names))
  native <- read_h5mu(h5mu, as = "HDF5MuData")
  expect_identical(r$obsmap, native$obsmap)
  expect_identical(r$varmap, native$varmap)
  for (mod in as.list(r$mod)) expect_s3_class(mod, "ReticulateAnnData")
  airr <- quietly(r[["airr"]]$obsm$airr)
  expect_s3_class(airr, "AwkwardArrowArray")
  expect_equal(
    as.data.frame(airr, unnest = TRUE),
    as.data.frame(get_awkward(native[["airr"]])$airr, unnest = TRUE)
  )
  expect_equal(r$uns$fixture, "scirpy-gex-airr")
  expect_output(print(r), "MuData object with n_obs")
})

test_that("Python update(), pull_obs() and copy() are used", {
  py <- py_fixture()
  r <- quietly(ReticulateMuData$new(py))
  cp <- r$copy()
  r$pull_obs(columns = "sample", mods = "gex")
  expect_equal(r$obs_keys(), "gex:sample")
  expect_length(cp$obs_keys(), 0)
  r[["prot"]] <- make_adata(c("c1", "zz"), c("p1", "p2"))
  expect_equal(r$mod_names(), c("gex", "airr", "prot"))
  expect_s3_class(r[["prot"]], "ReticulateAnnData")
  expect_equal(tail(r$obs_names, 1), "zz")
  r[["prot"]] <- NULL
  expect_equal(r$mod_names(), c("gex", "airr"))
  expect_error(r$obsmap <- list(), class = "mudatar_readonly_error")
  expect_s3_class(r$py_call("copy"), "ReticulateMuData")
})

test_that("native MuData objects convert to and from ReticulateMuData", {
  skip_if_no_reticulate_mudata()
  mem <- read_h5mu(h5mu)
  mem$uns <- list(note = "from R")
  r <- quietly(mem$as_ReticulateMuData())
  expect_s3_class(r, "ReticulateMuData")
  for (mod in as.list(r$mod)) expect_s3_class(mod, "ReticulateAnnData")
  # Python computes the same indices as mudataR
  expect_identical(r$obs_names, mem$obs_names)
  expect_identical(r$obsmap, mem$obsmap)
  expect_equal(r$uns$note, "from R")
  expect_s3_class(quietly(r[["airr"]]$obsm$airr), "AwkwardArrowArray")

  back <- quietly(r$as_InMemoryMuData())
  for (mod in as.list(back$mod)) expect_s3_class(mod, "InMemoryAnnData")
  expect_identical(back$obsmap, mem$obsmap)
  expect_equal(
    as.data.frame(as_AwkwardArrowArray(back[["airr"]]$obsm$airr), unnest = TRUE),
    as.data.frame(as_AwkwardArrowArray(mem[["airr"]]$obsm$airr), unnest = TRUE)
  )

  # a ReticulateMuData is written natively and converts like any MuData
  out <- withr::local_tempfile(fileext = ".h5mu")
  quietly(write_h5mu(r, out))
  expect_equal(read_h5mu(out)$obsmap, mem$obsmap)
  skip_if_not_installed("SingleCellExperiment")
  sce <- quietly(suppressMessages(r$as_SingleCellExperiment(main = "gex")))
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(SingleCellExperiment::altExpNames(sce), "airr")
})

test_that("ReticulateMuData rejects non-MuData input", {
  skip_if_no_reticulate_mudata()
  expect_error(ReticulateMuData$new(42), "mudata.MuData")
})
