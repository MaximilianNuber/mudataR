skip_if_no_hdf5()

h5mu <- fixture("gex_airr.h5mu")
expected <- scirpy_expected()

copy_fixture <- function(env = parent.frame()) {
  path <- withr::local_tempfile(fileext = ".h5mu", .local_envir = env)
  file.copy(h5mu, path)
  path
}

test_that("a Python-written H5MU is opened natively and lazily", {
  m <- read_h5mu(h5mu, as = "HDF5MuData")
  expect_true(is_HDF5MuData(m))
  expect_true(m$isbacked)
  expect_equal(m$filename, normalizePath(h5mu))
  expect_equal(m$axis, 0L)
  expect_equal(m$mod_names(), expected$mudata$mod_order |> unlist())
  expect_equal(m$obs_names, unlist(expected$mudata$obs_names))
  expect_equal(m$var_names, unlist(expected$mudata$var_names))
  for (mod in names(expected$mudata$obsmap)) {
    expect_equal(m$obsmap[[mod]], as.integer(unlist(expected$mudata$obsmap[[mod]])))
  }
  expect_equal(names(m$obsmap), c("gex", "airr"))
  expect_equal(m$uns$fixture, "scirpy-gex-airr")
})

test_that("modalities are read-only anndataR HDF5AnnData objects", {
  m <- read_h5mu(h5mu, as = "HDF5MuData")
  gex <- m[["gex"]]
  expect_s3_class(gex, "HDF5AnnData")
  expect_equal(gex$shape(), c(7L, 5L))
  expect_s4_class(gex$X, "dgRMatrix")
  expect_equal(levels(gex$obs$sample), c("s2", "s1"))
})

test_that("backed mode returns DelayedArray matrices and lazy Awkward buffers", {
  skip_if_not_installed("HDF5Array")
  m <- read_h5mu(h5mu, as = "HDF5MuData", backed = TRUE)
  gex <- m[["gex"]]
  expect_s4_class(gex$X, "DelayedArray")
  airr <- get_awkward(m[["airr"]])$airr
  expect_null(awk_field(airr, "state")$buffers)
  expect_length(as.data.frame(airr, unnest = TRUE)$locus, 13)
})

test_that("closing the parent or a modality never affects the others", {
  m <- read_h5mu(h5mu, as = "HDF5MuData")
  gex <- m[["gex"]]
  airr <- m[["airr"]]
  gex$close()
  expect_equal(airr$n_obs(), 7L)
  expect_equal(m$n_obs(), 9L)
  m$close()
  expect_false(m$is_open())
  expect_error(m$obs, class = "mudatar_closed_error")
  expect_error(m[["gex"]], class = "mudatar_closed_error")
  # modality objects obtained earlier hold their own handles
  expect_equal(airr$n_obs(), 7L)
  m$open()
  expect_equal(m$n_obs(), 9L)
})

test_that("modalities of an r+ H5MU are read-only; replacing them works", {
  path <- copy_fixture()
  m <- read_h5mu(path, as = "HDF5MuData", mode = "r+")
  gex <- m[["gex"]]
  obs <- gex$obs
  obs$extra <- 1
  expect_error(gex$obs <- obs, "read-only")
  mem <- gex$as_InMemoryAnnData()
  mem$obs <- obs
  m[["gex"]] <- mem
  expect_equal(m[["gex"]]$obs$extra, rep(1, 7))
  # the other modality and its Awkward Arrays are untouched
  expect_length(get_awkward(m[["airr"]])$airr, 7)
  # re-assigning a modality of the same store keeps its Awkward Arrays
  m[["airr"]] <- m[["airr"]]
  expect_length(get_awkward(m[["airr"]])$airr, 7)
})

test_that("read-only objects reject modification", {
  m <- read_h5mu(h5mu, as = "HDF5MuData")
  expect_error(m$uns <- list(a = 1), class = "mudatar_readonly_error")
  expect_error(m[["new"]] <- make_adata("c0", "x"), class = "mudatar_readonly_error")
  expect_error(m$axis <- 1, "read-only")
})

test_that("selected modalities limit the mapping without touching others", {
  m <- read_h5mu(h5mu, as = "HDF5MuData", modalities = "airr")
  expect_equal(m$mod_names(), "airr")
  expect_equal(names(m$obsmap), "airr")
  expect_error(m[["gex"]], "No modality")
  expect_error(read_h5mu(h5mu, as = "HDF5MuData", modalities = "nope"), "not in")
  mem <- read_h5mu(h5mu, modalities = "airr")
  expect_equal(mem$mod_names(), "airr")
})

test_that("r+ objects write global slots and modalities in place", {
  path <- copy_fixture()
  m <- read_h5mu(path, as = "HDF5MuData", mode = "r+")
  m$uns <- list(note = "edited in R")
  expect_equal(m$uns$note, "edited in R")
  obs <- m$obs
  obs$batch <- factor(rep(c("x", "y", "z"), 3))
  m$obs <- obs
  expect_equal(levels(m$obs$batch), c("x", "y", "z"))
  m[["prot"]] <- make_adata(c("c1", "zz"), c("p1", "p2"))
  expect_equal(m$mod_names(), c("gex", "airr", "prot"))
  expect_equal(m$obs_names, c(unlist(expected$mudata$obs_names), "zz"))
  # existing global order is kept, the new observation is appended
  expect_equal(m$obsmap$prot, c(0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L))
  expect_s3_class(m[["prot"]], "HDF5AnnData")
  # re-reading from disk gives the same state
  again <- read_h5mu(path)
  expect_equal(again$mod_names(), c("gex", "airr", "prot"))
  expect_equal(again$obs$batch[1:9], obs$batch)
  expect_true(is.na(again$obs$batch[10]))
})

test_that("non-MuData files are rejected", {
  skip_if_not_installed("anndataR")
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  anndataR::write_h5ad(make_adata(c("a", "b"), "g"), h5ad)
  expect_error(read_h5mu(h5ad, as = "HDF5MuData"), "not a MuData store") |>
    expect_warning("no MuData user block")
  expect_error(read_h5mu(tempfile(fileext = ".h5mu")), "does not exist")
})
