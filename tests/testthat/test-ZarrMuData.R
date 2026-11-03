skip_if_no_zarr()

expected <- scirpy_expected()

for (zarr_format in c(2L, 3L)) {
  store <- fixture_zarr(sprintf("gex_airr_v%d.zarr", zarr_format))

  test_that(sprintf("a Python-written MuData Zarr v%d store is read natively", zarr_format), {
    m <- read_mudata_zarr(store, as = "ZarrMuData")
    expect_true(is_ZarrMuData(m))
    expect_equal(m$zarr_format(), zarr_format)
    expect_equal(m$mod_names(), c("gex", "airr"))
    expect_equal(m$obs_names, unlist(expected$mudata$obs_names))
    expect_equal(m$obsmap$airr, as.integer(unlist(expected$mudata$obsmap$airr)))
    airr <- m[["airr"]]
    expect_s3_class(airr, "ZarrAnnData")
    expect_equal(airr$n_obs(), 7L)
    expect_length(get_awkward(airr)$airr, 7)
    m$close()
    expect_error(m$obs, class = "mudatar_closed_error")
  })

  test_that(sprintf("Zarr v%d modalities without Awkward Arrays are edited in place", zarr_format), {
    path <- withr::local_tempfile(fileext = ".zarr")
    write_mudata_zarr(read_mudata_zarr(store), path, zarr_format = zarr_format)
    m <- read_mudata_zarr(path, as = "ZarrMuData", mode = "r+")
    gex <- m[["gex"]]
    obs <- gex$obs
    obs$extra <- seq_len(nrow(obs))
    gex$obs <- obs
    expect_equal(read_mudata_zarr(path)[["gex"]]$obs$extra, seq_len(7))
    # modalities with Awkward Arrays are opened read-only
    airr <- m[["airr"]]
    expect_error(airr$obs <- airr$obs)
  })

  test_that(sprintf("MuData Zarr v%d stores are read into memory", zarr_format), {
    mem <- read_mudata_zarr(store)
    expect_true(is_InMemoryMuData(mem))
    expect_s3_class(mem[["gex"]], "InMemoryAnnData")
    expect_equal(dim(mem), c(9L, 5L))
  })
}
