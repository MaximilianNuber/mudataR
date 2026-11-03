# Index and map construction compared with Python mudata (oracle cases from
# tools/fixtures/generate_mudata_fixtures.py)

oracle <- jsonlite::fromJSON(
  test_path("fixtures", "expected_update.json"),
  simplifyVector = FALSE
)

build_mods <- function(spec) {
  lapply(spec, function(m) make_adata(unlist(m[[1]]), unlist(m[[2]])))
}

expect_matches_snapshot <- function(mdata, snap) {
  expect_identical(mdata$obs_names, as.character(unlist(snap$obs_names)))
  expect_identical(mdata$var_names, as.character(unlist(snap$var_names)))
  expect_identical(names(mdata$obsmap), names(snap$obsmap))
  for (m in names(snap$obsmap)) {
    expect_identical(mdata$obsmap[[m]], as.integer(unlist(snap$obsmap[[m]])), label = paste("obsmap", m))
  }
  for (m in names(snap$varmap)) {
    expect_identical(mdata$varmap[[m]], as.integer(unlist(snap$varmap[[m]])), label = paste("varmap", m))
  }
  for (col in names(snap$obs_columns)) {
    expected <- vapply(snap$obs_columns[[col]], function(v) if (is.null(v)) NA_real_ else as.numeric(v), numeric(1))
    expect_equal(mdata$obs[[col]], expected, label = paste("obs column", col))
  }
}

for (case in oracle$cases) {
  test_that(paste("MuData indices match Python for", case$name), {
    mods <- build_mods(case$mods)
    if (case$name == "axis0_shared_feature_names") {
      # Python also warns that var_names are not unique
      expect_warning(mdata <- MuData(mods, axis = case$axis), "not unique")
    } else {
      mdata <- MuData(mods, axis = case$axis)
    }
    expect_equal(mdata$axis, case$axis)
    expect_matches_snapshot(mdata, case$initial)

    if (!is.null(case$final)) {
      for (col in names(case$obs_column)) {
        obs <- mdata$obs
        obs[[col]] <- as.numeric(unlist(case$obs_column[[col]]))
        mdata$obs <- obs
      }
      for (m in names(case$replace)) {
        mdata[[m]] <- make_adata(unlist(case$replace[[m]][[1]]), unlist(case$replace[[m]][[2]]))
      }
      expect_matches_snapshot(mdata, case$final)
    }
  })
}
