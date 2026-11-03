test_that("MuData() builds an InMemoryMuData of anndataR modalities", {
  mdata <- MuData(list(
    rna = make_adata(c("a", "b", "c"), c("g1", "g2")),
    prot = make_adata(c("b", "c", "d"), "p1")
  ))
  expect_true(is_mudata(mdata))
  expect_true(is_InMemoryMuData(mdata))
  expect_s3_class(mdata[["rna"]], "InMemoryAnnData")
  expect_identical(mdata$mod$prot, mdata[["prot"]])
  expect_equal(dim(mdata), c(4L, 3L))
  expect_equal(nrow(mdata), 4L)
  expect_equal(rownames(mdata), c("a", "b", "c", "d"))
  expect_equal(colnames(mdata), c("g1", "g2", "p1"))
  expect_equal(names(mdata$mod), c("rna", "prot"))
  expect_length(mdata$mod, 2)
  expect_equal(mdata$n_mod(), 2L)
  # Python records modality membership in obsm/varm
  expect_equal(mdata$obsm$prot, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(mdata$varm$rna, c(TRUE, TRUE, FALSE))
  expect_output(print(mdata), "MuData object with n_obs . n_vars = 4 . 3")
})

test_that("modality order follows insertion order", {
  mdata <- MuData(list(z = make_adata("a", "x"), a = make_adata("a", "y")))
  expect_equal(mdata$mod_names(), c("z", "a"))
  expect_equal(mdata$var_names, c("x", "y"))
})

test_that("invalid modality collections are rejected", {
  expect_error(MuData(list()), "non-empty")
  expect_error(MuData(list(make_adata("a", "x"))), "named")
  expect_error(MuData(list(`a/b` = make_adata("a", "x"))), "must not contain")
  expect_error(MuData(list(x = data.frame())), "anndataR")
  ad <- make_adata("a", "x")
  expect_error(MuData(stats::setNames(list(ad, ad), c("m", "m"))), "unique")
  expect_error(MuData(list(m = ad), axis = 2), "axis")
})

test_that("supplied maps are validated, not silently recomputed", {
  mods <- list(rna = make_adata(c("a", "b"), "g"), prot = make_adata("b", "p"))
  ok <- MuData(
    mods,
    obs_names = c("b", "a"),
    obsmap = list(rna = c(2L, 1L), prot = c(1L, 0L)),
    var_names = c("g", "p"),
    varmap = list(rna = c(1L, 0L), prot = c(0L, 1L))
  )
  expect_equal(ok$obs_names, c("b", "a"))
  expect_error(
    MuData(mods, obs_names = c("a", "b"), obsmap = list(rna = c(1L, 3L), prot = c(0L, 1L))),
    "points beyond"
  )
  expect_error(
    MuData(mods, obs_names = c("a", "b"), obsmap = list(rna = c(2L, 1L), prot = c(0L, 1L))),
    "does not match"
  )
  expect_error(
    MuData(mods, obs_names = c("a", "b"), obsmap = list(rna = c(1L, 2L))),
    "missing modalities"
  )
})

test_that("a supplied global order is respected", {
  mdata <- MuData(
    list(rna = make_adata(c("a", "b", "c"), "g")),
    obs = data.frame(row.names = c("c", "a", "b"), score = 3:1)
  )
  expect_equal(mdata$obs_names, c("c", "a", "b"))
  expect_equal(mdata$obsmap$rna, c(3L, 1L, 2L))
  expect_equal(mdata$obs$score, 3:1)
  expect_equal(rownames(mdata$obs), c("c", "a", "b"))
})

test_that("modalities can be added, replaced and removed", {
  mdata <- MuData(list(rna = make_adata(c("a", "b"), "g")))
  mdata[["prot"]] <- make_adata(c("b", "c"), "p")
  expect_equal(mdata$obs_names, c("a", "b", "c"))
  expect_equal(mdata$obsmap$prot, c(0L, 1L, 2L))
  mdata$mod[["prot"]] <- NULL
  expect_equal(mdata$mod_names(), "rna")
  expect_null(mdata$obsm$prot)
  expect_equal(mdata$obsmap, list(rna = c(1L, 2L)))
  expect_error(mdata[["missing"]], "No modality")
  expect_error(mdata$mod <- list(), "Assign individual")
})

test_that("global slots are validated", {
  mdata <- MuData(list(rna = make_adata(c("a", "b"), c("g1", "g2"))))
  expect_error(mdata$obs <- data.frame(x = 1:3), "2 rows")
  expect_error(mdata$obsm <- list(x = matrix(1, 3, 2)), "aligned")
  expect_error(mdata$obsp <- list(x = matrix(1, 2, 3)), "aligned")
  mdata$obsm <- list(x = matrix(1:4, 2, 2), rna = c(TRUE, TRUE))
  expect_equal(mdata$obsm_keys(), c("x", "rna"))
  expect_error(mdata$axis <- 1, "read-only")
  expect_error(mdata$isbacked <- TRUE, "read-only")
  expect_false(mdata$isbacked)
  expect_null(mdata$filename)
  expect_error(mdata[1, ], class = "mudatar_unsupported")
})

test_that("duplicated names within a modality are reported, not approximated", {
  # anndataR cannot hold duplicated obs_names yet (scverse/anndataR#256), so
  # exercise the index builder directly
  expect_error(
    mudataR:::compute_global_attr(list(rna = c("a", "a")), "obs", 0L),
    class = "mudatar_unsupported"
  )
})
