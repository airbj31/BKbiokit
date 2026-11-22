ref_of <- function(x) levels(.pick_reference(x, "group")$grp)[1]

test_that("factor: first level is the reference", {
  x <- factor(c("B", "A", "B", "A"), levels = c("B", "A", "C"))
  expect_equal(ref_of(x), "B")
})

test_that("character: common control words are detected", {
  expect_equal(ref_of(c("KO", "KO", "WT", "WT")), "WT")
  expect_equal(ref_of(c("treat", "DMSO_vehicle")), "DMSO_vehicle")
  expect_equal(ref_of(c("Ctrl-24h", "drug-24h")), "Ctrl-24h")
  expect_equal(ref_of(c("siTP53", "siNC")), "siNC")
})

test_that("short words do not match inside other words", {
  expect_error(.pick_reference(c("Concentration", "antiPD1"), "group"),
               "factor")
})

test_that("stops when both groups match", {
  expect_error(.pick_reference(c("WT_treated", "WT_vehicle"), "group"),
               "WT_treated, WT_vehicle")
})

test_that("stops when there are not exactly two groups", {
  expect_error(.pick_reference(c("a", "b", "c"), "group"), "2")
  expect_error(.pick_reference(c("a", "a"), "group"), "2")
})

test_that("run_deseq checks meta columns and samples before fitting", {
  cnt  <- matrix(1:8, nrow = 2, dimnames = list(c("g1", "g2"), paste0("s", 1:4)))
  meta <- data.frame(sample_id = paste0("s", 1:4), group = c("WT", "WT", "KO", "KO"))

  expect_error(run_deseq(cnt, meta, "sample_id", "grp"), "grp")
  expect_error(run_deseq(cnt, meta[1:3, ], "sample_id", "group"), "s4")
})

test_that("run_deseq returns the expected structure", {
  skip_if_not_installed("DESeq2")
  set.seed(1)
  cnt <- matrix(rnbinom(400, mu = 100, size = 5), nrow = 50,
                dimnames = list(paste0("g", 1:50), paste0("s", 1:8)))
  meta <- data.frame(sample_id = paste0("s", 1:8),
                     group = rep(c("vehicle", "drug"), each = 4))

  res <- suppressWarnings(run_deseq(cnt, meta, "sample_id", "group", verbose = FALSE))

  expect_named(res, c("dds", "norm_cnt", "deseq", "comparison"))
  expect_equal(res$comparison$control, "vehicle")
  expect_equal(res$comparison$case, "drug")
  expect_true("gene_id" %in% names(res$deseq))
})
