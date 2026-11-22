ref_of <- function(x) levels(.pick_reference(x, "group")$grp)[1]

test_that("factor: first level is the reference", {
  x <- factor(c("B", "A", "B", "A"), levels = c("B", "A", "C"))
  expect_equal(ref_of(x), "B")
})

test_that("factor: warns when a control-looking level is not the reference", {
  x <- factor(c("vehicle", "drug"))          # alphabetical -> drug first
  expect_warning(r <- ref_of(x), "relevel")
  expect_equal(r, "drug")                     # factor order is still respected
  expect_no_warning(ref_of(factor(c("vehicle", "drug"), levels = c("vehicle", "drug"))))
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
  meta <- data.frame(group = c("WT", "WT", "KO", "KO"), row.names = paste0("s", 1:4))

  expect_error(run_deseq(cnt, meta, "grp"), "grp")
  expect_error(run_deseq(cnt, meta[1:3, , drop = FALSE], "group"), "s4")
  # tibble / no row names -> helpful error
  expect_error(run_deseq(cnt, data.frame(group = meta$group), "group"), "column_to_rownames")
})

test_that("run_deseq reports reordered and extra meta samples", {
  # messages are emitted before DESeq2 runs, so try() is enough here
  cnt  <- matrix(1:8, nrow = 2, dimnames = list(c("g1", "g2"), paste0("s", 1:4)))
  meta <- data.frame(group = c("KO", "KO", "WT", "WT", "WT"),
                     row.names = c("s4", "s3", "s2", "s1", "s5"))
  msgs <- character()
  withCallingHandlers(
    try(run_deseq(cnt, meta, "group"), silent = TRUE),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )
  expect_true(any(grepl("s5", msgs)))
  expect_true(any(grepl("reordered", msgs)))
})

test_that("run_deseq returns the expected structure", {
  skip_if_not_installed("DESeq2")
  set.seed(1)
  cnt <- matrix(rnbinom(400, mu = 100, size = 5), nrow = 50,
                dimnames = list(paste0("g", 1:50), paste0("s", 1:8)))
  meta <- data.frame(group = rep(c("vehicle", "drug"), each = 4),
                     row.names = paste0("s", 8:1))  # different order on purpose

  res <- suppressWarnings(run_deseq(cnt, meta, "group", verbose = FALSE))

  expect_named(res, c("dds", "norm_cnt", "size_factors", "deseq", "comparison"))
  expect_equal(names(res$size_factors), colnames(cnt))
  expect_equal(res$comparison$control, "vehicle")
  expect_equal(res$comparison$case, "drug")
  expect_true("gene_id" %in% names(res$deseq))
  expect_equal(colnames(res$norm_cnt), colnames(cnt))
})

fake_result <- function() {
  structure(list(
    dds = NULL,
    norm_cnt = matrix(1, nrow = 40, ncol = 8),
    size_factors = setNames(c(0.8, 1, 1.2, 0.9, 1.1, 1, 1, 1.05), paste0("s", 1:8)),
    deseq = data.frame(gene_id = paste0("g", 1:6),
                       log2FoldChange = c(2, -1, 3, 0.5, -4, 1.5),
                       padj = c(0.001, 0.02, NA, 0.5, 0.08, 0.0001)),
    comparison = list(case = "drug", control = "vehicle", design = "~group",
                      n = c(vehicle = 4, drug = 4),
                      genes = c(input = 50, tested = 40), cnt_thres = 10,
                      id = "gene_id")
  ), class = c("run_deseq_result", "list"))
}

test_that("print shows a summary and top genes, not the full tables", {
  out <- capture.output(print(fake_result()))
  expect_true(any(grepl("40 genes x 8 samples", out)))
  expect_true(any(grepl("2 up / 1 down", out)))
  expect_true(any(grepl("0.80 - 1.20", out)))
  # top genes are sorted by padj: g6 (1e-4) before g1 (1e-3)
  expect_lt(grep("g6", out)[1], grep("g1", out)[1])
  expect_false(any(grepl("\\bg4\\b", out)))   # not significant -> hidden

  out0 <- capture.output(print(fake_result(), n = 0))
  expect_false(any(grepl("Top", out0)))
})

test_that("summary counts DE genes at several cutoffs", {
  s <- summary(fake_result())
  expect_s3_class(s, "summary.run_deseq_result")
  expect_equal(s$de_table$total, c(2, 3, 4))      # padj < 0.01, 0.05, 0.1
  expect_equal(s$n_na_padj, 1)
  expect_equal(summary(fake_result(), lfc = 1.8)$de_table$total, c(1, 1, 2))
  expect_output(print(s), "DESeq2 summary")
})

test_that("as.data.frame and as_tibble return the result table", {
  x <- fake_result()
  expect_equal(as.data.frame(x), x$deseq)
  expect_s3_class(tibble::as_tibble(x), "tbl_df")
  expect_equal(nrow(tibble::as_tibble(x)), 6)
})
