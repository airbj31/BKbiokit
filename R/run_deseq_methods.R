# Methods for run_deseq_result objects


# ---- helpers ----------------------------------------------------------------

#' @noRd
.is_sig <- function(res, padj, lfc = 0) {
  !is.na(res$padj) & res$padj < padj & abs(res$log2FoldChange) > lfc
}

#' @noRd
.id_col <- function(x) {
  id <- x$comparison$id
  if (is.null(id) || !id %in% names(x$deseq)) names(x$deseq)[1] else id
}

#' @noRd
.field <- function(name, value) {
  cli::cat_line("  ", cli::col_grey(format(name, width = 12)), value)
}


# ---- print ------------------------------------------------------------------

#' Print a run_deseq result
#'
#' Shows a short summary and the top DE genes instead of the full tables.
#' The data are still in `x$deseq`, `x$norm_cnt`, `x$size_factors`, `x$dds`
#' and `x$comparison`.
#'
#' @param x A `run_deseq_result` returned by [run_deseq()].
#' @param padj Adjusted p-value cutoff used to count DE genes.
#' @param n Number of top DE genes (by padj) to show. Use 0 to hide them.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.run_deseq_result <- function(x, padj = 0.05, n = 5, ...) {
  cmp <- x$comparison
  res <- x$deseq
  sf  <- x$size_factors
  sig <- .is_sig(res, padj)
  up   <- sum(sig & res$log2FoldChange > 0)
  down <- sum(sig & res$log2FoldChange < 0)

  cli::cat_rule(left = cli::style_bold("DESeq2 result"),
                right = paste(cmp$case, "vs", cmp$control))

  .field("comparison", paste0(cli::col_red(cmp$case), " vs ",
                              cli::col_blue(cmp$control), " (reference)"))
  .field("design", cmp$design)
  .field("samples", paste0(sum(cmp$n), " (", cmp$case, " n = ", cmp$n[[cmp$case]],
                           ", ", cmp$control, " n = ", cmp$n[[cmp$control]], ")"))
  .field("genes", paste0(cmp$genes[["tested"]], " tested / ", cmp$genes[["input"]],
                         " input (total count >= ", cmp$cnt_thres, ")"))
  .field("norm_cnt", paste0(nrow(x$norm_cnt), " genes x ", ncol(x$norm_cnt), " samples"))
  if (!is.null(sf)) {
    .field("size factor", sprintf("%.2f - %.2f (median %.2f)",
                                  min(sf), max(sf), stats::median(sf)))
  }
  .field("DE genes", paste0(cli::col_red(paste(up, "up")), " / ",
                            cli::col_blue(paste(down, "down")),
                            " in ", cmp$case, " (padj < ", padj, ")"))

  n_show <- min(n, up + down)
  if (n_show > 0) {
    top <- res[sig, , drop = FALSE]
    top <- top[order(top$padj), , drop = FALSE][seq_len(n_show), , drop = FALSE]
    tbl <- data.frame(
      gene   = top[[.id_col(x)]],
      log2FC = sprintf("%+.2f", top$log2FoldChange),
      padj   = formatC(top$padj, format = "e", digits = 2)
    )
    cli::cat_line()
    cli::cat_line(cli::style_bold(paste0("  Top ", n_show, " DE genes (by padj)")))
    out <- utils::capture.output(print(tbl, row.names = FALSE, right = FALSE))
    cli::cat_line(paste0("  ", out))
  }

  cli::cat_line()
  cli::cat_line(cli::col_grey(
    "  $deseq $norm_cnt $size_factors $dds $comparison | summary() as_tibble()"))
  cli::cat_rule()
  invisible(x)
}


# ---- summary ----------------------------------------------------------------

#' Summarise a run_deseq result
#'
#' Counts up- and down-regulated genes at several adjusted p-value cutoffs,
#' and describes the log2 fold change and size factor distributions.
#'
#' @param object A `run_deseq_result` returned by [run_deseq()].
#' @param cutoffs Adjusted p-value cutoffs to count DE genes at.
#' @param lfc Minimum absolute log2 fold change for a gene to count as DE.
#' @param ... Ignored.
#' @return A `summary.run_deseq_result` object (a list) with `de_table`,
#'   `lfc_quantiles`, `size_factors`, `n_tested`, `n_na_padj`, `lfc` and
#'   `comparison`.
#' @export
summary.run_deseq_result <- function(object, cutoffs = c(0.01, 0.05, 0.1),
                                     lfc = 0, ...) {
  res <- object$deseq

  de_table <- do.call(rbind, lapply(cutoffs, function(p) {
    sig <- .is_sig(res, p, lfc)
    data.frame(padj = p,
               up   = sum(sig & res$log2FoldChange > 0),
               down = sum(sig & res$log2FoldChange < 0))
  }))
  de_table$total <- de_table$up + de_table$down

  out <- list(
    comparison    = object$comparison,
    de_table      = de_table,
    lfc_quantiles = stats::quantile(res$log2FoldChange,
                                    c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
    size_factors  = if (!is.null(object$size_factors)) summary(object$size_factors),
    n_tested      = nrow(res),
    n_na_padj     = sum(is.na(res$padj)),
    lfc           = lfc
  )
  class(out) <- "summary.run_deseq_result"
  out
}

#' @export
print.summary.run_deseq_result <- function(x, ...) {
  cmp <- x$comparison

  cli::cat_rule(left = cli::style_bold("DESeq2 summary"),
                right = paste(cmp$case, "vs", cmp$control))

  cli::cat_line()
  cli::cat_line(cli::style_bold("  DE genes"),
                if (x$lfc > 0) paste0(" (|log2FC| > ", x$lfc, ")"))
  tbl <- x$de_table
  tbl$padj <- paste("<", tbl$padj)
  out <- utils::capture.output(print(tbl, row.names = FALSE))
  cli::cat_line(paste0("  ", out))

  cli::cat_line()
  cli::cat_line(cli::style_bold("  log2 fold change"))
  q <- x$lfc_quantiles
  cli::cat_line("  ", paste(c("min", "Q1", "median", "Q3", "max"),
                            sprintf("%+.2f", q), sep = " ", collapse = "  "))

  if (!is.null(x$size_factors)) {
    cli::cat_line()
    cli::cat_line(cli::style_bold("  Size factors"))
    s <- x$size_factors
    cli::cat_line("  ", paste(c("min", "median", "mean", "max"),
                              sprintf("%.2f", s[c("Min.", "Median", "Mean", "Max.")]),
                              sep = " ", collapse = "  "))
  }

  cli::cat_line()
  cli::cat_line(cli::col_grey(paste0(
    "  ", x$n_tested, " genes tested, ", x$n_na_padj,
    " with padj = NA (low count / outliers, removed by DESeq2)")))
  cli::cat_rule()
  invisible(x)
}


# ---- conversion -------------------------------------------------------------

#' Convert a run_deseq result to a table
#'
#' Returns the DESeq2 result table (`x$deseq`).
#'
#' @param x A `run_deseq_result` returned by [run_deseq()].
#' @param ... Passed to `as.data.frame()` or `tibble::as_tibble()`.
#' @return A data frame or tibble.
#' @name run_deseq_result-conversion
NULL

#' @rdname run_deseq_result-conversion
#' @export
as.data.frame.run_deseq_result <- function(x, ...) {
  as.data.frame(x$deseq, ...)
}

#' @rdname run_deseq_result-conversion
#' @exportS3Method tibble::as_tibble
as_tibble.run_deseq_result <- function(x, ...) {
  tibble::as_tibble(x$deseq, ...)
}
