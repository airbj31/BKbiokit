# Common control-group words (case-insensitive, matched as whole words only)
#' @noRd
.ref_words <- c(
  "control", "ctrl", "ctl", "con",
  "vehicle", "veh", "dmso", "pbs", "saline", "untreated", "mock", "placebo",
  "wt", "wildtype", "wild_type", "parental",
  "scramble", "scr", "sinc", "shnc", "nc", "sictrl", "shctrl", "sgnt", "nt",
  "ev", "vector", "sham",
  "normal", "healthy", "hc",
  "baseline", "pre", "0h", "day0", "d0", "t0",
  "uninfected", "naive", "unstimulated"
)

# Match only the whole name, or a word separated by _ - . or a space
#' @noRd
.ref_pattern <- paste0("(^|[-_. ])(", paste(.ref_words, collapse = "|"), ")($|[-_. ])")


# Decide the reference group -> list(grp = factor, how = description)
#' @noRd
.pick_reference <- function(grp, condition) {
  if (length(unique(grp)) != 2) {
    stop(condition, " must have exactly 2 groups, but has: ",
         paste(unique(grp), collapse = ", "), call. = FALSE)
  }

  if (is.factor(grp)) {
    grp <- droplevels(grp)
    lv  <- levels(grp)
    looks_ref <- lv[grepl(.ref_pattern, lv, ignore.case = TRUE)]
    if (length(looks_ref) == 1 && looks_ref != lv[1]) {
      warning(condition, " is a factor with reference '", lv[1], "', but '",
              looks_ref, "' looks like the control group.\n",
              "  Factor levels are alphabetical unless set on purpose. If '",
              looks_ref, "' should be the reference:\n",
              "  meta$", condition, " <- relevel(meta$", condition,
              ", ref = \"", looks_ref, "\")",
              call. = FALSE)
    }
    return(list(grp = grp,
                how = "factor level order (first level is the reference)"))
  }

  lv  <- unique(as.character(grp))
  ref <- lv[grepl(.ref_pattern, lv, ignore.case = TRUE)]

  if (length(ref) != 1) {
    stop("Could not detect the control group automatically (candidates: ",
         if (length(ref)) paste(ref, collapse = ", ") else "none", ").\n",
         "  Make ", condition, " a factor to set the reference, e.g.:\n",
         "  meta$", condition, " <- factor(meta$", condition,
         ", levels = c(\"", lv[1], "\", \"", lv[2], "\"))  # first level is the reference",
         call. = FALSE)
  }

  hit <- regmatches(ref, regexpr(.ref_pattern, ref, ignore.case = TRUE))
  hit <- gsub("[-_. ]", "", hit)

  list(grp = factor(grp, levels = c(ref, setdiff(lv, ref))),
       how = paste0("auto-detected (matched '", hit, "')"))
}


#' Run a two-group DESeq2 comparison
#'
#' Fits a DESeq2 model for a condition with exactly two groups, picks the
#' reference (control) group automatically, and returns the results as a
#' tibble together with normalized counts.
#'
#' The reference group is chosen as follows:
#' * If `meta[[condition]]` is a factor, its first level is the reference.
#' * Otherwise, the group whose name contains a common control word
#'   (e.g. control, ctrl, vehicle, DMSO, PBS, WT, siNC, scramble, sham, mock)
#'   as a whole word is the reference.
#' * If no group, or both groups, match, the function stops and asks you to
#'   supply `condition` as a factor.
#'
#' @param cnt Raw count matrix (or data frame) with genes in rows and samples
#'   in columns. Row names are gene IDs and column names are sample IDs.
#' @param meta Data frame of sample information with sample IDs as row names
#'   (matching `colnames(cnt)`, in any order). A tibble has no row names, so
#'   convert it first, e.g. `tibble::column_to_rownames(meta, "sample_id")`.
#' @param condition Name of the column in `meta` to compare. Must have exactly
#'   two groups.
#' @param covariates Optional character vector of columns in `meta` to adjust
#'   for (e.g. `"batch"`). They are added before `condition` in the design.
#' @param id Name of the gene ID column in the result table.
#' @param cnt_thres Genes with a total count below this value are removed
#'   before fitting.
#' @param verbose If `TRUE`, print the comparison setup and DESeq2 progress.
#'
#' @return A list with:
#' * `dds`: the fitted `DESeqDataSet`.
#' * `norm_cnt`: matrix of normalized counts.
#' * `size_factors`: named vector of DESeq2 size factors (one per sample).
#'   `norm_cnt` equals raw counts divided by these.
#' * `deseq`: tibble of DESeq2 results (case vs control).
#' * `comparison`: list describing the comparison (`condition`, `case`,
#'   `control`, `method`, `n`, `design`, `genes`, `cnt_thres`, `id`).
#'
#' The list has class `run_deseq_result`. Printing it shows a short summary
#' with the top DE genes; use `$` to get the full tables, [summary()] for DE
#' counts at several cutoffs, and `as.data.frame()` / `tibble::as_tibble()`
#' for the result table.
#'
#' @examples
#' \dontrun{
#' meta <- tibble::column_to_rownames(meta_tbl, "sample_id")
#'
#' res <- cnt_mat |> run_deseq(meta, condition = "group")
#'
#' # Set the reference yourself with a factor
#' meta$group <- factor(meta$group, levels = c("WT", "KO"))
#' res <- run_deseq(cnt_mat, meta, "group", covariates = "batch")
#' }
#'
#' @export
run_deseq <- function(cnt, meta, condition,
                      covariates = NULL, id = "gene_id",
                      cnt_thres = 10, verbose = TRUE) {

  ## ---- Basic checks ----
  if (!is.character(condition) || length(condition) != 1 || is.na(condition) || condition == "") {
    stop("condition must be a single column name, e.g. \"group\".", call. = FALSE)
  }
  if (!is.null(covariates) &&
      (!is.character(covariates) || anyNA(covariates) || any(covariates == ""))) {
    stop("covariates must be NULL or column names, e.g. c(\"batch\"). ",
         "Use NULL (not \"\") when there are no covariates.", call. = FALSE)
  }

  cols_hint <- paste0("  Columns in meta: ", paste(names(meta), collapse = ", "), "\n",
                      "  (Sample names belong in rownames(meta), not in these arguments.)")
  if (!condition %in% names(meta)) {
    stop("condition = \"", condition, "\" is not a column of meta.\n", cols_hint,
         call. = FALSE)
  }
  bad_cov <- setdiff(covariates, names(meta))
  if (length(bad_cov) > 0) {
    stop("covariates not found as columns of meta: ",
         paste0("\"", bad_cov, "\"", collapse = ", "), "\n", cols_hint, call. = FALSE)
  }

  ## ---- Align meta to cnt sample order (by row names) ----
  if (is.null(colnames(cnt))) {
    stop("cnt must have sample IDs as column names.", call. = FALSE)
  }

  not_in_meta <- setdiff(colnames(cnt), rownames(meta))
  if (length(not_in_meta) > 0) {
    stop("Sample(s) in cnt not found in rownames(meta): ",
         paste(not_in_meta, collapse = ", "), "\n",
         "  meta must have sample IDs as row names. If they are in a column, use e.g.:\n",
         "  meta <- tibble::column_to_rownames(meta, \"sample_id\")",
         call. = FALSE)
  }
  meta <- as.data.frame(meta)

  # Report anything that differs before aligning
  only_meta <- setdiff(rownames(meta), colnames(cnt))
  if (verbose && length(only_meta) > 0) {
    message("Note: ", length(only_meta), " sample(s) in meta but not in cnt were dropped: ",
            paste(only_meta, collapse = ", "))
  }
  kept <- rownames(meta)[rownames(meta) %in% colnames(cnt)]
  if (verbose && !identical(kept, colnames(cnt))) {
    message("Note: meta rows were reordered to match colnames(cnt).")
  }

  meta <- meta[colnames(cnt), , drop = FALSE]
  stopifnot(identical(rownames(meta), colnames(cnt)))

  ## ---- Set the comparison groups ----
  picked <- .pick_reference(meta[[condition]], condition)
  grp    <- picked$grp
  meta[[condition]] <- grp

  design <- stats::reformulate(c(covariates, condition))
  n      <- table(grp)

  comparison <- list(
    condition = condition,
    case      = levels(grp)[2],
    control   = levels(grp)[1],
    method    = picked$how,
    n         = c(n),
    design    = deparse(design)
  )

  if (verbose) {
    message(
      "-- DESeq2 comparison --\n",
      "  condition : ", condition, "\n",
      "  case      : ", comparison$case,    " (n = ", n[[2]], ")\n",
      "  control   : ", comparison$control, " (n = ", n[[1]], ")\n",
      "  reference : ", comparison$method, "\n",
      "  design    : ", comparison$design, "\n",
      "  log2FC > 0 means higher expression in ", comparison$case
    )
  }

  ## ---- DESeq2 ----
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(as.matrix(cnt)),
    colData   = meta,
    design    = design
  )
  n_input <- nrow(dds)
  dds <- dds[rowSums(DESeq2::counts(dds)) >= cnt_thres, ]
  comparison$genes <- c(input = n_input, tested = nrow(dds))
  comparison$cnt_thres <- cnt_thres
  comparison$id <- id
  dds <- DESeq2::DESeq(dds, quiet = !verbose)

  res <- DESeq2::results(dds, contrast = c(condition, comparison$case, comparison$control)) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = id) |>
    tibble::as_tibble()

  out <- list(
    dds        = dds,
    norm_cnt     = DESeq2::counts(dds, normalized = TRUE),
    size_factors = DESeq2::sizeFactors(dds),
    deseq        = res,
    comparison = comparison
  )
  class(out) <- c("run_deseq_result", "list")
  out
}
