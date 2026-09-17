# Reusable, local over-representation analysis for GMT files or named gene sets.
# Source this file; it does not execute an analysis automatically.

.ora_stop <- function(...) stop(..., call. = FALSE)

.ora_stable_unique <- function(x) x[!duplicated(x)]

.ora_normalize_ids <- function(x, label, allow_whole_number_numeric = FALSE) {
  if (is.factor(x) || !is.atomic(x) || is.list(x)) {
    .ora_stop(label, " must be an atomic character or integer vector.")
  }
  if (is.double(x)) {
    if (!allow_whole_number_numeric) {
      .ora_stop(
        label,
        " is numeric. Convert IDs explicitly with as.character(), or set ",
        "allow_whole_number_numeric = TRUE after checking for precision loss."
      )
    }
    if (any(!is.finite(x)) || any(x != floor(x)) || any(abs(x) > 2^53)) {
      .ora_stop(label, " contains unsafe numeric identifiers.")
    }
    x <- format(x, scientific = FALSE, trim = TRUE)
  } else if (is.integer(x)) {
    if (anyNA(x)) .ora_stop(label, " contains NA identifiers.")
    x <- as.character(x)
  } else if (!is.character(x)) {
    .ora_stop(label, " must be character or integer.")
  }
  if (anyNA(x)) .ora_stop(label, " contains NA identifiers.")
  if (any(!nzchar(x))) .ora_stop(label, " contains blank identifiers.")
  if (any(x != trimws(x))) {
    .ora_stop(label, " contains leading/trailing whitespace.")
  }
  if (any(grepl("[[:cntrl:];]", x))) {
    .ora_stop(label, " contains an unsafe semicolon or control delimiter.")
  }
  x
}

.ora_validate_pathway_names <- function(x) {
  if (is.null(names(x)) || length(names(x)) != length(x)) {
    .ora_stop("gene_sets must be a named list.")
  }
  if (anyNA(names(x)) || any(!nzchar(names(x)))) {
    .ora_stop("Pathway names must be nonblank.")
  }
  if (any(names(x) != trimws(names(x)))) {
    .ora_stop("Pathway names must not have leading/trailing whitespace.")
  }
  if (any(grepl("[[:cntrl:]]", names(x)))) {
    .ora_stop("Pathway names contain control delimiters.")
  }
  if (anyDuplicated(names(x))) .ora_stop("Pathway names must be unique.")
  invisible(TRUE)
}

.ora_assert_bk_parity <- function(raw_sets, bk_sets, pathway_names) {
  if (!is.list(bk_sets) || !identical(names(bk_sets), pathway_names)) {
    .ora_stop("BKbiokit::read_GMT() names/order disagree with the validated GMT rows.")
  }
  same_members <- vapply(
    seq_along(raw_sets),
    function(i) identical(as.character(bk_sets[[i]]), raw_sets[[i]]),
    logical(1L)
  )
  if (!all(same_members)) {
    .ora_stop("BKbiokit::read_GMT() membership disagrees with the raw GMT parser.")
  }
  invisible(TRUE)
}

.ora_read_gmt_checked <- function(path, allow_whole_number_numeric = FALSE) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) {
    .ora_stop("GMT file does not exist: ", paste(path, collapse = ", "))
  }
  if (!requireNamespace("BKbiokit", quietly = TRUE)) {
    .ora_stop("BKbiokit is required to read GMT files.")
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (!length(lines)) .ora_stop("GMT file is empty: ", path)
  fields <- lapply(lines, function(line) {
    z <- strsplit(paste0(line, "\t__ORA_END__"), "\t", fixed = TRUE)[[1L]]
    z[-length(z)]
  })
  widths <- lengths(fields)
  if (any(widths < 3L)) {
    .ora_stop("Every GMT row must contain pathway, description, and at least one gene.")
  }
  pathway_names <- vapply(fields, `[[`, character(1L), 1L)
  descriptions <- vapply(fields, `[[`, character(1L), 2L)
  if (any(descriptions != trimws(descriptions))) {
    .ora_stop("GMT descriptions must not have leading/trailing whitespace.")
  }
  raw_sets <- lapply(fields, function(z) z[-c(1L, 2L)])
  names(raw_sets) <- pathway_names
  .ora_validate_pathway_names(raw_sets)
  for (i in seq_along(raw_sets)) {
    raw_sets[[i]] <- .ora_normalize_ids(
      raw_sets[[i]],
      paste0("GMT pathway '", pathway_names[[i]], "'"),
      allow_whole_number_numeric
    )
    if (any(raw_sets[[i]] != trimws(raw_sets[[i]]))) {
      .ora_stop("GMT pathway '", pathway_names[[i]], "' contains whitespace-padded IDs.")
    }
  }

  bk_sets <- BKbiokit::read_GMT(path, simple = TRUE)
  .ora_assert_bk_parity(raw_sets, bk_sets, pathway_names)
  list(
    pathways = raw_sets,
    descriptions = descriptions,
    source_type = "gmt",
    source = normalizePath(path, mustWork = TRUE)
  )
}

.ora_prepare_gene_sets <- function(gene_sets, allow_whole_number_numeric = FALSE) {
  if (is.character(gene_sets) && length(gene_sets) == 1L && file.exists(gene_sets)) {
    return(.ora_read_gmt_checked(gene_sets, allow_whole_number_numeric))
  }
  if (!is.list(gene_sets)) {
    .ora_stop("gene_sets must be a GMT filepath or a named list of gene vectors.")
  }
  .ora_validate_pathway_names(gene_sets)
  normalized <- lapply(seq_along(gene_sets), function(i) {
    .ora_normalize_ids(
      gene_sets[[i]],
      paste0("pathway '", names(gene_sets)[[i]], "'"),
      allow_whole_number_numeric
    )
  })
  names(normalized) <- names(gene_sets)
  list(
    pathways = normalized,
    descriptions = rep(NA_character_, length(normalized)),
    source_type = "named_list",
    source = NA_character_
  )
}

.ora_fgsea_metadata <- function() {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    .ora_stop(
      "engine='fgsea' was requested, but fgsea is not installed in: ",
      paste(.libPaths(), collapse = "; ")
    )
  }
  desc <- utils::packageDescription("fgsea")
  built <- unname(desc[["Built"]] %||% "")
  built_r <- sub("^R ([0-9]+\\.[0-9]+).*$", "\\1", built)
  active_r <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][[1L]], sep = ".")
  if (!nzchar(built) || identical(built_r, built) || built_r != active_r) {
    .ora_stop(
      "Installed fgsea was not built for the active R major/minor version. Built='",
      built, "', active R='", active_r, "'."
    )
  }
  if (!("fora" %in% getNamespaceExports("fgsea"))) {
    .ora_stop("Installed fgsea does not export fora().")
  }
  list(
    version = as.character(utils::packageVersion("fgsea")),
    built = built,
    library_path = dirname(find.package("fgsea"))
  )
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x)) y else x

.ora_sha256_file <- function(path) {
  command <- Sys.which("shasum")
  args <- c("-a", "256", shQuote(normalizePath(path, mustWork = TRUE)))
  if (!nzchar(command)) {
    command <- Sys.which("sha256sum")
    args <- shQuote(normalizePath(path, mustWork = TRUE))
  }
  if (!nzchar(command)) .ora_stop("Neither shasum nor sha256sum is available.")
  output <- system2(command, args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status") %||% 0L
  if (status != 0L || !length(output)) .ora_stop("SHA-256 calculation failed.")
  hash <- tolower(strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]])
  if (!grepl("^[0-9a-f]{64}$", hash)) .ora_stop("Invalid SHA-256 output.")
  hash
}

.ora_sha256_ids <- function(ids) {
  path <- tempfile("ora_ids_", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(ids, path, sep = "\n", useBytes = TRUE)
  .ora_sha256_file(path)
}

.ora_test_status <- function(K, M, min_size, max_size) {
  if (K == 0L) return("empty_after_background")
  if (K < min_size) return("below_min_size")
  if (K == M) return("equals_background")
  if (K > max_size) return("above_max_size")
  "tested"
}

#' Run local over-representation analysis.
#'
#' @param gene_sets GMT filepath or named list of gene-ID vectors.
#' @param query_genes selected/foreground gene IDs.
#' @param background_genes all genes eligible to enter query_genes for the same
#'   assay and contrast after the same filtering and identifier mapping.
#' @param min_size,max_size pathway size limits after intersection with background.
#'   max_size is capped at the number of unique background genes minus one;
#'   NULL uses that cap. min_size must not exceed the effective max_size.
#' @param engine deterministic engine: "base" (default) or explicitly "fgsea".
#' @param allow_whole_number_numeric permit explicitly audited whole-number doubles.
#' @return list with results (tested pathways only), coverage, pathway_audit
#'   (all pathways, descriptions and filtering status), gene_audit, and parameters
#'   (shared analysis settings and engine metadata). Join results to pathway_audit
#'   by pathway or pathway_index; pathway_index preserves the original input index.
#' @export
run_ora <- function(
    gene_sets,
    query_genes,
    background_genes,
    min_size = 1L,
    max_size = NULL,
    engine = c("base", "fgsea"),
    allow_whole_number_numeric = FALSE
) {
  engine <- match.arg(engine)
  prepared <- .ora_prepare_gene_sets(gene_sets, allow_whole_number_numeric)

  query_raw <- .ora_normalize_ids(
    query_genes, "query_genes", allow_whole_number_numeric
  )
  background_raw <- .ora_normalize_ids(
    background_genes, "background_genes", allow_whole_number_numeric
  )
  query <- .ora_stable_unique(query_raw)
  background <- .ora_stable_unique(background_raw)
  if (!length(query)) .ora_stop("query_genes is empty after validation.")
  if (length(background) < 2L) .ora_stop("background_genes must contain at least two unique IDs.")
  outside_query <- setdiff(query, background)
  if (length(outside_query)) {
    .ora_stop(
      length(outside_query), " unique query ID(s) are outside background, for example: ",
      paste(utils::head(outside_query, 10L), collapse = ", ")
    )
  }

  M <- length(background)
  n_query <- length(query)
  if (length(min_size) != 1L || is.na(min_size) || min_size != as.integer(min_size)) {
    .ora_stop("min_size must be one integer.")
  }
  min_size <- as.integer(min_size)
  if (is.null(max_size)) max_size <- M - 1L
  if (length(max_size) != 1L || is.na(max_size) || max_size != as.integer(max_size)) {
    .ora_stop("max_size must be one integer.")
  }
  max_size <- min(as.integer(max_size), M - 1L)
  if (min_size < 1L || max_size < min_size) {
    .ora_stop("Require 1 <= min_size <= effective max_size (capped at background size - 1).")
  }

  pathways_raw <- prepared$pathways
  pathways_unique <- lapply(pathways_raw, .ora_stable_unique)
  pathways_bg <- lapply(pathways_unique, function(z) z[z %in% background])
  outside_bg <- Map(setdiff, pathways_unique, pathways_bg)
  raw_sizes <- lengths(pathways_raw)
  unique_sizes <- lengths(pathways_unique)
  K <- lengths(pathways_bg)
  overlap_genes <- lapply(pathways_bg, function(z) query[query %in% z])
  k <- lengths(overlap_genes)
  statuses <- vapply(K, .ora_test_status, character(1L), M, min_size, max_size)
  tested <- statuses == "tested"
  if (!any(tested)) .ora_stop("No pathways remain testable after background and size filtering.")

  pvalue <- rep(NA_real_, length(pathways_raw))
  fold <- rep(NA_real_, length(pathways_raw))
  expected <- rep(NA_real_, length(pathways_raw))
  pvalue[tested] <- stats::phyper(
    q = k[tested] - 1L,
    m = K[tested],
    n = M - K[tested],
    k = n_query,
    lower.tail = FALSE
  )
  expected[tested] <- n_query * K[tested] / M
  fold[tested] <- (k[tested] / n_query) / (K[tested] / M)
  padj <- rep(NA_real_, length(pathways_raw))
  padj[tested] <- stats::p.adjust(pvalue[tested], method = "BH")

  engine_version <- paste0("R-", as.character(getRversion()), "/stats")
  engine_built <- R.version.string
  engine_library_path <- normalizePath(.Library, mustWork = TRUE)
  if (engine == "fgsea") {
    meta <- .ora_fgsea_metadata()
    fora_result <- fgsea::fora(
      pathways = pathways_bg[tested],
      genes = query,
      universe = background,
      minSize = min_size,
      maxSize = max_size
    )
    fora_result <- as.data.frame(fora_result, stringsAsFactors = FALSE)
    required <- c("pathway", "pval", "padj", "foldEnrichment", "overlap", "size", "overlapGenes")
    if (!all(required %in% names(fora_result))) {
      .ora_stop("fgsea::fora() output schema is incompatible with this script.")
    }
    tested_names <- names(pathways_bg)[tested]
    idx <- match(tested_names, fora_result$pathway)
    if (anyNA(idx) || nrow(fora_result) != length(tested_names) || anyDuplicated(fora_result$pathway)) {
      .ora_stop("fgsea::fora() tested-pathway set disagrees with the base preprocessing.")
    }
    fr <- fora_result[idx, , drop = FALSE]
    close_num <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = FALSE))
    same_edges <- vapply(
      seq_along(tested_names),
      function(i) setequal(as.character(fr$overlapGenes[[i]]), overlap_genes[tested][[i]]),
      logical(1L)
    )
    if (
      !identical(as.integer(fr$size), as.integer(K[tested])) ||
      !identical(as.integer(fr$overlap), as.integer(k[tested])) ||
      !close_num(fr$pval, pvalue[tested]) ||
      !close_num(fr$padj, padj[tested]) ||
      !close_num(fr$foldEnrichment, fold[tested]) ||
      !all(same_edges)
    ) {
      .ora_stop("fgsea::fora() failed parity with the pinned base hypergeometric contract.")
    }
    pvalue[tested] <- fr$pval
    padj[tested] <- fr$padj
    fold[tested] <- fr$foldEnrichment
    engine_version <- paste0("fgsea-", meta$version)
    engine_built <- meta$built
    engine_library_path <- meta$library_path
  }

  background_nonquery_in_pathway <- K - k
  query_not_in_pathway <- n_query - k
  background_nonquery_not_in_pathway <- M - n_query - K + k
  if (any(background_nonquery_not_in_pathway < 0L)) {
    .ora_stop("Internal contingency-table invariant failed.")
  }

  results <- data.frame(
    pathway = names(pathways_raw),
    pathway_index = seq_along(pathways_raw),
    pathway_size_in_background_K = K,
    background_size_M = M,
    query_size_n = n_query,
    overlap_size_k = k,
    query_not_in_pathway = query_not_in_pathway,
    background_nonquery_in_pathway = background_nonquery_in_pathway,
    background_nonquery_not_in_pathway = background_nonquery_not_in_pathway,
    expected_overlap = expected,
    fold_enrichment = fold,
    pvalue = pvalue,
    padj = padj,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  results$overlap_genes <- I(overlap_genes)
  results <- results[tested, , drop = FALSE]
  rownames(results) <- NULL

  union_raw <- unlist(pathways_raw, use.names = FALSE)
  union_unique <- .ora_stable_unique(union_raw)
  union_tested <- .ora_stable_unique(unlist(pathways_bg[tested], use.names = FALSE))
  status_table <- table(factor(
    statuses,
    levels = c("tested", "empty_after_background", "below_min_size", "equals_background", "above_max_size")
  ))
  coverage <- data.frame(
    metric = c(
      "query_raw", "query_unique", "query_duplicates_removed",
      "background_raw", "background_unique", "background_duplicates_removed",
      "query_outside_background", "pathway_count", "pathway_tested",
      "pathway_empty_after_background", "pathway_below_min_size",
      "pathway_equals_background", "pathway_above_max_size",
      "pathway_membership_raw", "pathway_membership_duplicates_removed",
      "pathway_membership_outside_background", "pathway_union_unique",
      "pathway_union_in_background", "background_covered_by_any_pathway",
      "background_covered_by_tested_pathway", "query_covered_by_any_pathway",
      "query_covered_by_tested_pathway"
    ),
    value = c(
      length(query_raw), length(query), length(query_raw) - length(query),
      length(background_raw), M, length(background_raw) - M,
      0L, length(pathways_raw), unname(status_table[["tested"]]),
      unname(status_table[["empty_after_background"]]),
      unname(status_table[["below_min_size"]]),
      unname(status_table[["equals_background"]]),
      unname(status_table[["above_max_size"]]),
      length(union_raw), sum(raw_sizes - unique_sizes), sum(lengths(outside_bg)),
      length(union_unique), sum(union_unique %in% background),
      sum(background %in% union_unique), sum(background %in% union_tested),
      sum(query %in% union_unique), sum(query %in% union_tested)
    ),
    stringsAsFactors = FALSE
  )
  coverage$denominator <- c(
    length(query_raw), length(query), length(query_raw),
    length(background_raw), M, length(background_raw),
    length(query), rep(length(pathways_raw), 6L),
    length(union_raw), length(union_raw), sum(unique_sizes),
    length(union_unique), length(union_unique), M, M, length(query), length(query)
  )
  coverage$percent <- ifelse(
    coverage$denominator > 0,
    100 * coverage$value / coverage$denominator,
    NA_real_
  )

  pathway_audit <- data.frame(
    pathway = names(pathways_raw),
    pathway_index = seq_along(pathways_raw),
    description = prepared$descriptions,
    test_status = statuses,
    raw_pathway_size = raw_sizes,
    unique_pathway_size = unique_sizes,
    duplicate_memberships_removed = raw_sizes - unique_sizes,
    pathway_size_in_background = K,
    outside_background_count = lengths(outside_bg),
    stringsAsFactors = FALSE
  )
  pathway_audit$outside_background_genes <- I(outside_bg)

  gene_audit <- data.frame(
    gene_id = background,
    background_order = seq_along(background),
    raw_query_occurrence_count = vapply(
      background, function(id) sum(query_raw == id), integer(1L)
    ),
    raw_background_occurrence_count = vapply(
      background, function(id) sum(background_raw == id), integer(1L)
    ),
    in_query = background %in% query,
    in_background = TRUE,
    covered_by_any_pathway = background %in% union_unique,
    covered_by_tested_pathway = background %in% union_tested,
    stringsAsFactors = FALSE
  )
  if (
    sum(gene_audit$raw_query_occurrence_count) != length(query_raw) ||
    sum(gene_audit$raw_background_occurrence_count) != length(background_raw) ||
    nrow(gene_audit) != M ||
    sum(gene_audit$in_query) != n_query ||
    !all(query %in% gene_audit$gene_id)
  ) {
    .ora_stop("Gene audit reconciliation failed.")
  }

  query_unique_sha256 <- .ora_sha256_ids(query)
  background_unique_sha256 <- .ora_sha256_ids(background)

  parameters <- data.frame(
    key = c(
      "source_type", "source", "min_size", "max_size", "p_adjust_method",
      "requested_engine", "resolved_engine", "engine_version", "engine_built",
      "engine_library_path", "R_version", "query_unique_sha256",
      "background_unique_sha256", "id_hash_canonicalization",
      "background_definition"
    ),
    value = c(
      prepared$source_type, prepared$source, min_size, max_size, "BH",
      engine, if (engine == "base") "base::phyper" else "fgsea::fora",
      engine_version, engine_built, engine_library_path, R.version.string,
      query_unique_sha256, background_unique_sha256,
      "stable first-occurrence unique ID sequence; one UTF-8 ID per line with final LF",
      "Caller-supplied genes eligible to enter query after the same assay/contrast filtering and ID mapping"
    ),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      results = results,
      coverage = coverage,
      pathway_audit = pathway_audit,
      gene_audit = gene_audit,
      parameters = parameters
    ),
    class = c("ora_result", "list")
  )
}

.ora_serialize_list_column <- function(x) {
  vapply(x, function(z) paste(z, collapse = ";"), character(1L))
}

.ora_character_table <- function(x) {
  out <- as.data.frame(
    lapply(x, function(z) {
      y <- as.character(z)
      y[is.na(y)] <- ""
      y
    }),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(out) <- names(x)
  out
}

.ora_read_tsv_character <- function(path) {
  utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    colClasses = "character", na.strings = character(), check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Write all ORA tables with a common prefix, without overwriting by default.
write_ora_results <- function(x, output_prefix) {
  if (!inherits(x, "ora_result")) .ora_stop("x must be returned by run_ora().")
  if (length(output_prefix) != 1L || is.na(output_prefix) || !nzchar(output_prefix)) {
    .ora_stop("output_prefix must be one nonblank path.")
  }
  paths <- paste0(
    output_prefix,
    c("_ora.tsv", "_coverage.tsv", "_pathway_audit.tsv", "_gene_audit.tsv", "_parameters.tsv")
  )
  if (any(file.exists(paths))) {
    .ora_stop("Refusing to overwrite existing output: ", paste(paths[file.exists(paths)], collapse = ", "))
  }
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  tables <- list(x$results, x$coverage, x$pathway_audit, x$gene_audit, x$parameters)
  tables[[1L]]$overlap_genes <- .ora_serialize_list_column(tables[[1L]]$overlap_genes)
  tables[[3L]]$outside_background_genes <- .ora_serialize_list_column(
    tables[[3L]]$outside_background_genes
  )
  staged <- vapply(paths, function(p) tempfile(paste0(basename(p), ".stage."), dirname(p)), character(1L))
  on.exit(unlink(staged[file.exists(staged)]), add = TRUE)
  for (i in seq_along(paths)) {
    utils::write.table(
      tables[[i]], staged[[i]], sep = "\t", quote = FALSE, row.names = FALSE,
      na = "", fileEncoding = "UTF-8"
    )
    if (!identical(.ora_read_tsv_character(staged[[i]]), .ora_character_table(tables[[i]]))) {
      .ora_stop("Staged TSV character read-back failed: ", basename(paths[[i]]))
    }
  }
  staged_hashes <- unname(tools::md5sum(staged))
  committed <- logical(length(paths))
  success <- FALSE
  on.exit(if (!success) unlink(paths[committed & file.exists(paths)]), add = TRUE)
  fail_after <- getOption("ora.test_fail_after_commit", Inf)
  for (i in seq_along(paths)) {
    if (!file.rename(staged[[i]], paths[[i]])) {
      unlink(paths[committed])
      .ora_stop("Failed to commit ORA output transaction.")
    }
    committed[[i]] <- TRUE
    if (i >= fail_after) .ora_stop("Injected ORA commit failure for self-test.")
  }
  if (!identical(unname(tools::md5sum(paths)), staged_hashes)) {
    .ora_stop("Final ORA output hashes disagree with staged files.")
  }
  for (i in seq_along(paths)) {
    if (!identical(.ora_read_tsv_character(paths[[i]]), .ora_character_table(tables[[i]]))) {
      .ora_stop("Final TSV character read-back failed: ", basename(paths[[i]]))
    }
  }
  success <- TRUE
  invisible(paths)
}

run_ora_self_tests <- function() {
  expect_error <- function(expr) {
    observed <- FALSE
    tryCatch(force(expr), error = function(e) observed <<- TRUE)
    if (!observed) .ora_stop("Expected an error but none occurred.")
  }
  sets <- list(
    A = c("a", "b", "c", "c", "outside"),
    B = c("d"),
    C = c("a", "b", "c", "d", "e", "f"),
    D = c("outside"),
    E = c("b", "c"),
    F = c("a", "b", "c", "d", "e")
  )
  background <- c("a", "b", "c", "d", "e", "f", "f")
  query <- c("a", "d", "a")
  x <- run_ora(sets, query, background, min_size = 2L, max_size = 5L)
  stopifnot(
    inherits(x, "ora_result"), nrow(x$results) == 3L,
    identical(x$results$pathway, c("A", "E", "F")),
    identical(x$results$pathway_index, c(1L, 5L, 6L)),
    identical(x$pathway_audit$pathway, names(sets)),
    identical(x$pathway_audit$test_status, c("tested", "below_min_size", "equals_background", "empty_after_background", "tested", "tested")),
    all(is.na(x$pathway_audit$description)),
    x$results$pvalue[x$results$pathway == "A"] == stats::phyper(0L, 3L, 3L, 2L, lower.tail = FALSE),
    x$results$pvalue[x$results$pathway == "E"] == 1,
    x$results$fold_enrichment[x$results$pathway == "E"] == 0,
    x$pathway_audit$raw_pathway_size[x$pathway_audit$pathway == "A"] == 5L,
    x$pathway_audit$unique_pathway_size[x$pathway_audit$pathway == "A"] == 4L,
    x$pathway_audit$outside_background_count[x$pathway_audit$pathway == "A"] == 1L,
    identical(x$results$padj, stats::p.adjust(x$results$pvalue, "BH")),
    !any(c("description", "test_status", "raw_pathway_size", "unique_pathway_size",
           "outside_background_count", "min_size", "max_size", "requested_engine",
           "resolved_engine", "engine_version", "engine_built", "engine_library_path") %in%
           names(x$results)),
    all(c("min_size", "max_size", "requested_engine", "resolved_engine",
          "engine_version", "engine_built", "engine_library_path") %in% x$parameters$key),
    nrow(x$results) == x$coverage$value[x$coverage$metric == "pathway_tested"]
  )

  above <- run_ora(sets, query, background, min_size = 1L, max_size = 4L)
  stopifnot(
    above$pathway_audit$test_status[above$pathway_audit$pathway == "F"] == "above_max_size",
    !("F" %in% above$results$pathway)
  )

  capped <- run_ora(sets, query, background, min_size = 2L, max_size = 500L)
  stopifnot(
    identical(capped$results, x$results),
    capped$parameters$value[capped$parameters$key == "max_size"] == "5"
  )
  expect_error(run_ora(sets, query, background, min_size = 6L, max_size = 500L))

  dedup_free <- run_ora(
    lapply(sets, .ora_stable_unique), unique(query), unique(background),
    min_size = 2L, max_size = 5L
  )
  stopifnot(
    identical(x$results$pvalue, dedup_free$results$pvalue),
    identical(x$results$padj, dedup_free$results$padj)
  )

  tested_rows <- x$results
  for (i in seq_len(nrow(tested_rows))) {
    tab <- matrix(c(
      tested_rows$overlap_size_k[[i]],
      tested_rows$query_not_in_pathway[[i]],
      tested_rows$background_nonquery_in_pathway[[i]],
      tested_rows$background_nonquery_not_in_pathway[[i]]
    ), nrow = 2L)
    fisher_p <- stats::fisher.test(tab, alternative = "greater")$p.value
    stopifnot(isTRUE(all.equal(fisher_p, tested_rows$pvalue[[i]], tolerance = 1e-14)))
  }

  for (M0 in 2:10) for (n0 in seq_len(M0)) for (K0 in seq_len(M0 - 1L)) {
    lower <- max(0L, n0 - (M0 - K0))
    upper <- min(n0, K0)
    for (k0 in lower:upper) {
      p_hyper <- stats::phyper(k0 - 1L, K0, M0 - K0, n0, lower.tail = FALSE)
      tab <- matrix(c(k0, K0 - k0, n0 - k0, M0 - K0 - n0 + k0), nrow = 2L)
      p_fisher <- stats::fisher.test(tab, alternative = "greater")$p.value
      stopifnot(isTRUE(all.equal(p_hyper, p_fisher, tolerance = 1e-13)))
    }
  }

  expect_error(run_ora(sets, c("a", "not_in_background"), background))
  expect_error(run_ora(list(A = c("a"), A = c("b")), "a", background))
  expect_error(run_ora(list(A = c("a", "")), "a", background))
  expect_error(run_ora(list(A = c("a", NA_character_)), "a", background))
  expect_error(run_ora(list(A = c("a\tb")), "a", background))
  expect_error(run_ora(list(A = c("a ", "b")), "a", background))
  expect_error(run_ora(sets, "a ", background))
  expect_error(run_ora(sets, "a", c(background, "z ")))
  expect_error(run_ora(sets, "a", c(background, "   ")))
  expect_error(run_ora(list(A = c(1, 2)), c(1, 2), c(1, 2, 3)))
  expect_error(run_ora(list(A = c("a")), "a", background, min_size = 4L, max_size = 5L))

  td <- tempfile("ora_self_test_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  good_gmt <- file.path(td, "good.gmt")
  writeLines(c("A\tdesc A\ta\tb\tb", "B\tdesc B\tc\td", "C\tdesc excluded\toutside"), good_gmt, useBytes = TRUE)
  gmt_run <- run_ora(good_gmt, "a", c("a", "b", "c", "d"))
  list_run <- run_ora(list(A = c("a", "b", "b"), B = c("c", "d")), "a", c("a", "b", "c", "d"))
  stopifnot(
    identical(gmt_run$results$pvalue, list_run$results$pvalue),
    gmt_run$pathway_audit$duplicate_memberships_removed[[1L]] == 1L,
    identical(gmt_run$pathway_audit$description, c("desc A", "desc B", "desc excluded")),
    gmt_run$pathway_audit$test_status[[3L]] == "empty_after_background",
    identical(gmt_run$results$pathway, c("A", "B"))
  )
  write_bad <- function(name, line) {
    path <- file.path(td, name)
    writeLines(line, path, useBytes = TRUE)
    path
  }
  expect_error(run_ora(write_bad("trailing.gmt", "A\tdesc\ta\t"), "a", c("a", "b")))
  expect_error(run_ora(write_bad("short.gmt", "A\tdesc"), "a", c("a", "b")))
  expect_error(run_ora(write_bad("dup.gmt", c("A\tdesc\ta", "A\tdesc\tb")), "a", c("a", "b")))
  expect_error(run_ora(write_bad("blank_name.gmt", "\tdesc\ta"), "a", c("a", "b")))
  expect_error(run_ora(write_bad("space_id.gmt", "A\tdesc\t a"), "a", c("a", "b")))
  expect_error(.ora_assert_bk_parity(list(A = "a"), list(A = "b"), "A"))

  output_prefix <- file.path(td, "write_test")
  output_paths <- write_ora_results(list_run, output_prefix)
  gmt_paths <- write_ora_results(gmt_run, file.path(td, "write_gmt_test"))
  saved_results <- .ora_read_tsv_character(gmt_paths[[1L]])
  saved_audit <- .ora_read_tsv_character(gmt_paths[[3L]])
  saved_parameters <- .ora_read_tsv_character(gmt_paths[[5L]])
  stopifnot(
    identical(names(saved_results), names(gmt_run$results)),
    identical(saved_results$pathway, c("A", "B")),
    identical(saved_audit$description, gmt_run$pathway_audit$description),
    identical(saved_audit$test_status, gmt_run$pathway_audit$test_status),
    saved_parameters$value[saved_parameters$key == "resolved_engine"] == "base::phyper",
    length(output_paths) == 5L,
    all(file.exists(output_paths)),
    sum(list_run$gene_audit$raw_query_occurrence_count) == 1L,
    sum(list_run$gene_audit$raw_background_occurrence_count) == 4L,
    nrow(list_run$gene_audit) == 4L,
    sum(list_run$gene_audit$in_query) == 1L,
    grepl("^[0-9a-f]{64}$", list_run$parameters$value[list_run$parameters$key == "background_unique_sha256"])
  )
  expect_error(write_ora_results(list_run, output_prefix))
  rollback_prefix <- file.path(td, "rollback_test")
  old_option <- getOption("ora.test_fail_after_commit")
  options(ora.test_fail_after_commit = 1L)
  expect_error(write_ora_results(list_run, rollback_prefix))
  if (is.null(old_option)) options(ora.test_fail_after_commit = NULL) else options(ora.test_fail_after_commit = old_option)
  stopifnot(!any(file.exists(paste0(rollback_prefix, c("_ora.tsv", "_coverage.tsv", "_pathway_audit.tsv", "_gene_audit.tsv", "_parameters.tsv")))))

  if (requireNamespace("fgsea", quietly = TRUE)) {
    run_ora(sets, query, background, min_size = 2L, max_size = 5L, engine = "fgsea")
    message("fgsea parity test: PASS")
  } else {
    message("fgsea parity test: SKIP (fgsea is not installed in active .libPaths())")
  }
  message("ORA synthetic tests: PASS")
  invisible(TRUE)
}
