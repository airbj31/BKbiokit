# Official gene symbol -> ENTREZID using an installed Bioconductor OrgDb.
# Sourcing this file defines the function only; no downloads or installations.
# Install the required database separately, for example:
# BiocManager::install(c("AnnotationDbi", "org.Mm.eg.db"))
# Documentation: https://bioconductor.org/packages/AnnotationDbi
#
# Returns a named character vector by default, preserving input order, length,
# and duplicates. Missing/unmapped/ambiguous symbols yield NA, never a guessed ID.
# details = TRUE returns one audit row per input with all candidate IDs.
# Only exact, case-sensitive SYMBOL matches are used (no alias fallback).
# Leading/trailing whitespace is trimmed for lookup; original input is retained.
# Results depend on the installed annotation database version.
#
# Examples:
# dict <- Sym2EntrezId(c("Trp53", "Il6"), species = "mouse")
# audit <- Sym2EntrezId(c("Trp53", "", NA), species = "mouse", details = TRUE)
Sym2EntrezId <- function(OfficialGeneSymbol, species = c("mouse", "human", "rat"),
                         details = FALSE) {
  species <- match.arg(species)
  if (!is.character(OfficialGeneSymbol) && !is.factor(OfficialGeneSymbol)) {
    stop("OfficialGeneSymbol must be a character vector or factor.", call. = FALSE)
  }
  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("details must be TRUE or FALSE.", call. = FALSE)
  }
  package <- switch(species, mouse = "org.Mm.eg.db", human = "org.Hs.eg.db",
                    rat = "org.Rn.eg.db")
  for (dependency in c("AnnotationDbi", package)) {
    if (!requireNamespace(dependency, quietly = TRUE)) {
      stop("Missing package: ", dependency, ". Install with BiocManager::install(\"",
           dependency, "\").", call. = FALSE)
    }
  }
  db <- getExportedValue(package, package)
  original <- as.character(OfficialGeneSymbol)
  symbols <- trimws(original)
  valid <- !is.na(symbols) & nzchar(symbols)
  lookup <- intersect(unique(symbols[valid]), AnnotationDbi::keys(db, keytype = "SYMBOL"))
  candidates <- setNames(vector("list", length(lookup)), lookup)
  if (length(lookup)) {
    candidates <- AnnotationDbi::mapIds(
      db, keys = lookup, column = "ENTREZID", keytype = "SYMBOL", multiVals = "list"
    )
    candidates <- lapply(candidates, function(x) {
      x <- as.character(x)
      sort(unique(x[!is.na(x) & nzchar(x)]))
    })
  }
  ids <- rep(NA_character_, length(symbols))
  status <- ifelse(valid, "unmapped", "missing_symbol")
  candidate_text <- rep(NA_character_, length(symbols))
  for (i in which(valid)) {
    hits <- candidates[[symbols[[i]]]]
    if (!length(hits)) next
    candidate_text[[i]] <- paste(hits, collapse = ";")
    if (length(hits) == 1L) {
      ids[[i]] <- hits[[1L]]
      status[[i]] <- "mapped"
    } else {
      status[[i]] <- "ambiguous"
    }
  }
  if (details) {
    return(data.frame(
      OfficialGeneSymbol = original, lookup_symbol = symbols, ENTREZID = ids,
      status = status, candidate_ENTREZID = candidate_text,
      species = rep(species, length(symbols)),
      annotation_package = rep(package, length(symbols)),
      annotation_version = rep(as.character(utils::packageVersion(package)), length(symbols)),
      stringsAsFactors = FALSE
    ))
  }
  names(ids) <- original
  if (any(valid & status != "mapped")) {
    warning(sum(status == "unmapped"), " unmapped and ", sum(status == "ambiguous"),
            " ambiguous input row(s) returned as NA; use details = TRUE to inspect.",
            call. = FALSE)
  }
  ids
}
