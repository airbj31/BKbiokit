# NCBI Entrez ID ortholog mapping. Base R only; source() performs no I/O.
# https://ftp.ncbi.nlm.nih.gov/gene/DATA/README
# NCBI has no standalone orthology ID. Anchor keys below are provenance only.

#' Validate source and target taxonomy IDs
#'
#' Checks that two distinct NCBI taxonomy IDs are positive integer values and
#' converts them to character strings.
#'
#' @param from_taxid,to_taxid Single positive integer or digit string identifying
#'   the source and target species. The IDs must differ.
#' @return A character vector of length two, ordered as source then target.
#' @keywords internal
#' @noRd
.ortholog_taxa <- function(from_taxid, to_taxid) {
  check <- function(x) {
    if (length(x) != 1L || is.na(x) || !grepl("^[1-9][0-9]*$", as.character(x))) {
      stop("Taxonomy IDs must be single positive integers (e.g. 9606, 10090).")
    }
    as.character(x)
  }
  taxa <- c(check(from_taxid), check(to_taxid))
  if (taxa[1] == taxa[2]) stop("Choose two different species.")
  taxa
}

#' Read a species-pair reference from an NCBI ortholog file
#'
#' Reads a local NCBI `gene_orthologs` snapshot in chunks and returns unique
#' Entrez Gene ID pairs for two species. No network access is performed.
#'
#' @param path Path to an unmodified NCBI `gene_orthologs` tab-delimited file.
#'   Files whose names end in `.gz` are read as gzip-compressed text.
#' @param from_taxid,to_taxid Single positive integer or digit string containing
#'   an NCBI taxonomy ID. For example, human is `9606`, mouse is `10090`, and
#'   rat is `10116`. Source and target must differ; species names are not accepted.
#' @param method Relationship selection method. `"shared_anchor"` (default)
#'   includes direct pairs and pairs sharing the same NCBI primary gene.
#'   `"direct"` includes only pairs explicitly recorded in the file, accepting
#'   either orientation.
#'
#' @details
#' NCBI records a primary gene and the other members of its ortholog set rather
#' than all possible member pairs. With `method = "shared_anchor"`, two genes
#' are linked if they occur in the same primary-gene set. For example, human H
#' linked to mouse M and rat R permits a mouse M to rat R mapping. The function
#' does not follow arbitrary chains or compute a transitive closure.
#'
#' NCBI does not supply a standalone orthology/group ID in this file.
#' `anchor_keys` are provenance keys, not official orthology identifiers.
#' Coverage depends on the species and snapshot; an absent mapping is not proof
#' that a biological ortholog does not exist. Invalid headers or malformed
#' retained records cause an error.
#'
#' @return A data frame with one row per distinct source/target pair, sorted
#'   lexicographically by source and target ID. An absent species pair returns
#'   a zero-row data frame. Columns are:
#'   * `source_entrezid`, `target_entrezid`: character Entrez Gene IDs.
#'   * `evidence`: `"direct"` or `"shared_anchor"`; direct evidence takes
#'     precedence if both are available.
#'   * `anchor_keys`: semicolon-separated primary `tax_id:GeneID` keys;
#'     `NA` when `method = "direct"`.
#'
#'   Attributes record `from_taxid`, `to_taxid`, `method`, `source`, `source_url`,
#'   `source_md5` (the original file's MD5), and `processed_utc`. Processing time
#'   is not the database release date. Use `saveRDS()` to preserve attributes.
#' @seealso [load_orthologs()], [find_orthologs()], [hs2mm()]
#' @references
#'   <https://ftp.ncbi.nlm.nih.gov/gene/DATA/README>
#' @examples
#' # A small synthetic file demonstrates a shared primary gene without a download.
#' path <- tempfile(fileext = ".tsv")
#' writeLines(c(
#'   "#tax_id\tGeneID\trelationship\tOther_tax_id\tOther_GeneID",
#'   "9606\t1\tOrtholog\t10090\t11",
#'   "9606\t1\tOrtholog\t10116\t21"
#' ), path)
#' read_orthologs(path, from_taxid = 10090, to_taxid = 10116)
#' read_orthologs(path, 10090, 10116, method = "direct") # No direct pair
#' unlink(path)
#' @md
#' @export
read_orthologs <- function(path, from_taxid, to_taxid,
                           method = c("shared_anchor", "direct")) {
  taxa <- .ortholog_taxa(from_taxid, to_taxid)
  method <- match.arg(method)
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  header <- readLines(con, n = 1L)
  expected <- c("#tax_id", "GeneID", "relationship", "Other_tax_id", "Other_GeneID")
  if (length(header) != 1L ||
      !identical(strsplit(header, "\t", fixed = TRUE)[[1]], expected)) {
    stop("Unexpected gene_orthologs header.")
  }
  chunks <- list()
  # Only rows whose member belongs to either requested species can contribute.
  # A requested primary gene is included automatically with its relevant member.
  pattern <- paste0("^[^\t]+\t[^\t]+\tOrtholog\t(",
                    paste(taxa, collapse = "|"), ")\t")
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[grepl(pattern, lines)]
    if (!length(lines)) next
    x <- utils::read.delim(text = lines, header = FALSE, sep = "\t", quote = "",
                           comment.char = "", colClasses = "character")
    if (ncol(x) != 5L || anyNA(x) || any(x[[3]] != "Ortholog")) {
      stop("Malformed ortholog record.")
    }
    for (j in c(1L, 2L, 4L, 5L)) {
      if (any(!grepl("^[1-9][0-9]*$", x[[j]]))) stop("Malformed taxonomy/GeneID.")
    }
    chunks[[length(chunks) + 1L]] <- x
  }
  out <- data.frame(source_entrezid = character(), target_entrezid = character(),
                    evidence = character(), anchor_keys = character())
  if (length(chunks)) {
    x <- unique(do.call(rbind, chunks))
    names(x) <- c("tax", "gene", "relationship", "other_tax", "other_gene")
    forward <- x$tax == taxa[1] & x$other_tax == taxa[2]
    reverse <- x$tax == taxa[2] & x$other_tax == taxa[1]
    direct <- unique(rbind(
      data.frame(source_entrezid = x$gene[forward], target_entrezid = x$other_gene[forward]),
      data.frame(source_entrezid = x$other_gene[reverse], target_entrezid = x$gene[reverse])))
    if (method == "direct") {
      out <- direct
      out$evidence <- rep("direct", nrow(out))
      out$anchor_keys <- rep(NA_character_, nrow(out))
    } else {
      # NCBI stores a star: primary gene -> each member, not every member pair.
      # Join members of the SAME primary-gene set; do not traverse arbitrary chains.
      anchor <- paste(x$tax, x$gene, sep = ":")
      members <- unique(rbind(
        data.frame(anchor = anchor, tax = x$tax, gene = x$gene),
        data.frame(anchor = anchor, tax = x$other_tax, gene = x$other_gene)))
      a <- members[members$tax == taxa[1], c("anchor", "gene")]
      b <- members[members$tax == taxa[2], c("anchor", "gene")]
      joined <- merge(a, b, by = "anchor", suffixes = c("_source", "_target"))
      if (nrow(joined)) {
        out <- aggregate(anchor ~ gene_source + gene_target, joined,
                         function(z) paste(sort(unique(z)), collapse = ";"))
        names(out) <- c("source_entrezid", "target_entrezid", "anchor_keys")
        key <- function(d) paste(d$source_entrezid, d$target_entrezid, sep = ":")
        out$evidence <- ifelse(key(out) %in% key(direct), "direct", "shared_anchor")
      }
    }
  }
  out <- out[order(out$source_entrezid, out$target_entrezid), ]
  rownames(out) <- NULL
  attr(out, "from_taxid") <- taxa[1]
  attr(out, "to_taxid") <- taxa[2]
  attr(out, "method") <- method
  attr(out, "source") <- "NCBI gene_orthologs"
  attr(out, "source_url") <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz"
  attr(out, "source_md5") <- unname(tools::md5sum(path))
  attr(out, "processed_utc") <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  out
}

#' Download or reuse an NCBI ortholog reference
#'
#' Loads a species-pair reference from a local file or a cached NCBI snapshot.
#' A downloaded compressed snapshot is retained for reuse across species pairs.
#'
#' @inheritParams read_orthologs
#' @param cache_dir Directory for the compressed NCBI snapshot and filtered RDS
#'   caches. Relative paths are resolved against the current working directory.
#' @param refresh Single logical value. If `TRUE`, download a new snapshot and
#'   rebuild the requested pair cache. If `FALSE`, reuse the snapshot when
#'   present. Ignored when `source_file` is supplied.
#' @param source_file Optional path to a local NCBI file accepted by
#'   [read_orthologs()]. When supplied, bypass downloading and caching.
#'
#' @details
#' The download URL is
#' <https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz>.
#' The raw snapshot is stored as `gene_orthologs.gz` in `cache_dir`. Filtered
#' caches are keyed by source/target taxonomy IDs, method, and the raw-file MD5,
#' preventing reuse of a pair cache from a different snapshot. There is no
#' automatic periodic refresh. Network, parsing, and file-writing errors are
#' propagated to the caller.
#'
#' @return The reference data frame documented in [read_orthologs()]. When using
#'   the managed cache, the additional `snapshot_file` attribute records the
#'   absolute path to the retained compressed snapshot.
#' @seealso [read_orthologs()], [find_orthologs()], [hs2mm()]
#' @examples
#' \dontrun{
#' # These calls may download the full NCBI snapshot and write cache files.
#' ref <- load_orthologs(9606, 10090)
#' mouse_rat <- load_orthologs(10090, 10116)
#' refreshed <- load_orthologs(9606, 10090, refresh = TRUE)
#' }
#' @md
#' @export
load_orthologs <- function(from_taxid, to_taxid,
                           cache_dir = "data/reference/ncbi_orthologs",
                           refresh = FALSE, source_file = NULL,
                           method = c("shared_anchor", "direct")) {
  taxa <- .ortholog_taxa(from_taxid, to_taxid)
  method <- match.arg(method)
  if (!is.logical(refresh) || length(refresh) != 1L || is.na(refresh)) stop("Invalid refresh.")
  if (!is.null(source_file)) return(read_orthologs(source_file, taxa[1], taxa[2], method))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  raw <- file.path(cache_dir, "gene_orthologs.gz")
  if (refresh || !file.exists(raw)) {
    pending <- tempfile("download_", tmpdir = cache_dir, fileext = ".gz")
    on.exit(unlink(pending), add = TRUE)
    old <- getOption("timeout")
    options(timeout = max(600, old))
    on.exit(options(timeout = old), add = TRUE)
    message("Downloading NCBI reference; retained for reuse across species.")
    utils::download.file("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz",
                         pending, mode = "wb", method = "libcurl")
    out <- read_orthologs(pending, taxa[1], taxa[2], method)
    if (!file.rename(pending, raw)) stop("Cannot save reference snapshot.")
  }
  # The snapshot hash prevents reusing a different species pair's stale cache.
  md5 <- unname(tools::md5sum(raw))
  cache <- file.path(cache_dir, paste0("v1_", taxa[1], "_", taxa[2], "_", method, "_", md5, ".rds"))
  if (file.exists(cache) && !refresh) return(readRDS(cache))
  if (!exists("out", inherits = FALSE)) out <- read_orthologs(raw, taxa[1], taxa[2], method)
  attr(out, "snapshot_file") <- normalizePath(raw)
  pending_cache <- tempfile("cache_", tmpdir = cache_dir)
  on.exit(unlink(pending_cache), add = TRUE)
  saveRDS(out, pending_cache)
  if (!file.rename(pending_cache, cache)) stop("Cannot save pair cache.")
  out
}

# Multiple target candidates return NA by default; details=TRUE keeps all.
# Cardinality is based on the full supplied reference, not the input subset.
#' Find Entrez Gene ID orthologs between two species
#'
#' Maps source-species Entrez Gene IDs to target-species IDs using an NCBI
#' reference or a supplied mapping table. Input order and duplicates are retained.
#'
#' @param entrezid A character, numeric, or factor vector of source Entrez Gene
#'   IDs. Character input is recommended. Leading/trailing whitespace is trimmed
#'   for lookup, and only positive integer strings are accepted after conversion
#'   to character. Gene symbols, prefixed IDs, decimal strings, and obsolete IDs
#'   are not automatically converted. `NA` and blank strings are missing inputs.
#' @inheritParams read_orthologs
#' @inheritParams load_orthologs
#' @param details Single logical value. If `FALSE` (default), return a named
#'   character vector. If `TRUE`, return all candidates and input statuses in a
#'   detailed data frame.
#' @param orthologs Optional reference data frame, typically produced by
#'   [load_orthologs()] or [read_orthologs()]. Required columns are
#'   `source_entrezid` and `target_entrezid`, containing nonmissing positive
#'   integer IDs. Optional `evidence` and `anchor_keys` columns are carried into
#'   detailed output. Loader-generated species attributes are checked against
#'   the requested direction. The caller is responsible for the species,
#'   completeness, and provenance of custom tables.
#'
#' @details
#' When `orthologs` is supplied, its existing relationships are used;
#' `cache_dir`, `refresh`, and `method` do not reload or filter it. Without a
#' supplied reference, the function calls [load_orthologs()] only when at least
#' one input is a valid ID. Sourcing the script itself performs no I/O.
#'
#' A single target candidate is returned even if another source gene also maps
#' to that target. The default output therefore does not enforce reciprocal
#' one-to-one mapping. To select reciprocal unique pairs, use `details = TRUE`
#' and filter `mapping_type == "one_to_one"`. Cardinality is calculated across
#' the entire supplied reference, not only the input subset, and is not a
#' biological confidence score. Missing reference pairs do not establish that
#' a biological ortholog is absent.
#'
#' @return With `details = FALSE`, a character vector of the same length as
#'   `entrezid`, named using the original inputs. Missing, invalid, unmapped,
#'   and multiple-candidate inputs return `NA`. Unmapped and multiple-candidate
#'   inputs produce a summary warning.
#'
#'   With `details = TRUE`, a data frame with one row per input/candidate and
#'   one row for each input without candidates. Columns are:
#'   * `input_row`: one-based input position.
#'   * `input_entrezid`: original input converted to character.
#'   * `source_entrezid`: trimmed lookup ID.
#'   * `target_entrezid`: candidate ID, or `NA`.
#'   * `relationship`: `"Ortholog"` for a mapped pair, otherwise `NA`.
#'   * `status`: `"mapped"`, `"ambiguous"`, `"unmapped"`, `"missing_id"`, or
#'     `"invalid_id"`. Ambiguous means multiple distinct target candidates.
#'   * `n_target_candidates`: number of distinct target candidates for the input.
#'   * `mapping_type`: `"one_to_one"`, `"one_to_many"`, `"many_to_one"`, or
#'     `"many_to_many"`, interpreted source-to-target; `NA` without candidates.
#'   * `evidence`, `anchor_keys`: provenance from the reference, or `NA`.
#'   * `source_taxid`, `target_taxid`: requested taxonomy IDs as strings.
#'
#'   The `reference_metadata` attribute preserves reference attributes other
#'   than its names, row names, and class. Empty input returns an empty vector
#'   or a zero-row detailed table.
#' @seealso [hs2mm()], [load_orthologs()], [read_orthologs()]
#' @examples
#' # Synthetic reference: no network access or cache writes.
#' ref <- data.frame(
#'   source_entrezid = c("1", "2", "2"),
#'   target_entrezid = c("11", "12", "13")
#' )
#' find_orthologs("1", 9606, 10090, orthologs = ref)
#' find_orthologs(c("1", "2", "99", NA), 9606, 10090,
#'                orthologs = ref, details = TRUE)
#' \dontrun{
#' # Real NCBI lookup; may download and cache a snapshot.
#' find_orthologs("22059", from_taxid = 10090, to_taxid = 10116,
#'                details = TRUE)
#' }
#' @md
#' @export
find_orthologs <- function(entrezid, from_taxid, to_taxid, details = FALSE,
                           cache_dir = "data/reference/ncbi_orthologs",
                           refresh = FALSE, orthologs = NULL,
                           method = c("shared_anchor", "direct")) {
  taxa <- .ortholog_taxa(from_taxid, to_taxid)
  method <- match.arg(method)
  if (!(is.character(entrezid) || is.numeric(entrezid) ||
        is.factor(entrezid)) || !is.null(dim(entrezid))) {
    stop("source_entrezid must be a character or numeric vector.")
  }
  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("details must be TRUE or FALSE.")
  }
  original <- as.character(entrezid)
  ids <- trimws(original)
  missing <- is.na(ids) | ids == ""
  valid <- !missing & grepl("^[1-9][0-9]*$", ids)
  if (is.null(orthologs)) {
    if (any(valid)) {
      orthologs <- load_orthologs(taxa[1], taxa[2], cache_dir, refresh, method = method)
    } else {
      orthologs <- data.frame(source_entrezid = character(), target_entrezid = character())
    }
  }
  if (!is.null(attr(orthologs, "from_taxid")) &&
      (!identical(attr(orthologs, "from_taxid"), taxa[1]) ||
       !identical(attr(orthologs, "to_taxid"), taxa[2]))) {
    stop("Reference species/direction does not match requested taxonomy IDs.")
  }
  required <- c("source_entrezid", "target_entrezid")
  if (!is.data.frame(orthologs) || !all(required %in% names(orthologs))) {
    stop("orthologs must contain source_entrezid and target_entrezid columns.")
  }
  pairs <- unique(data.frame(source_entrezid = as.character(orthologs$source_entrezid),
                             target_entrezid = as.character(orthologs$target_entrezid)))
  if (anyNA(pairs) || any(!grepl("^[1-9][0-9]*$", pairs$source_entrezid)) ||
      any(!grepl("^[1-9][0-9]*$", pairs$target_entrezid))) stop("Invalid reference IDs.")
  lookup <- lapply(split(pairs$target_entrezid, pairs$source_entrezid), sort)
  reverse_count <- table(pairs$target_entrezid)
  result <- rep(NA_character_, length(ids))
  rows <- vector("list", length(ids))
  ambiguous <- unmapped <- 0L
  for (i in seq_along(ids)) {
    hits <- if (valid[i]) lookup[[ids[i]]] else NULL
    n <- length(hits)
    status <- if (missing[i]) "missing_id" else if (!valid[i]) "invalid_id" else
      if (n == 0L) "unmapped" else if (n == 1L) "mapped" else "ambiguous"
    if (n == 1L) result[i] <- hits
    ambiguous <- ambiguous + (status == "ambiguous")
    unmapped <- unmapped + (status == "unmapped")
    if (!details) next
    types <- if (!n) NA_character_ else ifelse(
      reverse_count[hits] > 1L,
      if (n > 1L) "many_to_many" else "many_to_one",
      if (n > 1L) "one_to_many" else "one_to_one"
    )
    rows[[i]] <- data.frame(
      input_row = i, input_entrezid = original[i], source_entrezid = ids[i],
      target_entrezid = if (n) hits else NA_character_,
      relationship = if (n) "Ortholog" else NA_character_,
      status = status, n_target_candidates = n, mapping_type = unname(types),
      stringsAsFactors = FALSE
    )
  }
  if (details) {
    out <- if (length(rows)) do.call(rbind, rows) else data.frame(
      input_row = integer(), input_entrezid = character(), source_entrezid = character(),
      target_entrezid = character(), relationship = character(), status = character(),
      n_target_candidates = integer(), mapping_type = character()
    )
    pair_key <- function(a, b) paste(a, b, sep = ":")
    index <- match(pair_key(out$source_entrezid, out$target_entrezid),
                   pair_key(orthologs$source_entrezid, orthologs$target_entrezid))
    for (field in c("evidence", "anchor_keys")) {
      out[[field]] <- if (field %in% names(orthologs)) orthologs[[field]][index] else
        rep(NA_character_, nrow(out))
    }
    out$source_taxid <- rep(taxa[1], nrow(out))
    out$target_taxid <- rep(taxa[2], nrow(out))
    rownames(out) <- NULL
    attr(out, "reference_metadata") <- attributes(orthologs)[
      setdiff(names(attributes(orthologs)), c("names", "row.names", "class"))]
    return(out)
  }
  names(result) <- original
  if (ambiguous + unmapped > 0L) warning(
    unmapped, " unmapped and ", ambiguous,
    " ambiguous input row(s) returned as NA; use details = TRUE.", call. = FALSE)
  result
}

# Convenience wrapper: human (9606) -> mouse (10090).
#' Convert human Entrez Gene IDs to mouse ortholog IDs
#'
#' Convenience wrapper around [find_orthologs()] with human (`9606`) as the
#' source species and mouse (`10090`) as the target species.
#'
#' @inheritParams find_orthologs
#' @param ... Additional arguments passed to [find_orthologs()], such as
#'   `details`, `cache_dir`, `refresh`, `orthologs`, and `method`. Do not supply
#'   `from_taxid` or `to_taxid`; they are fixed by this wrapper.
#' @return The named character vector or detailed data frame described in
#'   [find_orthologs()]. Multiple mouse candidates yield `NA` by default;
#'   `details = TRUE` retains all candidates.
#' @seealso [find_orthologs()], [load_orthologs()], [read_orthologs()]
#' @examples
#' # Synthetic reference for an offline example.
#' ref <- data.frame(source_entrezid = "1", target_entrezid = "11")
#' hs2mm("1", orthologs = ref)
#' hs2mm(c("1", NA), orthologs = ref, details = TRUE)
#' \dontrun{
#' hs2mm(c("7157", "1956", "3569"))
#' hs2mm(c("7157", "1956", "3569"), details = TRUE)
#' }
#' @md
#' @export
hs2mm <- function(entrezid, ...) {
  find_orthologs(entrezid, from_taxid = 9606, to_taxid = 10090, ...)
}
