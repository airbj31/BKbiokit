#' Preranked GSEA with per-subset contribution and NES
#'
#' Runs `fgsea::fgsea()` on `pathways` against `stats`, then, for every
#' pathway and every named gene subset in `subsets`, reports
#'
#'   1. how much of the pathway's enrichment is carried by the subset genes
#'      (proportions of pathway size, leading-edge size and leading-edge
#'      weight), and
#'   2. an independent ES / NES / p-value for the subset genes that fall in
#'      the pathway, treating `intersect(pathway, subset)` as its own gene set
#'      against the same ranked background.
#'
#' @param pathways Named list of gene sets (character vectors), as in fgsea.
#' @param stats Named numeric vector of gene-level statistics. The names are
#'   the background (universe).
#' @param subsets Named list of character vectors, e.g.
#'   `list(subsetA = c("123","234"), subsetB = c("134","156"))`.
#'   Every gene must be present in `names(stats)` (see `strictSubsets`).
#'   Subsets may overlap; in that case the per-pathway proportions of the
#'   subsets can sum to more than 1.
#' @param minSize,maxSize Pathway size filters passed to `fgsea::fgsea()`.
#' @param subsetMinSize Minimal size of `intersect(pathway, subset)` for which
#'   a subset-level ES/NES is computed. Smaller overlaps still get the
#'   proportion columns but NA for the subset enrichment columns.
#' @param gseaParam GSEA weight parameter (passed to fgsea and used here).
#' @param scoreType "std", "pos" or "neg" (passed to fgsea and used here).
#' @param strictSubsets If TRUE (default) an error is raised when a subset
#'   contains genes missing from `names(stats)`. If FALSE, those genes are
#'   dropped with a warning.
#' @param ... Further arguments for `fgsea::fgsea()` (e.g. `eps`, `nproc`,
#'   `sampleSize`, `nPermSimple`, `BPPARAM`, or `nperm` to force
#'   `fgseaSimple`). They are used for both the pathway-level and the
#'   subset-level runs.
#'
#' @return A `data.table` in long format with one row per
#'   (pathway, subset) pair. Pathway-level columns are the usual fgsea
#'   columns (`pathway`, `pval`, `padj`, `log2err`, `ES`, `NES`, `size`,
#'   `leadingEdge`), repeated for each subset. Subset-level columns:
#'   \itemize{
#'     \item `subset` -- name of the subset;
#'     \item `subsetSize` -- number of subset genes in the pathway
#'       (after restricting to the background);
#'     \item `subsetLESize` -- number of subset genes in the pathway's
#'       leading edge;
#'     \item `subsetPropSize` -- `subsetSize / size`;
#'     \item `subsetPropLE` -- `subsetLESize / length(leadingEdge)`;
#'     \item `subsetPropWeight` -- share of the leading-edge hit weight
#'       (sum of |stat|^gseaParam over leading-edge genes) that comes from
#'       subset genes. This is the weighted "effect in proportion";
#'     \item `subsetLeadingEdge` -- subset genes in the leading edge (list);
#'     \item `subsetES`, `subsetNES`, `subsetPval`, `subsetPadj`,
#'       `subsetLog2err` (or `subsetNMoreExtreme` when `nperm` is used) --
#'       fgsea results for `intersect(pathway, subset)` as a gene set.
#'       `subsetPadj` is BH-adjusted across all (pathway, subset) tests.
#'   }
#'   Rows are ordered by pathway then subset. Use `data.table::dcast()`
#'   for a wide layout.
#'
#' @examples
#' library(fgsea)
#' data(examplePathways); data(exampleRanks)
#' subs <- list(A = names(exampleRanks)[1:200], B = names(exampleRanks)[201:600])
#' res <- fgseaSubset(examplePathways[1:20], exampleRanks, subs, minSize = 15, maxSize = 500)
#' @export
fgseaSubset <- function(pathways,
                        stats,
                        subsets,
                        minSize       = 1,
                        maxSize       = length(stats) - 1,
                        subsetMinSize = 1,
                        gseaParam     = 1,
                        scoreType     = c("std", "pos", "neg"),
                        strictSubsets = TRUE,
                        ...) {
  if (!requireNamespace("fgsea", quietly = TRUE)) stop("Package 'fgsea' is required")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required")
  scoreType <- match.arg(scoreType)

  ## ---- validate inputs ---------------------------------------------------
  if (!is.list(pathways) || is.null(names(pathways)) || any(names(pathways) == ""))
    stop("`pathways` must be a named list of gene identifier vectors")
  if (anyDuplicated(names(pathways)))
    stop("Duplicated pathway names are not allowed")
  if (is.null(names(stats))) stop("`stats` must be a named vector")
  universe <- names(stats)

  subsets <- .checkSubsets(subsets, universe, strictSubsets)

  ## ---- 1. pathway-level fgsea --------------------------------------------
  mainRes <- fgsea::fgsea(pathways = pathways, stats = stats,
                          minSize = minSize, maxSize = maxSize,
                          gseaParam = gseaParam, scoreType = scoreType, ...)
  mainCols <- names(mainRes)
  if (nrow(mainRes) == 0) return(.emptySubsetResult(mainRes))

  ## ---- 2. per-(pathway, subset) proportions -------------------------------
  # Same ordering / weighting as fgsea uses internally (prepareStats):
  ranks <- sort(stats, decreasing = TRUE)
  rankNames <- names(ranks)
  rAdj <- abs(ranks) ^ gseaParam

  subsetIdx <- lapply(subsets, function(g) sort(fastmatch::fmatch(g, rankNames)))

  propList <- lapply(seq_len(nrow(mainRes)), function(i) {
    pw <- mainRes$pathway[i]
    S  <- unique(stats::na.omit(fastmatch::fmatch(pathways[[pw]], rankNames)))
    # Leading edge exactly as fgsea reported it (already in leading-edge order)
    LE <- fastmatch::fmatch(mainRes$leadingEdge[[i]], rankNames)
    leWeight <- sum(rAdj[LE])
    data.table::rbindlist(lapply(names(subsets), function(sn) {
      ov   <- intersect(S, subsetIdx[[sn]])
      ovLE <- intersect(LE, subsetIdx[[sn]])
      data.table::data.table(
        pathway           = pw,
        subset            = sn,
        subsetSize        = length(ov),
        subsetLESize      = length(ovLE),
        subsetPropSize    = if (length(S)) length(ov) / length(S) else NA_real_,
        subsetPropLE      = if (length(LE)) length(ovLE) / length(LE) else NA_real_,
        subsetPropWeight  = if (length(LE) && leWeight > 0) sum(rAdj[ovLE]) / leWeight else NA_real_,
        subsetLeadingEdge = list(rankNames[LE[LE %in% ovLE]]) # keep leading-edge order
      )
    }))
  })
  propTab <- data.table::rbindlist(propList)

  ## ---- 3. subset-level fgsea (intersect(pathway, subset) as gene set) -----
  keys <- data.table::CJ(pathway = mainRes$pathway, subset = names(subsets), sorted = FALSE)
  keys[, key := sprintf("k%06d", .I)]
  subPathways <- lapply(seq_len(nrow(keys)), function(i) {
    g <- intersect(intersect(pathways[[keys$pathway[i]]], universe), subsets[[keys$subset[i]]])
    if (length(g) >= max(1L, subsetMinSize)) g else NULL
  })
  names(subPathways) <- keys$key
  subPathways <- Filter(Negate(is.null), subPathways)

  # fgsea columns to carry over with a "subset" prefix (log2err for multilevel, nMoreExtreme for nperm)
  statCols <- setdiff(mainCols, c("pathway", "size", "leadingEdge"))
  prefixed <- function(x) paste0("subset", toupper(substr(x, 1, 1)), substring(x, 2))

  if (length(subPathways) > 0) {
    subRes <- fgsea::fgsea(pathways = subPathways, stats = stats,
                           minSize = max(1L, subsetMinSize), maxSize = length(stats) - 1,
                           gseaParam = gseaParam, scoreType = scoreType, ...)
    subRes <- data.table::as.data.table(subRes)
    data.table::setnames(subRes, "pathway", "key")
    subRes <- merge(keys, subRes, by = "key", all.y = TRUE, sort = FALSE)
    subRes[, c("key", "size", "leadingEdge") := NULL]
    old <- setdiff(names(subRes), c("pathway", "subset"))
    data.table::setnames(subRes, old, prefixed(old))
    propTab <- merge(propTab, subRes, by = c("pathway", "subset"), all.x = TRUE, sort = FALSE)
  } else {
    for (cc in prefixed(statCols)) propTab[, (cc) := NA_real_]
  }

  ## ---- 4. combine ----------------------------------------------------------
  out <- merge(data.table::as.data.table(mainRes), propTab, by = "pathway", all.x = TRUE, sort = FALSE)
  out[, subset := factor(subset, levels = names(subsets))]
  data.table::setorder(out, pathway, subset)
  out[, subset := as.character(subset)]
  data.table::setcolorder(out, c(mainCols, setdiff(names(out), mainCols)))
  out[]
}

## Validate the `subsets` argument. Returns the (possibly cleaned) list.
.checkSubsets <- function(subsets, universe, strict = TRUE) {
  if (!is.list(subsets) || length(subsets) == 0)
    stop("`subsets` must be a non-empty named list of gene identifier vectors, ",
         "e.g. list(subsetA = c(\"123\",\"234\"), subsetB = c(\"134\",\"156\"))")
  nms <- names(subsets)
  if (is.null(nms) || any(is.na(nms)) || any(nms == ""))
    stop("Every element of `subsets` must be named")
  if (any(duplicated(nms)))
    stop("Duplicated subset names: ", paste(unique(nms[duplicated(nms)]), collapse = ", "))

  for (sn in nms) {
    g <- subsets[[sn]]
    if (!is.atomic(g) || is.null(g) || length(g) == 0)
      stop("Subset '", sn, "' is empty")
    if (is.numeric(g)) g <- format(g, scientific = FALSE, trim = TRUE)  # avoid "1e+05"
    g <- as.character(g)
    if (any(is.na(g)) || any(g == ""))
      stop("Subset '", sn, "' contains NA or empty gene identifiers")
    g <- unique(g)
    missing <- setdiff(g, universe)
    if (length(missing) > 0) {
      msg <- sprintf("Subset '%s' has %d gene(s) not in the background (names(stats)): %s%s",
                     sn, length(missing),
                     paste(head(missing, 10), collapse = ", "),
                     if (length(missing) > 10) ", ..." else "")
      if (strict) stop(msg, "\nAll subset genes must be in names(stats). ",
                       "Set strictSubsets = FALSE to drop them with a warning.")
      warning(msg, " -- dropping them.")
      g <- setdiff(g, missing)
      if (length(g) == 0) stop("Subset '", sn, "' has no genes left in the background")
    }
    subsets[[sn]] <- g
  }
  subsets
}

## Zero-row result with the full column set.
.emptySubsetResult <- function(mainRes) {
  out <- data.table::as.data.table(mainRes)
  out[, `:=`(subset = character(), subsetSize = integer(), subsetLESize = integer(),
             subsetPropSize = numeric(), subsetPropLE = numeric(), subsetPropWeight = numeric(),
             subsetLeadingEdge = list())]
  statCols <- setdiff(names(mainRes), c("pathway", "size", "leadingEdge"))
  for (cc in paste0("subset", toupper(substr(statCols, 1, 1)), substring(statCols, 2))) out[, (cc) := numeric()]
  out[]
}
