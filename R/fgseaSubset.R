#' Preranked GSEA with per-subset share of each pathway's NES
#'
#' Runs `fgsea::fgsea()` ONCE on `pathways` against `stats`. Then, for every
#' pathway and every named gene subset in `subsets`, it attributes a share of
#' the pathway's enrichment to the subset genes. Nothing is re-tested: the
#' pathway NES is computed a single time and split among the subsets.
#'
#' How the split works: the GSEA enrichment score is the peak of a running
#' sum in which every leading-edge gene adds `|stat|^gseaParam`. The share of
#' that leading-edge weight carried by a subset's genes is the subset's share
#' of the enrichment. For a partition of the pathway (e.g. a subset and its
#' complement) the shares add up to 100 % and the NES contributions add up
#' to the pathway NES.
#'
#' A direction-aware variant (`pctNESdir`) additionally signs each
#' leading-edge gene's weight by whether the gene moves with or against the
#' pathway direction.
#'
#' @param pathways Named list of gene sets (character vectors), as in fgsea.
#' @param stats Named numeric vector of gene-level statistics. The names are
#'   the background (universe).
#' @param subsets Named list of character vectors, e.g.
#'   `list(subsetA = c("123","234"), subsetB = c("134","156"))`.
#'   Every gene must be present in `names(stats)` (see `strictSubsets`).
#'   Subsets may overlap; in that case shares can sum to more than 100 %.
#' @param minSize,maxSize Pathway size filters passed to `fgsea::fgsea()`.
#' @param gseaParam GSEA weight parameter (passed to fgsea and used here).
#' @param scoreType "std", "pos" or "neg" (passed to fgsea).
#' @param format "long" (default): one row per (pathway, subset), with a
#'   `subset` column. "wide": one row per pathway, subset columns prefixed
#'   with the subset name, e.g. `subsetA_pctNES`.
#' @param subsetTest If TRUE, additionally run an independent fgsea test on
#'   `intersect(pathway, subset)` as its own gene set (columns prefixed
#'   `indep`). Off by default because it costs one extra fgsea run.
#' @param subsetMinSize Minimal overlap size for the independent test.
#' @param strictSubsets If TRUE (default) an error is raised when a subset
#'   contains genes missing from `names(stats)`. If FALSE, those genes are
#'   dropped with a warning.
#' @param ... Further arguments for `fgsea::fgsea()` (e.g. `eps`, `nproc`,
#'   `nPermSimple`, `BPPARAM`, or `nperm` to force `fgseaSimple`).
#'
#' @return A `data.table`. Pathway-level columns are the usual fgsea columns
#'   (`pathway`, `pval`, `padj`, `log2err`, `ES`, `NES`, `size`,
#'   `leadingEdge`). Per-subset columns (long format; in wide format each is
#'   prefixed with `<subsetName>_` and the `subset` column is dropped):
#'   \itemize{
#'     \item `subset` -- name of the subset;
#'     \item `nGenes` -- number of subset genes in the pathway;
#'     \item `pctSize` -- `100 * nGenes / size`;
#'     \item `nLE` -- number of subset genes in the pathway's leading edge;
#'     \item `pctLE` -- `100 * nLE / length(leadingEdge)`;
#'     \item `pctNES` -- percentage of the pathway's enrichment carried by the
#'       subset: share of leading-edge weight `sum(|stat|^gseaParam)`;
#'     \item `NEScontrib` -- `NES * pctNES / 100`, the subset's contribution
#'       to the pathway NES (additive over a partition of the pathway);
#'     \item `EScontrib` -- `ES * pctNES / 100`;
#'     \item `nLEopp` -- number of subset leading-edge genes whose statistic
#'       has the opposite sign to the pathway's ES (moving against the
#'       pathway direction);
#'     \item `pctNESdir` -- direction-aware share. Each leading-edge gene's
#'       weight is signed: `+|stat|^gseaParam` if the gene moves in the
#'       pathway direction, `-|stat|^gseaParam` if against it. `pctNESdir`
#'       is the sum of the subset's signed weights over the total absolute
#'       leading-edge weight, times 100. Ranges from -100 to 100 and equals
#'       `pctNES` when no subset gene opposes the pathway direction;
#'     \item `NEScontribDir` -- `NES * pctNESdir / 100`;
#'     \item `meanStat` -- mean of the raw ranking statistic over the subset
#'       genes in the pathway (all of them, not only the leading edge). Its
#'       sign is the subset's own direction, its magnitude the amount;
#'     \item `genes` -- subset genes in the pathway (list column);
#'     \item `leGenes` -- subset genes in the leading edge (list column).
#'     \item with `subsetTest = TRUE`: `indepES`, `indepNES`, `indepPval`,
#'       `indepPadj`, `indepLog2err` (or `indepNMoreExtreme` with `nperm`)
#'       for `intersect(pathway, subset)` tested as its own gene set.
#'   }
#'
#' @examples
#' library(fgsea)
#' data(examplePathways); data(exampleRanks)
#' subs <- list(A = names(exampleRanks)[1:200], B = names(exampleRanks)[201:600])
#' fgseaSubset(examplePathways[1:20], exampleRanks, subs, minSize = 15, maxSize = 500)
#' fgseaSubset(examplePathways[1:20], exampleRanks, subs, minSize = 15, maxSize = 500,
#'             format = "wide")
#' @export
fgseaSubset <- function(pathways,
                        stats,
                        subsets,
                        minSize       = 1,
                        maxSize       = length(stats) - 1,
                        gseaParam     = 1,
                        scoreType     = c("std", "pos", "neg"),
                        format        = c("long", "wide"),
                        subsetTest    = FALSE,
                        subsetMinSize = 1,
                        strictSubsets = TRUE,
                        ...) {
    if (!requireNamespace("fgsea", quietly = TRUE)) stop("Package 'fgsea' is required")
    if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required")
    scoreType <- match.arg(scoreType)
    format    <- match.arg(format)

    ## ---- validate inputs ---------------------------------------------------
    if (!is.list(pathways) || is.null(names(pathways)) || any(names(pathways) == ""))
        stop("`pathways` must be a named list of gene identifier vectors")
    if (anyDuplicated(names(pathways)))
        stop("Duplicated pathway names are not allowed")
    if (is.null(names(stats))) stop("`stats` must be a named vector")
    universe <- names(stats)

    subsets <- .checkSubsets(subsets, universe, strictSubsets)
    subsetNames <- names(subsets)

    ## ---- 1. pathway-level fgsea: the ONLY NES computation --------------------
    mainRes  <- data.table::as.data.table(
        fgsea::fgsea(pathways = pathways, stats = stats,
                     minSize = minSize, maxSize = maxSize,
                     gseaParam = gseaParam, scoreType = scoreType, ...))
    mainCols <- names(mainRes)
    statCols <- setdiff(mainCols, c("pathway", "size", "leadingEdge"))
    indepCols <- if (subsetTest) .indepNames(statCols) else character()

    if (nrow(mainRes) == 0)
        return(.emptySubsetResult(mainRes, indepCols, format, subsetNames))

    ## ---- 2. split each pathway's enrichment among the subsets ----------------
    # Same weighting as fgsea's running sum: |stat|^gseaParam per leading-edge gene.
    ranks     <- sort(stats, decreasing = TRUE)
    rankNames <- names(ranks)
    weight    <- abs(ranks) ^ gseaParam
    subsetIdx <- lapply(subsets, function(g) fastmatch::fmatch(g, rankNames))

    shareTab <- data.table::rbindlist(lapply(seq_len(nrow(mainRes)), function(i) {
        pw <- mainRes$pathway[i]
        S  <- unique(stats::na.omit(fastmatch::fmatch(pathways[[pw]], rankNames)))
        LE <- fastmatch::fmatch(mainRes$leadingEdge[[i]], rankNames)  # as fgsea reported it
        leWeight <- sum(weight[LE])
        # Pathway direction: +1 up-regulated, -1 down-regulated (from ES; NES has the same sign)
        dir <- sign(mainRes$ES[i])
        # Signed weight: +|stat|^p when the gene moves with the pathway, -|stat|^p against it
        aligned <- dir * sign(ranks) * weight
        data.table::rbindlist(lapply(subsetNames, function(sn) {
            inPw <- S[S %in% subsetIdx[[sn]]]
            inLE <- LE[LE %in% subsetIdx[[sn]]]           # keeps leading-edge order
            pct  <- if (length(LE) && leWeight > 0) 100 * sum(weight[inLE])  / leWeight else NA_real_
            pctD <- if (length(LE) && leWeight > 0) 100 * sum(aligned[inLE]) / leWeight else NA_real_
            data.table::data.table(
                pathway       = pw,
                subset        = sn,
                nGenes        = length(inPw),
                pctSize       = if (length(S)) 100 * length(inPw) / length(S) else NA_real_,
                nLE           = length(inLE),
                pctLE         = if (length(LE)) 100 * length(inLE) / length(LE) else NA_real_,
                pctNES        = pct,
                NEScontrib    = mainRes$NES[i] * pct / 100,
                EScontrib     = mainRes$ES[i]  * pct / 100,
                nLEopp        = sum(sign(ranks[inLE]) == -dir),
                pctNESdir     = pctD,
                NEScontribDir = mainRes$NES[i] * pctD / 100,
                meanStat      = if (length(inPw)) mean(ranks[inPw]) else NA_real_,
                genes         = list(rankNames[inPw]),
                leGenes       = list(rankNames[inLE]))
        }))
    }))

    ## ---- 3. optional independent test of intersect(pathway, subset) --------
    if (subsetTest) {
        indep <- .independentSubsetTest(mainRes$pathway, pathways, universe, subsets, stats,
                                        subsetMinSize, gseaParam, scoreType, statCols, ...)
        shareTab <- merge(shareTab, indep, by = c("pathway", "subset"), all.x = TRUE, sort = FALSE)
    }

    ## ---- 4. combine ----------------------------------------------------------
    out <- merge(mainRes, shareTab, by = "pathway", all.x = TRUE, sort = FALSE)
    # Order: pathway, then subsets in the order given by the user.
    # data.table::set()/setorderv() are used (not `[` with `:=`) so the function
    # also works when sourced from an environment that is not data.table-aware.
    data.table::set(out, j = ".subsetOrd", value = match(out$subset, subsetNames))
    data.table::setorderv(out, c("pathway", ".subsetOrd"))
    data.table::set(out, j = ".subsetOrd", value = NULL)
    data.table::setcolorder(out, c(mainCols, setdiff(names(out), mainCols)))

    if (format == "wide") out <- .toWide(out, mainCols, subsetNames)
    out
}

## Names of the independent-test columns derived from the fgsea columns present.
.indepNames <- function(statCols) paste0("indep", toupper(substr(statCols, 1, 1)), substring(statCols, 2))

## Run fgsea on every intersect(pathway, subset) as its own gene set.
.independentSubsetTest <- function(pathwayNames, pathways, universe, subsets, stats,
                                   subsetMinSize, gseaParam, scoreType, statCols, ...) {
    keys <- data.table::CJ(pathway = pathwayNames, subset = names(subsets), sorted = FALSE)
    data.table::set(keys, j = "key", value = sprintf("k%06d", seq_len(nrow(keys))))
    subPathways <- lapply(seq_len(nrow(keys)), function(i) {
        g <- intersect(intersect(pathways[[keys$pathway[i]]], universe), subsets[[keys$subset[i]]])
        if (length(g) >= max(1L, subsetMinSize)) g else NULL
    })
    names(subPathways) <- keys$key
    subPathways <- Filter(Negate(is.null), subPathways)

    pairs <- data.table::data.table(pathway = keys$pathway, subset = keys$subset)
    empty <- data.table::copy(pairs)
    for (cc in .indepNames(statCols)) data.table::set(empty, j = cc, value = NA_real_)
    if (length(subPathways) == 0) return(empty)

    res <- data.table::as.data.table(
        fgsea::fgsea(pathways = subPathways, stats = stats,
                     minSize = max(1L, subsetMinSize), maxSize = length(stats) - 1,
                     gseaParam = gseaParam, scoreType = scoreType, ...))
    data.table::setnames(res, "pathway", "key")
    res <- merge(keys, res, by = "key", all.y = TRUE, sort = FALSE)
    data.table::set(res, j = c("key", "size", "leadingEdge"), value = NULL)
    old <- setdiff(names(res), c("pathway", "subset"))
    data.table::setnames(res, old, .indepNames(old))
    res <- merge(pairs, res, by = c("pathway", "subset"), all.x = TRUE, sort = FALSE)
    for (cc in .indepNames(statCols)) if (!cc %in% names(res)) data.table::set(res, j = cc, value = NA_real_)
    res
}

## Long -> wide: one row per pathway, subset columns prefixed with the subset name.
.toWide <- function(out, mainCols, subsetNames) {
    subCols <- setdiff(names(out), c(mainCols, "subset"))
    wide <- unique(.selectDT(out, seq_len(nrow(out)), mainCols), by = "pathway")
    for (sn in subsetNames) {
        part <- .selectDT(out, which(out$subset == sn), c("pathway", subCols))
        data.table::setnames(part, subCols, paste0(sn, "_", subCols))
        wide <- merge(wide, part, by = "pathway", all.x = TRUE, sort = FALSE)
    }
    wide
}

## Row/column subset of a data.table without `[`, so it works from any environment.
.selectDT <- function(dt, rows, cols) {
    data.table::setDT(lapply(as.list(dt)[cols], function(col) col[rows]))
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
.emptySubsetResult <- function(mainRes, indepCols, format, subsetNames) {
    out <- data.table::copy(mainRes)
    template <- list(subset = character(), nGenes = integer(), pctSize = numeric(),
                     nLE = integer(), pctLE = numeric(), pctNES = numeric(),
                     NEScontrib = numeric(), EScontrib = numeric(),
                     nLEopp = integer(), pctNESdir = numeric(), NEScontribDir = numeric(),
                     meanStat = numeric(), genes = list(), leGenes = list())
    for (cc in names(template)) data.table::set(out, j = cc, value = template[[cc]])
    for (cc in indepCols) data.table::set(out, j = cc, value = numeric())
    if (format == "wide") out <- .toWide(out, names(mainRes), subsetNames)
    out
}
