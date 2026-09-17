#' Detect and Correct Typos in Column Names
#'
#' Matches a given name against a vector of actual names using case-insensitive
#' exact matching first, then falls back to fuzzy matching via edit distance
#' (Levenshtein). Useful for forgiving column name lookups in data frames.
#'
#' @param x An unquoted name to search for (captured via \code{\link[rlang]{enquo}}).
#' @param y A character vector of valid names to search within (e.g. \code{colnames(df)}).
#' @param maxTypos An integer specifying the maximum edit distance allowed for
#'   fuzzy matching. Defaults to \code{3}.
#' @param silent Logical. If \code{TRUE}, suppresses the correction warning.
#'   Defaults to \code{FALSE}.
#'
#' @return A single character string from \code{y} that best matches \code{x}.
#'
#' @details
#' The function first attempts a case-insensitive exact match. If no exact match
#' is found, it computes the Levenshtein edit distance between \code{x} and every
#' element in \code{y} via \code{\link{adist}}, and accepts the closest match only
#' if its distance is within \code{maxTypos}. Unless \code{silent = TRUE}, a warning
#' is issued reporting the original input and its resolved name. An error is thrown
#' if no suitable match is found.
#'
#' @examples
#' cols <- c("PatientID", "SampleDate", "GeneExpression")
#'
#' # Exact match, different case
#' detect_typo(patientid, cols)
#'
#' # Fuzzy match with a small typo
#' detect_typo(GeneEpression, cols)
#'
#' # Suppress the correction warning
#' detect_typo(SmpleDate, cols, silent = TRUE)
#'
#' @seealso \code{\link{adist}} for the underlying string distance computation.
#'
#' @importFrom rlang enquo as_label
#'
#' @export
detect_typo <- function(x, y, maxTypos = 3, silent = FALSE) {
  CNAME <- gsub('"', '', rlang::as_label(rlang::enquo(x)))
  actual_cols <- y
  match_idx <- which(tolower(actual_cols) == tolower(CNAME))
  dst <- as.vector(adist(tolower(CNAME), tolower(actual_cols)))
  min_idx <- which.min(dst)
  if (dst[min_idx] <= maxTypos) match_idx <- min_idx

  if (length(match_idx) == 0) stop("we can't find ", CNAME, " from given vector.")
  if (!silent) message("[INFO] typo fix: ", CNAME, " -> ", actual_cols[match_idx])
  return(actual_cols[match_idx])
}
