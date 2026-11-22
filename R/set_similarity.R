#' Pairwise set similarity: overlap coefficient or Jaccard index (Rcpp)
#'
#' Computes a similarity score for every pair of sets in a list.
#'
#' * `"overlap"`: overlap coefficient (Szymkiewicz-Simpson coefficient)
#'   \deqn{O(A, B) = \frac{|A \cap B|}{\min(|A|, |B|)}}
#' * `"jaccard"`: Jaccard index
#'   \deqn{J(A, B) = \frac{|A \cap B|}{|A \cup B|} = \frac{|A \cap B|}{|A| + |B| - |A \cap B|}}
#'
#' For each set, the shared-item counts with all later sets are accumulated
#' through an item-to-set inverted index, and the scores are written straight
#' into a dense symmetric matrix. No sparse matrix product and no dense
#' intermediate matrices are created, and rows can be processed in parallel
#' with OpenMP.
#'
#' @param list_data A list of atomic vectors, each representing one set.
#'   Names of the list are used as row and column names of the result.
#'   Duplicate elements and `NA` within a set are ignored.
#' @param method Similarity to compute: `"overlap"` (default) or `"jaccard"`.
#' @param nthreads Integer. Number of OpenMP threads. Ignored if the code was
#'   compiled without OpenMP support.
#'
#' @return A symmetric `length(list_data)` x `length(list_data)` numeric
#'   matrix of similarities in \[0, 1\]. The diagonal is 1 for non-empty
#'   sets. Pairs involving an empty set, and pairs with no shared elements,
#'   are 0.
#'
#' @seealso The overlap coefficient is always greater than or equal to the
#'   Jaccard index for the same pair; it equals 1 whenever one set is a subset
#'   of the other, while the Jaccard index penalises differences in set size.
#'
#' @export
#'
#' @examples
#' sets <- list(
#'   a = c("x", "y", "z"),
#'   b = c("y", "z"),
#'   c = c("w"),
#'   d = character(0)
#' )
#' set_similarity(sets)                     # overlap coefficient
#' set_similarity(sets, method = "jaccard") # Jaccard index
set_similarity <- function(list_data,
                           method = c("overlap", "jaccard"),
                           nthreads = 1L) {
  stopifnot(is.list(list_data))
  method <- match.arg(method)

  v <- unlist(list_data, use.names = FALSE)
  lev <- unique(v[!is.na(v)])  # NA is not an element; match() gives NA -> skipped

  res <- set_similarity_cpp(
    items    = match(v, lev),
    lens     = lengths(list_data),
    n_items  = length(lev),
    method   = if (method == "overlap") 0L else 1L,
    nthreads = as.integer(nthreads)
  )

  dimnames(res) <- list(names(list_data), names(list_data))
  res
}
