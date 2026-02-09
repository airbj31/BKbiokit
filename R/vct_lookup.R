#' vector lookup / Dictionary Mapping
#'
#' Maps elements of a vector using a named vector (dictionary-style lookup).
#' Elements present in the dictionary keys are replaced with their corresponding values.
#' Elements not found in the dictionary remain unchanged.
#'
#' @param query A vector to be mapped (Left Hand Side for the operator).
#' @param ref A named vector serving as the lookup dictionary (Right Hand Side for the operator).
#'   The names of the vector are the keys, and the values are the replacements.
#'
#' @return A vector of the same class as \code{query} with mapped values.
#' @export
#'
#' @examples
#' # 1. Define data
#' my_genes <- c("Gnas1", "Tp53", "Unknown")
#' conversion_table <- c("Gnas1" = "Gnas", "Tp53" = "TP53")
#'
#' # 2. Use the standard function
#' vct_lookup(my_genes, conversion_table)
#'
#' # 3. Use the infix operator
#' my_genes %map_by% conversion_table
#'
vct_lookup <- function(query, ref) {
  # Input validation (optional but recommended)
  if (is.null(names(ref))) {
    stop("The reference vector (ref) must be a named vector.")
  }

  # Logical indexing for fast vectorization
  matches <- query %in% names(ref)
  query[matches] <- ref[query[matches]]
  
  return(query)
}

#' @rdname hash_lookup
#' @export
`%vlookup%` <- vct_lookup