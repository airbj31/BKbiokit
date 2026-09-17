#' Standard genetic table
#'
#' A dataset containing triplet codon(codon), 1-letter amino acids (AA1), 3-letter amino acids (AA3), and property (property)
#'
#' @format A data frame with 64 rows and 3 variables
#' \describe{
#'  \item{codon}{triplet codon AAA to TTT}
#'  \item{AA1}{1 letter amino acid code}
#'  \item{AA3}{3 letter amino acid code}
#'  \item{property}{chemical property}
#'  }
"standard_genetic_code"

#' UCSC GRCm38 to GRCm39 gene symbol changes
#'
#' A dataset containing triplet codon(codon), 1-letter amino acids (AA1), 3-letter amino acids (AA3), and property (property)
#'
#' @format A data frame with 64 rows and 3 variables
#' \describe{
#'  \item{GRCm38_igenome}{UCSC/igenome gene symbol}
#'  \item{GRCm39}{mm10 gene symbol}
#'  }
"GRCm38_to_GRCm39"

#' Gene Ontology (GO) Evidence Code Metadata
#'
#' A dataset containing the Gene Ontology evidence codes, their broad categories,
#' and detailed descriptions. This metadata is useful for filtering or annotating
#' enrichment analysis results based on the reliability of the evidence.
#'
#' @format A data frame with 3 columns:
#' \describe{
#'   \item{code}{Short abbreviation of the GO evidence code (e.g., EXP, IDA, IEA).}
#'   \item{category}{The high-level classification of the evidence (e.g., Experimental,
#'   Computational Analysis, Author Statement).}
#'   \item{description}{A full definition of what the specific evidence code represents.}
#' }
#' @source \url{http://geneontology.org/docs/guide-go-evidence-codes/}
#' @usage data(GO_evidence_meta)
"go_evidence_meta"
