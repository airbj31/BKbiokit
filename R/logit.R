#' logit transformation
#'
#' `logit` returns inverse of logistic transformation
#'
#' @author Byungju Kim (byungju@bu.edu)
#' @param x data
#' @return inverse of logistic transformation of a given data x
#' @examples
#' gexpit(runif(100,0,1))
#' @export
#'
logit <- function(x) {
  log(x/(1 - x))
}

#' expit transformation
#'
#' `expit` returns logistic transformation of a given data x
#'
#' @author Byungju Kim (byungju@bu.edu)
#' @param x data
#' @return logistic transformation of x
#' @examples
#' gexpit(runif(100,0,1))
#' @export
expit <- function(x) {
  return(exp(x)/(1 + exp(x)))
}
