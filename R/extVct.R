#' extVct
#'
#' return character vector ranging index x to k-x for given vector. this function is developed to k-mer extraction for a given vector.
#'
#' @param x index. numeric value
#' @param y chracter vector to manipulate.
#'
#' @return vector of characters. `y[x:length(y)-(k-x)]`
#' @examples
#' testvector<-sample(letters,1000,replace=T)
#' extVct(1,testvector,10)[1:10]
#' extVct(2,testvector,10)[1:10]
#' @export
extVct <- function (x,y,k) {return(y[x:(length(y)-(k-x))])}
