#' BWTcode
#'
#' This function is introduced from Burrow-wheeler Transform to count repetitive motifs correctly.
#' the development aim
#'
#' please use vBWTcode
#' @param str a string or vector.
#' @param rc if you want to count repetitive motifs from the
#' @param lexi_one print first entry of lexicographical sorting. if FALSE, it prints every entry.
#' @return downloaded file
#' @import dplyr
#' @examples
#' BWTcode("TCACA")
#' BWTcode("TCACA",rc=TRUE)
#' BWTcode("TCACA",rc=TRUE,lexi_one=FALSE)
#' @export
#'
BWTcode <-function(str="ACACT",rc=TRUE,lexi_one=TRUE) {

  #  str<-paste0("^",str,"|")
  FINALVC<-str
  for(i in 1:(stringi::stri_length(str)-1))
  {
    #FINALVC<-c(FINALVC,paste0(substr(str,i+1,stringi::stri_length(str)),substr(str,1,i)))
    FINALVC<-c(FINALVC,paste0(substr(str,stringi::stri_length(str)-i+1,stringi::stri_length(str)),substr(str,1,stringi::stri_length(str)-i)))
  }

  if(isTRUE(rc))
  {
  str<-Biostrings::reverseComplement(Biostrings::DNAString(str)) %>% as.character()
  FINALVC<-c(FINALVC,str)
  for(i in 1:(stringi::stri_length(str)-1))
  {
    #FINALVC<-c(FINALVC,paste0(substr(str,i+1,stringi::stri_length(str)),substr(str,1,i)))
    FINALVC<-c(FINALVC,paste0(substr(str,stringi::stri_length(str)-i+1,stringi::stri_length(str)),substr(str,1,stringi::stri_length(str)-i)))
  }

  }
  if(isTRUE(lexi_one)) {return(sort(FINALVC)[1])}
  else {return(sort(FINALVC))}

}

vBWTcode<-Vectorize(BWTcode,vectorize.args="str")
