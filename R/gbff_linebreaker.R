#' gbff_linebreaker()
#'
#' Redirect texts into gbff multiline outputs in metadata section
#'
#' param
#'
gbff_linebreaker <- function(vc,sep=" ",linelimit=67,header="DEFINITION") {
  vc     <- as.vector(stringr::str_split_fixed(vc," ",n=Inf))
  len_vc <- stringi::stri_length(vc) +1

  LINES<-c() ## output.
  i<-1
  while(length(vc)>0) {
    point<-sum(cumsum(len_vc) %/% linelimit == 0) ## number of position of each line.
    LINES[i]<-paste0(vc[1:point],collapse=" ")
    if(sum(len_vc) %/% linelimit > 0) { vc<-vc[(point+1):length(vc)];len_vc <- stringi::stri_length(vc) +1}
    else {vc<-c()}
    i<-i+1
  }

  for(i in 1:length(LINES)) {
    LINES[i] <- ifelse(i==1,paste0("DEFINITION  ",LINES[i],"\n",collapse=""),paste0("            ",LINES[i],"\n",collapse=""))
  }
  return(LINES)
}
