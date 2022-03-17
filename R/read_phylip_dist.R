#' read_phylip_dist: read Phylip style distance matrix as distance matrix in R
#'
#' @name read.phylip_dist
#' @author Byungju Kim (bjk@inu.ac.kr)
#' @param file a file name,
#' @param rm.list list for sample exclusion.
#' @param skip a numeric valueof line skipping.
#' @export
#' @examples
#' read.dist(matrix.tsv)

read_phylip_dist <- function(file,rm.list=c(),skip=1)
{
  nskip=skip
  ncol<-as.numeric(readLines(file,n=1))
  tmat<-read_table2(file,col_names=paste0("V",seq(1,ncol+1,by=1)),skip=nskip)
  CNAME<-tmat$V1
  tmat<-tmat %>% select(-V1) %>% as.matrix()
  colnames(tmat)<-CNAME
  row.names(tmat)<-CNAME
  tmat[upper.tri(tmat)] <- t(tmat)[upper.tri(tmat)]

  if(length(rm.list)==0) {return(as.dist(tmat))}
  else {return(as.dist(tmat[!row.names(tmat) %in% rm.list,!colnames(tmat) %in% rm.list]))}
}
