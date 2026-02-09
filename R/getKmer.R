#' getKmer
#'
#' This function is R version feature counting for FFP by vector rotating.
#'
#' @param x DNA or character vector
#' @param k k-mer
#' @return vector of values for given xpath.
#' @import xml2
#' @examples
#' xml_doc<-xml2::read_xml("<animals><pet>Cat</pet><pet>Dog</pet><insect>Ant</insect><insect>Bee</insect></animals>")
#' xml2vector(xml_doc,"//pet")
#' xml2vector(xml_doc,"//insect")
#' @export
getKmer <- function(df,k=18,count=T,norm=T) {
  df <- df
  tmp<-list();
  k<-k;
  #for(j in 1:k) { tmp[[j]] <- x[j:(length(x)-(k-j))] };
  tmp <- lapply(1:k,FUN=extVct,df,k)
  tmp <- data.frame((sapply(tmp,c))) %>% as_tibble()
  tmp <- tmp %>% unite(feature,all_of(colnames(tmp)),sep="")
  if(isTRUE(count) && isTRUE(norm)) { tmp <- tmp %>% group_by(feature) %>% summarize(n=n()/dim(tmp)[1]); }
  else if(isTRUE(count)) { tmp <- tmp %>% group_by(feature) %>% summarize(n=n()); }
  return(tmp)
}

