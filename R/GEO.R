#' makeGEOPath
#'
#' @author Byungju Kim (byungju@bu.edu)
#' @param GSE GSE accession
#' @param type type of file format you download
makeGEOPath <- function(GSE,type=c("matrix","suppl","xml",""),prefix="GEO",...) {
  DIR  <- paste0("GSE",substr(GSE,4,stringi::stri_length(GSE)-3),"nnn")
  FULLPATH <- paste0(prefix,"/",DIR,"/",GSE,"/",type,"/")
  return(FULLPATH)
}


#' read_geo
#'
#' `logit` returns inverse of logistic transformation
#'
#' @author Byungju Kim (byungju@bu.edu)
#' @param GSE GSE accession
#' @param prefix the directory which you saved all of GSE structure.
#' @return inverse of logistic transformation of a given data x
read_geov2 <- function(GSE,prefix="GEO",...) {
  DIR  <- paste0("GSE",substr(GSE,4,stringi::stri_length(GSE)-3),"nnn")
  FULLPATH <- paste0(prefix,"/",DIR,"/",GSE,"/matrix/")
  FILES <- list.files(FULLPATH)
  FILES <- FILES[grepl(".chrs.v2$",FILES)]
  output<-tibble()
  for(i in FILES) {
    message(paste("Reading ",i))
    tmp1  <- read_tsv(paste0(FULLPATH,i),...)
    tmp2  <- read_tsv(paste0(FULLPATH,gsub(".chrs.v2",".chrs.uniq",i)),...)
    tmp2  <- tmp2 |> dplyr::filter(!is.na(column))
    CNAME <- tmp2$column
    tmp2  <- tmp2$value |> t() |> as.data.frame() |> as_tibble()
    colnames(tmp2)<-CNAME
    tmp <- left_join(tmp1,tmp2,by="GSE")
    output <- bind_rows(output,tmp)
  }
  CNAME <- colnames(output)
  remNAME <- CNAME[unique(c(grep("protocol",CNAME),grep("contact",CNAME),grep("count",CNAME)))]
  output <- output |> select(-all_of(remNAME))
  if("age" %in% colnames(output)) { output <- output |> mutate(age=as.character(age))}
  return(output)
}

#` download_GEO
#'
#' download GEO matrix from GEO ftp.
#'
#' @author Byungju Kim (byungju@bu.edu)
#' @param GSE list of GSEs for downloading
#' @param prefix base directory
#' @param type file type to download
#' @examples
#' gexpit(runif(100,0,1))
download_GEO <- function(GSE,prefix="GEO",...) {
 ## GSE


}
