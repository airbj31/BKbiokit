#' @title get full PATH of a given GSE record.
#' @author KIM, BYUNGJU
#' @description parse GSE accession and build GSE directory.
#' @param GSE single GSE accession
#' @param prefix base directory which have GSEs.
#' @return single vector of PATH

getGSEdir <- function(GSE,prefix="data/GEO") {
  res <- tryCatch({
    if(!grepl("^GSE\\d+$", GSE)) {stop("The given GSE is not GSE accession");return(TRUE)}
  },error=function(e) {message("\033[31m","error: ",e$message,"\033[0m")})
  DIR  <- paste0("GSE",substr(GSE,4,stringi::stri_length(GSE)-3),"nnn")
  FULLPATH <- paste0(prefix,"/",DIR,"/",GSE,"/matrix/")
  tryCatch({
    if(dir.exists(file.path(FULLPATH))) {return(FULLPATH)} else {stop("The given GSE is not accessible")}
  },error=function(e) {message("\033[31m","error: ",e$message,"\033[0m")})

}


#' read_GSE
#'
#' read GSE documents downloaded by BJK
#' @param GSE GSE accession
#' @param prefix base directory
#' @param keep vector of categories you want to keep in the returned object.
#' @return list object containing description about the GSE accession

read_GSE <- function(GSE,prefix="GEO",keep=c("Series_title","Series_summary","Series_overall_design","Series_type"),...) {
  FULLPATH <- getGSEdir(GSE,prefix=prefix)
  FILES <- list.files(FULLPATH)
  FILES <- FILES[grepl(".gse.xz$",FILES)]
  output<-tibble()
  i <- FILES[1]
  message(paste("Reading ",i))
  tmp1  <- read_tsv(paste0(FULLPATH,i),col_names=c("column","description"),show_col_types = FALSE)
  tmp1  <- tmp1[tmp1$column %in% keep,]
  tmp1 <- tmp1 |> dplyr::filter(!is.na(description)) |> dplyr::filter(description != "This SuperSeries is composed of the SubSeries listed below.")
  tmp1 <- tmp1 |> group_by(column) |> mutate(description=paste0(description,collapse="; ")) |> unique()
  CN <- tmp1$column
  tmp1 <- t(tmp1)
  colnames(tmp1) <-CN
  return(as.list(tmp1[-1,]))
}

#' @title read_GSM
#' @author KIM, BYUNGJU
#' @param GSE GSE accession
#' @param prefix base directory. We read data from prefix/GSEXXXnnn/GSEXXXXXX

read_GSM <- function(GSE,prefix="GEO",...) {
  FULLPATH <- getGSEdir(GSE,prefix=prefix)
  FILES <- list.files(FULLPATH)
  FILES <- FILES[grepl(".chrs.v2$",FILES)]
  output<-tibble()
  for(i in FILES) {
    message(paste("Reading ",i))
    tmp1  <- read_tsv(paste0(FULLPATH,i),...)
    tmp2  <- read_tsv(paste0(FULLPATH,gsub(".chrs.v2",".chrs.uniq",i)),...)
    tmp2  <- tmp2 |> dplyr::filter(!is.na(column))
    tmp2  <- tmp2 |> dplyr::filter(!is.na(value))
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
