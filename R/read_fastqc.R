#' read_fastqc
#'
#' this function read a given fastqc output (.zip files within specific directory) and saving them into
#'
#' @param fastqc a zip file for fastqc output (FILE_fastqc.zip)
#' @return list object which contains summary and reports for denoted qc reports.
#' @import dplyr
#' @import tidyverse
#' @import utils
#' @import tibble
#' @import readr
#' @examples
#' read_fastqc("demo_fastqc.zip")
#' @export
read_fastqc<-function(fastqc="K34_1_fastqc.zip")
{
  ## getting the file names within zip files and do not unarchive it.
  content_lists<-utils::unzip(zip = fastqc,list=TRUE)
  content_lists<-content_lists$Name[grep(".txt",content_lists$Name)]

  ## read_lines
  ## summary_tbl  <- read_tsv(unz(FILE,file=content_lists[1]))      ## Not needed.
  FILE2 <- unz(fastqc,file=content_lists[2])
  data_LINES   <- readLines(FILE2)
  on.exit(close(FILE2))

  ## every two lines are the separators for each dataset.
  ranges       <- grep(">>",data_LINES)
  fastQC_OBJ     <- list()
  fastQC_Summary <- tibble()

  for(i in seq(1,length(ranges),by=2)) {
    metric<-gsub(">>","",data_LINES[ranges[i]]) %>% stringr::str_split("\t") %>% unlist()
    header <- metric[1]
    result <- metric[2]
    fastQC_Summary<-bind_rows(fastQC_Summary,tibble(file=fastqc,metric=header,result=result))

    if(ranges[i]+1 == ranges[i+1]) {fastQC_OBJ[[header]] <-tibble() }
    else if(header=="Sequence Duplication Levels") {
      metric<-gsub(">>","",data_LINES[ranges[i]+1]) %>% stringr::str_split("\t") %>% unlist()
      header <- metric[1]
      result <- metric[2]
      fastQC_Summary<-bind_rows(fastQC_Summary,tibble(file=fastqc,metric=header,result=result))
      fastQC_OBJ[[header]] <- data_LINES[(ranges[i]+2):(ranges[i+1]-1)] %>% paste0(concat="\n") %>% read_tsv(show_col_types=FALSE) %>% mutate(file=fastqc)
    }
    else {fastQC_OBJ[[header]] <- data_LINES[(ranges[i]+1):(ranges[i+1]-1)] %>% paste0(concat="\n") %>% read_tsv(show_col_types=FALSE) %>% mutate(file=fastqc)}
  }
  out<-list(fastQC_OBJ,fastQC_Summary)
  names(out)<-c("fastQC","summary")
  return(out)
}
