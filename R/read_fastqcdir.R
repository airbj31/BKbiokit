#' read_fastqcdir
#'
#' this function read all fastqc output (.zip files within specific directory) and saving them into
#'
#' @param DIR a directory which stores fastqc results.
#' @return list object which contains aggregated all fastqc reports within the directory.
#' @examples
#' read_fastqcdir(DIR)
#' @import tidyverse
#' @export
read_fastqcdir <- function(DIR)
{
  vc<-list.files(DIR)
  vc<-vc[grep("_fastqc.zip",vc)]
  full_report<-list()

  pb <- progress::progress_bar$new(format = " processing [:bar] :percent eta: :eta",
                                   total=length(vc),clear=FALSE, width=60
  )

  for(i in 1:length(vc))
  {
    r1<-read_fastqc(FILE=paste0(DIR,"/",vc[i]))
    full_report[["summary"]]          <- bind_rows(full_report[["summary"]],r1[["summary"]])
    for(J in names(r1[["fastQC"]]))
    {
      full_report[["fastQC"]][[J]] <- bind_rows(full_report[["fastQC"]][[J]],r1[["fastQC"]][[J]])
    }
    pb$tick()
  }
  return(full_report)
}
