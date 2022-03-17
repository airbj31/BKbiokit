#' trimmomatic
#'
#' Trimmomatic: A flexible read trimming tool for Illumina NGS data
#'
#'
#' @param x sample fastq file basename
#' @param ext extension for fastq file. if the length of ext vector is 1, we consider it as single-end reads.
#' @param out1
#' @param out2
#' @threads
#'

trimmomatic <- function(x,
                    ext=c("_1.fastq.gz","_2.fastq.gz"),
                    out1=c("_1.trimmed.fastq.gz","_1.unpaired_trimmed.fastq.gz"),
                    out2=c("_2.trimmed.fastq.gz","_2.unpaired_trimmed.fastq.gz"),
                    threads=1,
                    logfile="",
                    summary="",
                    phred=33,
                    AdapFile="",
                    AdapAllowMaxMM=2,
                    AdapPEQ=30,
                    AdapSEQ=10,
                    leading=15,
                    trailing=15,
                    command="") {
  ## x
  ## check the dependency
  suppressWarnings(CMD <- system("which trimmomatic",intern=TRUE))
  suppressWarnings(JAVA <- system("which JAVA",intern=TRUE))
  if(length(CMD)==0) {stop("You should install trimmomatic or the path is not")}


}
