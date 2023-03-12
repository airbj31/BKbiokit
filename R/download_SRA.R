#' download_SRA
#'
#' this function downloads fastq file using aws client into target folder, decrypt downloaded file into fastq files, and then remove the original file
#' aws client and fasterq-dump in sra-toolkit are required to run this command.
#'
#' please see the following links to download required packages.
#' aws client: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html
#' fasterq-dump: https://github.com/ncbi/sra-tools
#' @param x a SRA Run file accession.
#' @param dest target forder to download. if the folder is not exists, the function will make it.
#' @return downloaded file
#' @import dplyr
#' @examples
#' download_SRA("SRR17041298")
#' @export

download_SRA <-function(x, dest=".//",credential="no-sign-request",...) {
   if (!dir.exists(dest) & dest != ".//") {dir.create(dest)}
   exec<-paste0("aws s3 cp s3://sra-pub-run-odp/sra/",x,"/",x," ",dest," --",credential)
   system(exec)
   exec <-paste0("fasterq-dump --split-3 ",dest,"/",x)
   system(exec)
   exec <-paste0("rm ",dest,"/",x)
   system(exec)
}
