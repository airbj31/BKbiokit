#' construct BWA table for specific k-mer
#'
#' Redirect texts into gbff multiline outputs in metadata section
#'
#' @param x k-mer
#' @param format file format
#'
make_bwt_tbl <- function(x) {
  RESULT<-c("A","C","G","T")
    for(i in 1:x) {
      RESULT<-as.vector(outer(RESULT,c("A","C","G","T"),paste,sep=""))
    }
    ## reducing vector size by reverse complement
    RESULT_rc<-reverseComplement(DNAStringSet(RESULT)) |> as.character()
    RESULT <- tibble(RESULT,RESULT_rc) |> mutate(new=ifelse(RESULT<RESULT_rc,RESULT,RESULT_rc)) |> pull(new) |> unique()
    rm(RESULT_rc)

    hash_tbl<-tibble()
    while(length(RESULT) > 0) {
      BWToutput <- BWTcode(RESULT[1],lexi_one = FALSE,rc=FALSE) |> unique()
      BWTinput <- rep(RESULT[1],length(BWToutput))
      tmp_df <- tibble(Hash=BWToutput,value=BWTinput)
      hash_tbl<-bind_rows(hash_tbl,tmp_df)
      RESULT <- RESULT[! RESULT %in% BWToutput]
    }
    return(hash_tbl)
  }

make_bwt_tbl2 <- function(x) {
  RESULT<-c("A","C","G","T")
  for(i in 1:x) {
    RESULT<-as.vector(outer(RESULT,c("A","C","G","T"),paste,sep=""))
  }
  hash_tbl<-tibble()
  while(length(RESULT) > 0) {
    BWToutput <- BWTcode(RESULT[1],lexi_one = FALSE,rc=TRUE) |> unique()
    BWTinput <- rep(RESULT[1],length(BWToutput))
    tmp_df <- tibble(Hash=BWToutput,value=BWTinput)
    hash_tbl<-bind_rows(hash_tbl,tmp_df)
    RESULT <- RESULT[! RESULT %in% BWToutput]
  }
  return(hash_tbl)
}


## check and remove
## appending all vector

## make_bwt_tbl <- function(x) {
##  RESULT<-c("A","C","G","T") {
##    for(i in 1:x) {
##      RESULT<-as.vector(outer(RESULT,c("A","C","G","T"),paste,sep=""))
##    }
##    hash_tbl<-tibble()
##    while(length(RESULT) > 0) {
##      BWToutput <- BWTcode(RESULT[1],lexi_one = FALSE,rc=TRUE) |> unique()
##      BWTinput <- rep(RESULT[1],length(BWToutput))
##      tmp_df <- tibble(Hash=BWToutput,value=BWTinput)
##      hash_tbl<-bind_rows(hash_tbl,tmp_df)
##      RESULT <- RESULT[! RESULT %in% BWToutput]
##    }
##    return(hash_tbl)
##  }
