#' BamDiversity
#' UNDER DEV
#' K-mer diversity of mapped reads.
#' returns data frame containing 3 columns
#'  1. idx     : alignment starting position of the reference sequence
#'  2. Sequence: Nucleotide sequence including INDEL.
#'  3. count   : # of sequence
#'
#' @param bam  (required) bam file to read
#' @param kmer (required) kmer length
#' @param chrom (required) sequence ids for test.
#' @param clip including clip or not
#' @return list or vector object.
#' @import tidyverse
#' @import Rsamtools
#' @import stringi
#'
#'
BamDiversity <- function(param,kmer,chrom) {
  BamRecords <- Rsamtools::BamFile(bam,yieldSize = 1000000)
  open(BamRecords)
  result<-tibble()
  while(length(chunk <- scanBam(BamRecords)[[1]])) {
    chunk$qual <- as.character(chunk$qual)
    chunk$seq <- as.character(chunk$seq)
    chunk   <- as_tibble(chunk)
    chunk   <- chunk |> mutate(cigar_ops=str_extract_all(cigar,"[MIDNSHP=X]"),cigar_len=str_extract_all(cigar,"\\d+"))

    ## divide reads by the cases
    sftclp  <- chunk |> dplyr::filter(grepl("S",cigar))
    sftclp_ID  <- sftclp |> dplyr::filter(grepl("I",cigar)) |> dplyr::filter(grepl("D",cigar))
    sftclp_I   <- sftclp |> dplyr::filter(grepl("I",cigar)) |> dplyr::filter(!grepl("D",cigar))
    sftclp_D   <- sftclp |> dplyr::filter(grepl("D",cigar)) |> dplyr::filter(!grepl("I",cigar))

    if(nrow(sftclp_ID)==0) { rm(sftclp_ID)}
    if(nrow(sftclp_I)==0) { rm(sftclp_I)}
    if(nrow(sftclp_D)==0) { rm(sftclp_D)}
    sftclp <- sftclp |> dplyr::filter(!grepl("I",cigar),!grepl("D",cigar))

    ## divide reads by the cases
    chunk     <- chunk |> dplyr::filter(!grepl("S",cigar))
    chunk_ID  <- chunk |> dplyr::filter(grepl("I",cigar)) |> dplyr::filter(grepl("D",cigar))
    chunk_I   <- chunk |> dplyr::filter(grepl("I",cigar)) |> dplyr::filter(!grepl("D",cigar))
    chunk_D   <- chunk |> dplyr::filter(grepl("D",cigar)) |> dplyr::filter(!grepl("I",cigar))
    chunk <- chunk |> dplyr::filter(!grepl("I",cigar),!grepl("D",cigar))

    sftclp  <- sftclp |> mutate(clip.direction=ifelse(grepl("^[0-9]+S",cigar),"Forward",ifelse(grepl("[0-9]+S$",cigar),"Backward",NA)),clip.length=as.numeric(gsub("S","",str_extract(cigar,"[0-9)]+S"))))
    sftclp  <- sftclp |> mutate(clipped.seq=case_when(clip.direction=="Forward" ~ substr(seq,1,clip.length),
                                                  clip.direction=="Backward" ~ substr(seq,stringi::stri_length(sftclp$seq)-clip.length+1,stringi::stri_length(sftclp$seq))))
    sftclp  <- sftclp |> mutate(remained.seq=case_when(clip.direction=="Forward" ~ substr(seq,clip.length+1,stringi::stri_length(seq)),
                                                       clip.direction=="Backward" ~ substr(seq,1,stringi::stri_length(seq)-clip.length)))

    ## Now chunk do not has soft-clipping



  }
  close(BamRecords)
  return(result)
}

addDeletion  <- function() {}


addInsertion <- function() {}
