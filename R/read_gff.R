#' read_gff
#'
#' read gff file into wide-form tibble Object.
#'
#' @param x gff file
#' @return data frame which contained cleanned
#' @import tidyverse
#' @import stringr
#' @export

read_gff <- function(x,comment="#") {
  message("reading ",x)
  gff <- read_tsv(x,comment=comment,col_names=c("seqid","source","type","start","end","score","strand","phase","attributes"),col_types=cols())

  OUT <- list()
  var<-unique(gff$type)

  for(i in var)
  {
  message(paste0("Extracting ",i),appendLF = FALSE)
  OUT[[i]]     <- gff %>% dplyr::filter(type==i)
  OUT[[i]]$idx <- 1:dim(OUT[[i]])[1]
  lmax<-max(stringi::stri_count(OUT[[i]]$attributes,fixed=";"))
  OUT[[i]] <- OUT[[i]] %>% separate(attributes,into=paste0("attr",1:(lmax+1)),sep=";",fill="right") %>%
    gather(var,val,paste0("attr",1:(lmax+1))) %>% dplyr::filter(!is.na(val)) %>% select(-var) %>%
    separate(val,into=c("attr","value"),sep="=") %>% spread(attr,value,fill=NA) %>% select(-type,-idx)
  v2<-colnames(OUT[[i]])
  v2<-v2[9:length(v2)]

  for(j in v2)
  {
    #OUT[[i]][,j] <- gsub("%3B",";",OUT[[i]] %>% pull(j))
    OUT[[i]][,j] <- stringi::stri_replace_all_regex(OUT[[i]] %>% pull(j),
                                                    pattern     = c("%3B","%3D","%26","%2C"),
                                                    replacement = c(";","=","&",","),
                                                    vectorize   = FALSE)
  }
  message("... Done")
  }
  return(OUT)
}
