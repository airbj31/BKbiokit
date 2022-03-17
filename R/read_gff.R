#' read_gff
#'
#' read gff file into wide-form tibble Object.
#'
#' @param x gff file
#' @return data frame which contained cleanned
#' @import tidyverse
#' @import stringr
#'
#'


## x<-"../../Research/bee-genome/refseq/genome_assemblies_genome_gff/ncbi-genomes-2022-03-14/GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz"
##x <-"../../Research/bee-genome/refseq/
read_gff <- function(x,comment="#") {
  message("reading ",x)
  gff <- read_tsv(x,comment=comment,col_names=c("seqid","source","type","start","end","score","strand","phase","attributes"),col_types=cols())
  ##message("reading ... done")
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
  message("... Done")
  }


  return(OUT)
}


#
# library(tidyverse)
# selected_region <-
#   c(
#     "NC_001566.1",
#     ## MT
#     "NC_037638.1",
#     ## LG01
#     "NC_037639.1",
#     ## LG02
#     "NC_037640.1",
#     ## LG03
#     "NC_037641.1",
#     ## LG04
#     "NC_037642.1",
#     ## LG05
#     "NC_037643.1",
#     ## LG06
#     "NC_037644.1",
#     ## LG07
#     "NC_037645.1",
#     ## LG08
#     "NC_037646.1",
#     ## LG09
#     "NC_037647.1",
#     ## LG10
#     "NC_037648.1",
#     ## LG11
#     "NC_037649.1",
#     ## LG12
#     "NC_037650.1",
#     ## LG13
#     "NC_037651.1",
#     ## LG14
#     "NC_037652.1",
#     ## LG15
#     "NC_037653.1"  ## LG16
#   )
#
# ## showing Mitochondria
#
# gff$gene %>% dplyr::filter(seqid == "NC_001566.1") %>%
#   ggplot(aes(5, start)) +
#   annotate(
#     "rect",
#     xmin = 4,
#     xmax = 6,
#     ymin = 1,
#     ymax = 16323,
#     alpha = 0.2,
#     color = "black"
#   ) +
#   geom_segment(data = tibble(x = rep(6, 33), y = seq(1, 16323, by = 500)),
#                aes(
#                  x = x,
#                  xend = x + 0.5,
#                  y = y,
#                  yend = y
#                ),
#                color = "black") +
#   geom_segment(data = tibble(x = rep(6, 164), y = seq(1, 16323, by = 100)),
#                aes(
#                  x = x,
#                  xend = x + 0.2,
#                  y = y,
#                  yend = y
#                ),
#                color = "black") +
#   geom_segment(
#     data = gff$gene %>% dplyr::filter(seqid == "NC_001566.1"),
#     aes(
#       x = 5,
#       xend = 5,
#       y = 1,
#       yend = 16343
#     ),
#     size = 2,
#     alpha = 0.3
#   ) +
#   geom_rect(
#     data = gff$gene %>% dplyr::filter(seqid == "NC_001566.1", strand == "+"),
#     aes(
#       xmin = 5,
#       xmax = 6,
#       ymin = start,
#       ymax = end,
#       fill = gene
#     ),
#     size = 1
#   ) +
#   #geom_gene_arrow(data=gff$gene %>% dplyr::filter(seqid =="NC_001566.1",strand=="+"),
#   #             aes(xmin=5.5,xmax=5.5,y=start,fill=gene)) +
#
#   geom_rect(
#     data = gff$gene %>% dplyr::filter(seqid == "NC_001566.1", strand == "-"),
#     aes(
#       xmin = 4,
#       xmax = 5,
#       ymin = start,
#       ymax = end,
#       fill = gene
#     ),
#     size = 1
#   ) +
#
#   coord_polar(theta = "y",
#               start = 0,
#               clip = "on") +
#   xlim(c(-2, 6.5)) + ylim(0, 16500)
