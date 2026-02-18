#' gffcmpr
#'
#' The function compares two gff list loaded from read_gff function in BKbiokit
#' This function is desgined to compare two or more `liftover` or `annotation result`
#'  @param x list object loaded from read_gff
#'  @param y list object loaded from read_gff
#'  @return merged list output with some notation.
#'  @examples



gffcmpr <-function(x,y) {

  tmp <- bind_rows(x,y) %>% arrange(seqid,start,end)
  tmp1 <- tmp[1:(dim(tmp1)[1]-1),]
  colnames(tmp1)<-paste0(colnames(tmp1),".x")
  tmp2 <- tmp2[2:dim(tmp2)[1],]
  colnames(tmp2)<-paste0(colnames(tmp2),".y")

  tmp<-cbind(tmp1,tmp2)  %>% mutate(idx=1:(dim(tmp)[1]-1)) %>% select(idx,all_of(paste0(rep(colnames(gene),each=2),c(".x",".y")))) %>% dplyr::filter(source.x!=source.y,strand.x==strand.y)

  tmp %>%
    mutate(desc=case_when(
      start.x == end.y && start.y==end.y ~ "same",
      start.x =="+" && end.x < start.y  ~ "diff1",
      start.x =="-" && end.y < start.x  ~ "diff2",
      TRUE ~ "overlap"
    )) %>% select(seqid.x,idx,start.x,start.y,end.x,end.y,desc,strand.x,strand.y,source.x,source.y)


}
