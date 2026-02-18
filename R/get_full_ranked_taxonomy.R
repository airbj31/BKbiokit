#' get_full_ranked_taxonomy
#'
#' download and read get_full_ranked_taxonomy
#'
#' @param locus (required)
#'

get_full_ranked_taxonomy <-function() {
  new_taxdump<-"https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip"
  download.file("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip")
  content_lists <- utils::unzip(zip = new_taxdump, list = TRUE)
}
## https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump/

## out1<-xml_children(nodes) |> xml_attr("attribute_name")
## out2<-xml_children(nodes) |> xml_text("attribute_name")
