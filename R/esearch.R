#' esearch
#'
#' search NCBI query NCBU using NCBI eutility API.
#'
#' @param query NCBI string query. for more information please read
#' @param db NCBI DBs described in details
#' @param api_key NCBI api_key. you can get your own api_key from NCBI setting page after log-in.
#' @param usehistory T/F. save result in your search history. you can also get WebEnv information. (default=T)
#' @param retstart number of skipping from the search  (Default=0)
#' @param retmax total number of unique identifiers from the retrived set to be shown in the output (default=20)
#' @param rettype (optional) "uilist" or "count" (default="uilist")
#' @param datetype (optional) one from c("crdt","edat","pdat","mhda")
#' @return vector of IDs.
#' @import tidyverse
#' @import xml2
#' @export
#' @details
#' available DBs :
#'   - https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
#'   - or use the function einfo()
#'
#' Followings are examples of search field tags from https://pubmed.ncbi.nlm.nih.gov/help/#proximity-searching.
#' You may also get the information by `einfo(db="pubmed")$Field`
#'
#' | Affiliation [ad]      | All Fields [all]          | Article Identifier [aid]  |
#' | --------------------- | ------------------------- | ------------------------- |
#' | Author [au]           | Author Identifier [auid]  | Last Author Name [lastau] |
#' | PMID [pmid]           | Title/Abstract [tiab]     | Grant Number [gr]         |
#' | Publication Date [dp] | Title [ti]                | MeSH Terms [mh]           |
#' | Text Words [tw]       | Title [ti]                | MeSH Major Topic [majr]   |

esearch <- function(query,
                    db="pubmed",
                    api_key="",
                    usehistory=TRUE,
                    retmax=20,
                    retstart=0,
                    rettype="uilist"
                    ) {
  ## error


  ## query manipulation
  query <- query |> gsub(' ','+') |> gsub('#','%23') |> gsub('"','%22')

  ## db
  if(!db %in% c("bioproject","biosample","books","cdd","gap",
                "dbvar","gene","genome","gds","geoprofiles",
                "homologene","mesh","toolkit","nlmcatalog","nuccore",
                "popset","probe","protein","proteinclusters","pcassay",
                "pccompound","pcsubstance","pubmed","pmc","snp",
                "sra","structure","taxonomy")) {stop("db isn't involved in the allowed list. please see the manual of the function.")}
  URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  URL <- paste0(URL,"db=",db)
  if(query != "") {URL <- paste0(URL,"&term=",query)}
  if(api_key != "") {URL <- paste0(URL,"&api_key=",api_key)}
  if(isTRUE(usehistory)) {URL <-paste0(URL,"&usehistory=y")}
  if(retmax != "") {URL <- paste0(URL,"&retmax=",retmax)}
  if(retstart != "") {URL <- paste0(URL,"&retstart=",api_key)}
  if(rettype != "") {URL <- paste0(URL,"&rettype=",rettype)}

  search_result<-xml2::read_xml(URL)
  total<-search_result |> xml2::xml_find_all(xpath = ".//Count") |> xml_text(trim=TRUE) |> as.numeric()
  findings <- search_result %>% xml2::xml_find_all(xpath = ".//IdList/Id") |> xml2::xml_text(trim=TRUE) |> as.numeric()
  translated_query <- search_result %>% xml2::xml_find_all(xpath = ".//QueryTranslation") |> xml2::xml_text(trim=TRUE)

  message(paste0(length(findings)," among ", total, " entries were found."))
  message(paste0("the query is translated into ",translated_query))
  message("The followings are PMIDs")
  return(findings)
}
