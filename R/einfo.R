#' einfo
#'
#'
#' Provides the number of records indexed in each field of a given database, the date of the last update of the database, and the available links from the database to other Entrez databases.
#' NCBI E-utility wrapper series for R users.
#'
#' @param db database name. if this field is "", it shows all of the available databases.
#' @import xml2
#' @export
#' @details
#' https://www.ncbi.nlm.nih.gov/books/NBK25497/
einfo <- function(db="") {
  URL<-"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"

  if(db=="") {
    search_result<-xml2::read_xml(URL)
    return(search_result |> xml2::xml_find_all(xpath = "//DbList/DbName") |> xml2::xml_text())
  }  else {
    URL<-paste0(URL,"?db=",db)
    search_result<-xml2::read_xml(URL)
    Output <- list()
    Output$DbName      <- search_result |> xml2::xml_find_first(xpath=".//DbName")      |> xml2::xml_text()
    Output$MenuName    <- search_result |> xml2::xml_find_first(xpath=".//MenuName")    |> xml2::xml_text()
    Output$Description <- search_result |> xml2::xml_find_first(xpath=".//Description") |> xml2::xml_text()
    Output$DbBuilds    <- search_result |> xml2::xml_find_first(xpath=".//DbBuild")     |> xml2::xml_text()
    Output$Count       <- search_result |> xml2::xml_find_first(xpath=".//Count")       |> xml2::xml_text()
    Output$LastUpdate  <- search_result |> xml2::xml_find_first(xpath=".//Count")       |> xml2::xml_text()

    FieldColumns<-xml2::xml_name(xml2::xml_children(search_result |> xml2::xml_find_all(xpath=".//FieldList/Field"))) |> unique()
    LinkColumns<-xml2::xml_name(xml2::xml_children(search_result |> xml2::xml_find_all(xpath=".//LinkList/Link"))) |> unique()

    Output$Link  <- do.call(cbind,lapply(setNames(paste0(".//LinkList/Link/",LinkColumns),LinkColumns),function(x) BKbiokit::xml2vector(search_result,x))) |> as.data.frame() |> tibble::as_tibble()
    Output$Field <- do.call(cbind,lapply(setNames(paste0(".//FieldList/Field/",FieldColumns),FieldColumns),function(x) BKbiokit::xml2vector(search_result,x))) |> as.data.frame() |> tibble::as_tibble()

    Output$Link  <- Output$Link %>% mutate_if(is.character, list(~na_if(.,"")))
    Output$Field <- Output$Field %>% mutate_if(is.character, list(~na_if(.,""))) %>% mutate(IsDate = (IsDate=="Y"),IsNumerical=IsNumerical=="Y",SingleToken=SingleToken=="Y",Hierarchy=Hierarchy=="Y",IsHidden=IsHidden=="Y")

        colnames(Output$Field)<-FieldColumns
    colnames(Output$Link)<-LinkColumns

    return(Output)
  }
}
