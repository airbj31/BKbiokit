#' xml2vector
#'
#' this function extract spacific path(xpath) node values and convert them into vectors
#'
#' @param xml xml object
#' @param path xpath
#' @return vector of values for given xpath.
#' @import xml2
#' @examples
#' xml_doc<-xml2::read_xml("<animals><pet>Cat</pet><pet>Dog</pet><insect>Ant</insect><insect>Bee</insect></animals>")
#' xml2vector(xml_doc,"//pet")
#' xml2vector(xml_doc,"//insect")
#' @export

xml2vector<-function(xml,path) {
  xml |> xml2::xml_find_all(xpath=path) |>xml2::xml_text()
}
