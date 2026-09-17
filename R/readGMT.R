#' Read/Write MySigDb gmt file format
#'
#' @name read_gmt
#' @author KIM, BYUNGJU
#' @description Read GMT files into list object
#' @param x gmt file name
#' @return list object which containing all gene sets in a given gmt file.
#' @export

read_GMT <- function(x,simple=TRUE) {
  result_list <- list()
  GMT <- readLines(x)
  for (line in GMT) {
    split_line <- unlist(strsplit(line, "\t"))
    if (length(split_line) > 0) {
      element_name <- split_line[1]
      element_vector <- split_line[-c(1,2)]
      ##element_vector <- as.numeric(element_vector)
      if(simple) {
        result_list[[element_name]]  <- element_vector
      } else {
        result_list[[element_name]]$url   <- split_line[2]
        result_list[[element_name]]$genes <- element_vector
      }
    }
  }
  return(result_list)
}

#' write GMT file from list vector
#'
#' @rdname read_gmt
#' @param x list file.
#' @param file output file
#' @param url url
#' @return saving file.
#' @export
write_gmt <- function (x, file, url = "BKbiokit")
{
  tmp<-list()
  text <- c()
  for (i in names(x)) {
    if ("url" %in% names(x[[i]])) {
      current_url <- x[[i]][["url"]]
    } else {
      current_url <- url
    }
    if (!"genes" %in% names(x[[i]])) {
      text <- c(text, paste0(paste0(i, "\t", current_url,
                                    collapse = "\t"), "\t", paste0(x[[i]],
                                                                   collapse = "\t"), collapse = "\t"))
    } else {
    text <- c(text, paste0(paste0(i, "\t", current_url,
                                  collapse = "\t"), "\t", paste0(x[[i]][["genes"]],
                                                                 collapse = "\t"), collapse = "\t"))
    }
  }
  cat(trimws(text), file = file, sep = "\n")
}
