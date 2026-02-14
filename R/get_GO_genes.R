#' get GO_IDs from full
#' @param go_id 문자열(Character). (예: "GO:0002376").
#' @param org_db dbname. (e.g:org.Mm.eg.db - mouse db)
#' @param evidence_filter = c("EXP","IDA","IEA)
#' @export
get_GO_genes <- function(go_id, org_db = org.Mm.eg.db, evidence_filter = c("EXP", "IDA", "IEA")) {
  tryCatch({
    go_data <- AnnotationDbi::select(
      x = org_db,
      keys = go_id,
      columns = c("EVIDENCE", "SYMBOL"),
      keytype = "GOALL"
    )
    if (!is.null(evidence_filter)) {
      go_data <- go_data[go_data$EVIDENCE %in% evidence_filter, ] |> pull(SYMBOL)
    } else {
      go_data <- go_data |> pull(SYMBOL)
    }

    unique_entrez <- unique(go_data)

    return(go_data)

  }, error = function(e) {
    message("Error: ", e$message)
    return(NULL)
  })
}

#' get GO_IDs from full
#' @param go_id 문자열(Character). (예: "GO:0002376").
#' @param ontology CC or BP or MF
#' @export
get_direct_children <- function(go_id, ontology = "CC") {
  map_obj <- switch(ontology,
                    "CC" = GOCCCHILDREN,
                    "BP" = GOBPCHILDREN,
                    "MF" = GOMFCHILDREN)
  children_ids_list <- mget(go_id, map_obj, ifnotfound = NA)
  children_ids <- children_ids_list[[1]]
  if (any(is.na(children_ids))) {
    message(paste("No children found for:", go_id, "(It might be a leaf node)"))
    return(NULL)
  }
  children_names <- Term(children_ids)
  result_df <- data.frame(
    GO_ID = names(children_names),
    Term_Name = children_names,
    stringsAsFactors = FALSE
  )
  rownames(result_df) <- NULL

  return(result_df)
}

clock_gen, download_GEO, efetch, expit, get_GO_genes, get_direct_children, logit, makeGEOPath, read_GSE, read_GSM, read_geov2
