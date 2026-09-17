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

library(AnnotationDbi)
library(org.Mm.eg.db)
library(GO.db)
library(dplyr)

get_gogenes_recursive <- function(go_id, org_db = org.Mm.eg.db, evidence_filter = c("EXP", "IDA")) {
  # 1. Ontology 종류 판별 (BP, CC, MF)
  ont_type <- tryCatch(AnnotationDbi::Ontology(go_id), error = function(e) return(NULL))
  if (is.null(ont_type)) return(NULL)
  # 2. 해당 Ontology에 맞는 Offspring DB 선택
  offspring_map <- switch(ont_type,
                          "BP" = GOBPOFFSPRING,
                          "CC" = GOCCOFFSPRING,
                          "MF" = GOMFOFFSPRING)
  # 3. 하위 ID 추출 (자신 포함)
  # mget은 ID가 없을 경우 에러 대신 list(NULL)을 반환하여 안전함
  offspring_ids <- unlist(mget(go_id, offspring_map, ifnotfound = NA))
  all_ids <- unique(c(go_id, offspring_ids))
  all_ids <- all_ids[!is.na(all_ids)]
  # 4. 유전자 추출 (Direct Mapping 'GO' 사용)
  # keys가 유효한지 사전에 검사하여 에러 방지
  valid_keys <- intersect(all_ids, keys(org_db, keytype = "GO"))
  if (length(valid_keys) == 0) {
    return(NULL)
  }
  tryCatch({
    go_data <- AnnotationDbi::select(x = org_db,
                                     keys = valid_keys,
                                     keytype = "GO",
                                     columns = c("EVIDENCE", "SYMBOL"))
    # 5. Evidence 필터링
    if (!is.null(evidence_filter)) {
      filtered_symbols <- go_data %>%
        filter(EVIDENCE %in% evidence_filter) %>%
        pull(SYMBOL) %>%
        unique() %>%
        na.omit()
      return(as.character(filtered_symbols))
    } else {
      return(unique(na.omit(go_data$SYMBOL)))
    }
  }, error = function(e) {
    message("Error in select for ", go_id, ": ", e$message)
    return(NULL)
  })
}
