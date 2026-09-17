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

    # 1. 맵 객체 설정 (Children 대신 Offspring을 원하신다면 OFFSPRING으로 변경 가능)
    # 패키지 내부에 이 객체들이 로드되어 있어야 합니다.
    map_obj <- switch(ontology,
                      CC = GO.db::GOCCCHILDREN,
                      BP = GO.db::GOBPCHILDREN,
                      MF = GO.db::GOMFCHILDREN)

    # 2. mget 대신 AnnotationDbi::exists와 get을 조합하여 안전하게 추출
    # mget은 환경 객체 조건이 까다롭지만, get은 Bimap 객체도 잘 처리합니다.
    if (!AnnotationDbi::exists(go_id, map_obj)) {
      message(paste("No children found for:", go_id, "(It might be a leaf node or invalid ID)"))
      return(NULL)
    }

    children_ids <- AnnotationDbi::get(go_id, map_obj)

    # NA 처리
    children_ids <- children_ids[!is.na(children_ids)]

    if (length(children_ids) == 0) {
      return(NULL)
    }

    # 3. Term 정보 가져오기
    children_names <- AnnotationDbi::Term(children_ids)

    # 4. 데이터프레임 생성
    result_df <- data.frame(
      GO_ID = names(children_names),
      Term_Name = as.character(children_names),
      stringsAsFactors = FALSE
    )

    rownames(result_df) <- NULL
    return(result_df)
  }

library(AnnotationDbi)
library(org.Mm.eg.db)
library(GO.db)
library(dplyr)

#' Extract Genes for a GO Term with Evidence-Based Recursive Mapping
#'
#' This function retrieves a list of gene symbols associated with a specific GO ID
#' and all its descendant (offspring) terms. Unlike the standard `GOALL` mapping,
#' this function allows for strict filtering based on Evidence Codes before
#' aggregating genes across the GO hierarchy.
#'
#' @param go_id Character. The parent Gene Ontology ID (e.g., "GO:0002376").
#' @param org_db An AnnotationDb object. The organism-specific database
#'   (default: `org.Mm.eg.db`).
#' @param evidence_filter Character vector. A list of GO evidence codes to include
#'   (e.g., `c("EXP", "IDA")` for experimental levels, or `c("EXP", "IDA", "IEA")`
#'   to include electronic annotations). If NULL, all evidence codes are included.
#'
#' @return A character vector of unique gene symbols. Returns `NULL` if no genes
#'   are found or if an error occurs.
#'
#' @details
#' The standard `keytype = "GOALL"` approach relies on a pre-computed table that
#' merges all evidence types, which can lead to "data dilution" by low-confidence
#' annotations (like IEA) across all levels.
#'
#' This function ensures data integrity by:
#' 1. Identifying all offspring nodes using the `GO.db` hierarchy.
#' 2. Querying each node individually using `keytype = "GO"` (Direct Annotation).
#' 3. Filtering by the user-defined evidence codes at each step.
#'
#' This is particularly useful for creating multi-omics feature sets (GMT files)
#' categorized by different biological confidence levels.
#'
#' @examples
#' \dontrun{
#' # Level 1: Experimental evidence only
#' genes_lvl1 <- get_filtered_recursive_genes("GO:0002376",
#'                                             evidence_filter = c("EXP", "IDA"))
#'
#' # Level 4: Including electronic annotations (IEA)
#' genes_lvl4 <- get_filtered_recursive_genes("GO:0002376",
#'                                             evidence_filter = c("EXP", "IDA", "IEA"))
#' }
#' @export
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
