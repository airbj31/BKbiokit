#' 식별자(DOI, PMID, PMCID)로 BibTeX 가져오기
#'
#' DOI, PubMed ID(PMID), PubMed Central ID(PMCID) 중 하나를 받아
#' DOI로 변환한 뒤, doi.org content negotiation으로 BibTeX 항목을 가져옵니다.
#' 저자는 전부 그대로 두고, 필드마다 한 줄씩 보기 좋게 정리합니다.
#'
#' @param id `character(1)`. 논문 식별자. 다음 형식을 지원합니다.
#'   * DOI: `"10.1038/s41586-020-2649-2"` 또는 `"https://doi.org/10.1038/..."`
#'   * PMID: `"32939066"`
#'   * PMCID: `"PMC7759461"` (대소문자 무관)
#' @param key_type `character(1)`. BibTeX citation key 형식.
#'   * `"author_year"` (기본값): 제1저자 성 + 연도, 예: `harris2020`
#'   * `"id"`: Crossref가 준 key를 그대로 사용, 예: `Harris_2020`
#'   * 제1저자가 단체면 단체 이름에서 key를 만듭니다,
#'     예: The 1000 Genomes Project Consortium -> `1000genomes2015`
#' @param consortium `character(1)`. 제1저자가 컨소시엄 같은 단체인 논문의 저자 처리.
#'   * `"name_only"` (기본값): 단체 이름 하나만 저자로 씁니다,
#'     예: `author = {{The 1000 Genomes Project Consortium}}`.
#'     이런 논문은 보통 이렇게 인용합니다.
#'   * `"members"`: 단체 이름 뒤에 Crossref에 등록된 소속 연구자를 모두 씁니다.
#'     "Steering committee", 기관명 같은 그룹 제목은 빼고, 중복된 이름은 한 번만 씁니다.
#'
#'   제1저자가 사람인 논문에서는 중간에 나오는 컨소시엄 이름을 그대로 둡니다
#'   (예: gnomAD 논문의 "Genome Aggregation Database Consortium").
#'
#' @return `bibtex_str` 클래스가 붙은 `character(1)`. BibTeX 항목 문자열이며,
#'   콘솔에서는 `print.bibtex_str()`로 그대로 출력됩니다.
#'   파일로 저장하려면 `save_bibtex()`를 쓰세요.
#'
#' @details
#' * PMCID는 PMC ID Converter API로 DOI를 찾습니다.
#' * PMID는 PubMed E-utilities(`esummary`)로 DOI를 찾으므로
#'   PMC에 수록되지 않은 논문도 처리됩니다.
#' * DOI가 없는 논문(오래된 논문, 일부 학회 초록 등)은 에러가 납니다.
#' * 인터넷 연결이 필요하며, `url(headers = )`를 쓰므로 R 3.6.0 이상이 필요합니다.
#'
#' @seealso `save_bibtex()`,
#' [PMC ID Converter API](https://pmc.ncbi.nlm.nih.gov/tools/id-converter-api),
#' [Crossref content negotiation](https://citation.crosscite.org/docs.html)
#'
#' @examples
#' \dontrun{
#' get_bibtex("10.1038/s41586-020-2649-2")
#' get_bibtex("PMC7759461")
#' get_bibtex("32939066", key_type = "id")
#'
#' # 컨소시엄 논문
#' get_bibtex("10.1038/nature15393")                         # 단체 이름만
#' get_bibtex("10.1038/nature15393", consortium = "members") # 소속 연구자까지
#' }
#'
#' @importFrom jsonlite fromJSON
#' @export
get_bibtex <- function(id, key_type = c("author_year", "id"),
                       consortium = c("name_only", "members")) {
  key_type   <- match.arg(key_type)
  consortium <- match.arg(consortium)
  clean_id <- trimws(as.character(id))
  clean_id <- gsub("^https?://(dx\\.)?doi\\.org/", "", clean_id)
  clean_id <- sub("^doi:\\s*", "", clean_id, ignore.case = TRUE)

  if (grepl("^10\\.", clean_id)) {
    type <- "doi";   doi <- clean_id
  } else if (grepl("^PMC\\d+$", clean_id, ignore.case = TRUE)) {
    type <- "pmcid"; doi <- pmcid_to_doi(toupper(clean_id))
  } else if (grepl("^\\d+$", clean_id)) {
    type <- "pmid";  doi <- pmid_to_doi(clean_id)
  } else {
    stop("유효하지 않은 식별자 형식입니다. DOI, PMID, PMCID 형태를 확인해 주세요.")
  }

  if (is.null(doi) || is.na(doi) || doi == "") {
    stop(paste0("[", type, "] ", clean_id, " 에 해당하는 DOI를 찾을 수 없습니다."))
  }

  bib <- fetch_crossref_bibtex(doi)

  # Crossref BibTeX는 단체 저자(컨소시엄 등)를 빈칸으로 빼먹으므로 JSON에서 다시 만듦
  authors <- crossref_author_names(fetch_crossref_authors(doi), consortium)
  if (length(authors)) bib <- set_bib_field(bib, "author", paste(authors, collapse = " and "))

  if (key_type == "author_year") bib <- customize_bibkey(bib)
  bib <- format_bibtex(bib)

  class(bib) <- c("bibtex_str", "character")
  bib
}


#' 여러 식별자의 BibTeX를 .bib 파일로 저장
#'
#' 각 식별자에 `get_bibtex()`를 적용해 결과를 하나의 `.bib` 파일에 씁니다.
#' 실패한 식별자는 경고만 내고 건너뜁니다. 같은 citation key가 생기면
#' (예: 같은 저자의 같은 해 논문) 뒤에 `a`, `b`, ...를 붙여 구분합니다.
#'
#' @param ids `character`. DOI, PMID, PMCID를 섞어서 넣어도 됩니다.
#' @param file `character(1)`. 저장할 `.bib` 파일 경로.
#' @param key_type,consortium `get_bibtex()`의 같은 이름 인자와 같습니다.
#' @param append `logical(1)`. `TRUE`면 기존 파일 뒤에 추가합니다.
#'   이때 기존 파일의 key와도 겹치지 않게 처리합니다.
#'
#' @return 저장한 BibTeX 항목들(`character`)을 보이지 않게(invisible) 반환.
#'   실패한 식별자는 `attr(, "failed")`에 담깁니다.
#'
#' @examples
#' \dontrun{
#' ids <- c("10.1038/s41586-020-2649-2", "PMC7759461", "32939066")
#' save_bibtex(ids, "refs.bib")
#'
#' # 나중에 더 추가
#' save_bibtex("10.1093/nar/gks1195", "refs.bib", append = TRUE)
#' }
#'
#' @export
save_bibtex <- function(ids, file, key_type = c("author_year", "id"),
                        consortium = c("name_only", "members"), append = FALSE) {
  key_type   <- match.arg(key_type)
  consortium <- match.arg(consortium)
  ids <- unique(trimws(as.character(ids)))

  entries <- character(0)
  failed  <- character(0)
  for (id in ids) {
    if (id != ids[1]) Sys.sleep(1)   # API 요청 제한(429)을 피하려고 간격을 둠
    res <- tryCatch(get_bibtex(id, key_type, consortium), error = function(e) {
      warning(id, ": ", conditionMessage(e), call. = FALSE)
      NULL
    })
    if (is.null(res)) failed <- c(failed, id) else entries <- c(entries, unclass(res))
  }

  # 이미 파일에 있는 key, DOI 수집 (append일 때)
  used <- character(0); seen_doi <- character(0)
  if (append && file.exists(file)) {
    old  <- readLines(file, warn = FALSE, encoding = "UTF-8")
    hits <- regmatches(old, regexec("^\\s*@[A-Za-z]+\\{([^,]+),", old))
    used <- vapply(hits[lengths(hits) > 1], `[`, character(1), 2)
    dois <- regmatches(old, regexec("^\\s*doi\\s*=\\s*\\{([^}]+)\\}", old, ignore.case = TRUE))
    seen_doi <- tolower(vapply(dois[lengths(dois) > 1], `[`, character(1), 2))
  }

  # 같은 논문(같은 DOI)은 한 번만: PMID와 DOI로 같은 논문을 넣은 경우 등
  entry_doi <- vapply(entries, function(e) {
    s <- bib_field_span(e, "doi")
    if (is.null(s)) NA_character_ else tolower(substr(e, s[1], s[2]))
  }, character(1), USE.NAMES = FALSE)
  dup <- !is.na(entry_doi) & (duplicated(entry_doi) | entry_doi %in% seen_doi)
  if (any(dup)) message("중복 DOI ", sum(dup), "개는 건너뜀")
  entries <- entries[!dup]

  # 중복 key에 a, b, c ... 붙이기
  for (i in seq_along(entries)) {
    key <- sub("^@[A-Za-z]+\\{([^,]+),.*$", "\\1", entries[i])
    new_key <- key; k <- 0
    while (new_key %in% used) { k <- k + 1; new_key <- paste0(key, letters[k]) }
    if (new_key != key) {
      entries[i] <- sub("^(@[A-Za-z]+\\{)[^,]+,", paste0("\\1", new_key, ","), entries[i])
    }
    used <- c(used, new_key)
  }

  if (length(entries)) {
    con <- file(file, open = if (append) "a" else "w", encoding = "UTF-8")
    on.exit(close(con))
    writeLines(paste0(entries, "\n"), con)
  }

  message(length(entries), "개 저장 -> ", file,
          if (length(failed)) paste0(" (실패 ", length(failed), "개: ",
                                     paste(failed, collapse = ", "), ")"))
  attr(entries, "failed") <- failed
  invisible(entries)
}


#' PMCID를 DOI로 변환
#'
#' PMC ID Converter API를 호출해 PMCID에 대응하는 DOI를 찾습니다.
#' PMC에 수록된 논문만 변환됩니다.
#'
#' @param pmcid `character(1)`. `"PMC"` 접두사가 붙은 PMCID.
#' @return DOI 문자열, 찾지 못하면 `NULL`.
#' @keywords internal
#' @noRd
pmcid_to_doi <- function(pmcid) {
  url <- paste0("https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/?ids=",
                pmcid, "&format=json&tool=r_bibtex_fetcher")
  res <- tryCatch(jsonlite::fromJSON(read_url(url, tries = 3L)), error = function(e) NULL)
  if (!is.null(res) && "doi" %in% names(res$records)) res$records$doi[1] else NULL
}


#' PMID를 DOI로 변환
#'
#' PubMed E-utilities `esummary`의 `articleids`에서 DOI를 찾습니다.
#' PMC 수록 여부와 상관없이 PubMed에 있는 논문이면 동작합니다.
#'
#' @param pmid `character(1)`. 숫자로만 된 PMID.
#' @return DOI 문자열, 찾지 못하면 `NULL`.
#' @keywords internal
#' @noRd
pmid_to_doi <- function(pmid) {
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?",
                "db=pubmed&retmode=json&id=", pmid)
  res <- tryCatch(jsonlite::fromJSON(read_url(url, tries = 3L)), error = function(e) NULL)
  ids <- res$result[[pmid]]$articleids
  if (is.null(ids)) return(NULL)
  doi <- ids$value[ids$idtype == "doi"]
  if (length(doi)) doi[1] else NULL
}


#' DOI로 BibTeX 받기
#'
#' `https://doi.org/<doi>`에 `Accept: application/x-bibtex` 헤더로 요청해
#' 등록기관(Crossref, DataCite 등)이 제공하는 BibTeX를 받습니다.
#'
#' @param doi `character(1)`. `https://doi.org/` 접두사가 없는 DOI.
#' @return 앞뒤 공백을 제거한 BibTeX 문자열.
#' @keywords internal
#' @noRd
fetch_crossref_bibtex <- function(doi) {
  bib <- tryCatch(
    read_url(paste0("https://doi.org/", doi),
             headers = c(Accept = "application/x-bibtex; charset=utf-8")),
    error = function(e) stop("BibTeX 정보를 불러오는 데 실패했습니다: ", conditionMessage(e), call. = FALSE)
  )
  trimws(bib)
}


#' URL 내용을 문자열로 읽기 (재시도 포함)
#'
#' 실패하면(예: Crossref의 429 Too Many Requests) 2초, 4초, 8초 기다렸다가
#' 다시 시도합니다.
#'
#' @param url `character(1)`.
#' @param headers 이름 붙은 `character`. HTTP 요청 헤더.
#' @param tries `integer(1)`. 최대 시도 횟수.
#' @return 응답 본문 `character(1)`.
#' @keywords internal
#' @noRd
read_url <- function(url, headers = NULL, tries = 4L) {
  for (k in seq_len(tries)) {
    out <- tryCatch({
      con <- url(url, method = "libcurl", headers = headers)
      on.exit(close(con), add = TRUE)
      paste(suppressWarnings(readLines(con, warn = FALSE, encoding = "UTF-8")), collapse = "\n")
    }, error = function(e) e)
    if (!inherits(out, "error")) return(out)
    if (k < tries) Sys.sleep(2^k)
  }
  stop(conditionMessage(out), call. = FALSE)
}


#' Citation key를 저자성+연도 형식으로 변경
#'
#' `author` 필드의 첫 성과 `year` 필드를 읽어 key를 `harris2020` 형태로 바꿉니다.
#' 악센트 문자는 ASCII로 음역합니다(예: `Müller` -> `muller`).
#'
#' @param bib_str `character(1)`. BibTeX 항목 문자열.
#' @return key가 바뀐 BibTeX 문자열. 저자나 연도를 못 찾으면
#'   각각 `unknown`, `0000`을 씁니다.
#' @keywords internal
#' @noRd
customize_bibkey <- function(bib_str) {
  span  <- bib_field_span(bib_str, "author")
  first <- if (is.null(span)) "" else split_authors(substr(bib_str, span[1], span[2]))[1]
  first <- iconv(trimws(first), "UTF-8", "ASCII//TRANSLIT")
  if (is.na(first)) first <- ""

  if (startsWith(first, "{")) {
    # 단체 저자: 흔한 단어를 빼고 앞 두 단어 (The 1000 Genomes Project Consortium -> 1000genomes)
    words <- strsplit(gsub("[^A-Za-z0-9 ]", " ", first), "\\s+")[[1]]
    stop_words <- c("the", "of", "for", "and", "on", "consortium", "project", "group",
                    "collaboration", "collaborative", "investigators", "study", "network")
    words <- words[words != "" & !(tolower(words) %in% stop_words)]
    last  <- paste(utils::head(words, 2), collapse = "")
  } else {
    last <- gsub("[^A-Za-z]", "", sub(",.*$", "", first))
  }
  if (last == "") last <- "Unknown"

  y <- regmatches(bib_str, regexec("year\\s*=\\s*\\{?(\\d{4})", bib_str, ignore.case = TRUE))[[1]]
  year <- if (length(y) > 1) y[2] else "0000"

  sub("^\\s*@([A-Za-z]+)\\{[^,]+,", paste0("@\\1{", tolower(last), year, ","), bib_str)
}


#' Crossref JSON에서 저자 목록 받기
#'
#' `https://api.crossref.org/works/<doi>`의 `author` 배열을 받습니다.
#' Crossref에 없는 DOI(DataCite 등)면 `NULL`을 반환하고,
#' 그때는 doi.org가 준 BibTeX 저자를 그대로 씁니다.
#'
#' @param doi `character(1)`.
#' @return Crossref 저자 목록(list), 실패하면 `NULL`.
#' @keywords internal
#' @noRd
fetch_crossref_authors <- function(doi) {
  url <- paste0("https://api.crossref.org/works/", utils::URLencode(doi))
  res <- tryCatch(jsonlite::fromJSON(read_url(url, tries = 3L), simplifyVector = FALSE),
                  error = function(e) NULL)
  res$message$author
}


#' Crossref 저자 목록을 BibTeX 이름으로 변환
#'
#' 사람은 `성, 이름`, 단체는 BibTeX가 한 이름으로 보도록 `{단체 이름}`으로 씁니다.
#'
#' @param authors `fetch_crossref_authors()`의 결과.
#' @param consortium `"name_only"` 또는 `"members"`. `get_bibtex()` 참고.
#' @return BibTeX 이름 문자열 벡터. 저자 정보가 없으면 `character(0)`.
#' @keywords internal
#' @noRd
crossref_author_names <- function(authors, consortium = "name_only") {
  if (!length(authors)) return(character(0))

  is_org <- vapply(authors, function(a) is.null(a$family) && !is.null(a$name), logical(1))
  org_nm <- vapply(authors, function(a) if (is.null(a$name)) "" else a$name, character(1))
  names_ <- vapply(authors, function(a) {
    if (!is.null(a$family)) {
      if (!is.null(a$given)) paste0(a$family, ", ", a$given) else a$family
    } else if (!is.null(a$name)) {
      paste0("{", a$name, "}")
    } else NA_character_
  }, character(1))
  names_ <- gsub("(?<!\\\\)&", "\\\\&", names_, perl = TRUE)   # LaTeX에서 깨지는 & 처리

  # 제1저자가 단체 -> 컨소시엄 논문
  if (is_org[1] && consortium == "name_only") return(names_[1])

  # 단체 항목 중 "Steering committee", 기관명 같은 그룹 제목은 빼고
  # 컨소시엄/연구자 그룹 이름만 남김 (제1저자 단체는 항상 남김)
  team_re <- "consortium|collaborat|investigators|study group|working group|network|initiative"
  keep <- !is.na(names_) &
    (!is_org | seq_along(names_) == 1L | grepl(team_re, org_nm, ignore.case = TRUE))
  unique(names_[keep])
}


#' BibTeX 필드 값의 위치 찾기
#'
#' `field = {...}`에서 중괄호 안 값의 시작·끝 위치를 괄호 짝으로 찾습니다.
#'
#' @param bib_str `character(1)`. BibTeX 항목 문자열.
#' @param field `character(1)`. 필드 이름(대소문자 무관).
#' @return `c(start, end)` 정수 벡터, 필드가 없거나 괄호가 안 맞으면 `NULL`.
#'   값이 비어 있으면 `end = start - 1`.
#' @keywords internal
#' @noRd
bib_field_span <- function(bib_str, field) {
  m <- regexpr(paste0("(?<![A-Za-z])", field, "\\s*=\\s*\\{"), bib_str,
               ignore.case = TRUE, perl = TRUE)
  if (m == -1) return(NULL)
  start <- as.integer(m + attr(m, "match.length"))

  chars <- strsplit(substring(bib_str, start), "")[[1]]
  depth <- 1L
  for (i in seq_along(chars)) {
    if (chars[i] == "{") depth <- depth + 1L
    else if (chars[i] == "}") { depth <- depth - 1L; if (depth == 0L) return(c(start, start + i - 2L)) }
  }
  NULL
}


#' BibTeX 필드 값 바꾸기
#'
#' 필드가 있으면 값을 바꾸고, 없으면 key 바로 다음에 새로 넣습니다.
#'
#' @param bib_str `character(1)`. BibTeX 항목 문자열.
#' @param field,value `character(1)`. 필드 이름과 새 값(중괄호 제외).
#' @return 수정된 BibTeX 문자열.
#' @keywords internal
#' @noRd
set_bib_field <- function(bib_str, field, value) {
  span <- bib_field_span(bib_str, field)
  if (is.null(span)) {
    return(sub("^(\\s*@[A-Za-z]+\\{[^,]+,)", paste0("\\1 ", field, "={", value, "},"), bib_str))
  }
  paste0(substr(bib_str, 1L, span[1] - 1L), value, substring(bib_str, span[2] + 1L))
}


#' BibTeX 저자 필드를 저자별로 나누기
#'
#' 중괄호 밖의 ` and `에서만 나눕니다. `{Indian Telugu ... and Sri Lankan ...}`처럼
#' 단체 이름 안의 and는 나누지 않습니다.
#'
#' @param field `character(1)`. `author` 필드 값.
#' @return 저자 문자열 벡터.
#' @keywords internal
#' @noRd
split_authors <- function(field) {
  out <- strsplit(trimws(field), "(\\{(?:[^{}]++|(?1))*\\})(*SKIP)(*FAIL)|\\s+and\\s+", perl = TRUE)[[1]]
  out[trimws(out) != ""]
}


#' BibTeX를 필드당 한 줄로 정리
#'
#' doi.org가 주는 BibTeX는 항목 전체가 한 줄이라 읽기 어렵습니다.
#' 최상위 쉼표(중괄호·따옴표 밖)를 기준으로 필드를 나눠
#' 한 줄에 하나씩, 두 칸 들여쓰기로 다시 씁니다. 내용은 바꾸지 않습니다.
#'
#' @param bib_str `character(1)`. BibTeX 항목 문자열.
#' @return 정리된 BibTeX 문자열. 형식을 알 수 없으면 원본을 그대로 반환.
#' @keywords internal
#' @noRd
format_bibtex <- function(bib_str) {
  bib_str <- trimws(bib_str)
  head <- regmatches(bib_str, regexpr("^@[A-Za-z]+\\{[^,]+,", bib_str))
  if (!length(head)) return(bib_str)

  chars <- strsplit(substring(bib_str, nchar(head) + 1L), "")[[1]]
  depth <- 0L; in_quote <- FALSE
  cuts <- 0L; end <- NA_integer_
  for (i in seq_along(chars)) {
    ch <- chars[i]
    if (ch == "\"" && depth == 0L) in_quote <- !in_quote
    else if (in_quote) next
    else if (ch == "{") depth <- depth + 1L
    else if (ch == "}") {
      if (depth == 0L) { end <- i; break }       # 항목을 닫는 }
      depth <- depth - 1L
    } else if (ch == "," && depth == 0L) cuts <- c(cuts, i)
  }
  if (is.na(end)) return(bib_str)

  starts <- cuts + 1L
  stops  <- c(cuts[-1] - 1L, end - 1L)
  fields <- trimws(mapply(function(a, b) paste(chars[seq_len(b - a + 1L) + a - 1L], collapse = ""),
                          starts, stops))
  fields <- fields[fields != ""]
  fields <- sub("^([A-Za-z_-]+)\\s*=\\s*", "\\1 = ", fields)

  paste0(head, "\n", paste0("  ", fields, collapse = ",\n"), "\n}")
}


#' BibTeX 문자열 출력
#'
#' `get_bibtex()` 결과를 따옴표나 `\n` 이스케이프 없이 그대로 콘솔에 출력합니다.
#'
#' @param x `bibtex_str` 객체.
#' @param ... 사용하지 않음.
#' @return `x`를 보이지 않게(invisible) 반환.
#' @export
print.bibtex_str <- function(x, ...) {
  cat(x, sep = "\n")
  invisible(x)
}
