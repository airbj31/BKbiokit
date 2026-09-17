#' Categorize a Column Based on Conditions
#'
#' Creates a new categorical column in a data frame by applying a set of
#' ordered conditions (rules) to an existing column. Supports both numeric
#' comparison rules and categorical value-matching rules. Optionally converts
#' the result to an ordered factor.
#'
#' @param data A data frame containing the column to categorize.
#' @param input_col A character string naming the column to evaluate.
#' @param output_col A character string naming the new categorized column to create.
#' @param conditions A named list of rule objects. Each rule is itself a list with
#'   either:
#'   \itemize{
#'     \item \code{op}, \code{threshold}, and \code{label} for numeric comparisons
#'       (e.g. \code{list(op = ">", threshold = 50, label = "High")}), or
#'     \item \code{values} and \code{label} for categorical matching
#'       (e.g. \code{list(values = c("Y", "Yes"), label = "Positive")}).
#'   }
#'   Rules are evaluated in order; the first match wins.
#' @param lvl An optional character vector of factor levels. If provided, the
#'   output column is coerced to a \code{\link{factor}} with levels in the given
#'   order. Useful for downstream plotting or modelling that requires ordered
#'   categories. Defaults to \code{NULL} (plain character column).
#' @param default The value assigned when no condition is matched. Defaults to
#'   \code{NA_character_}.
#'
#' @return The input \code{data} frame with one additional column named
#'   \code{output_col}. If \code{lvl} is supplied the column is a factor,
#'   otherwise a character vector.
#'
#' @details
#' Internally, the conditions are compiled into a \code{\link[dplyr]{case_when}}
#' expression via a helper \code{build_case_when()}. Two validation checks are
#' performed before mutation:
#' \itemize{
#'   \item A numeric rule (\code{op} + \code{threshold}) applied to a non-numeric
#'     column raises an \strong{error}.
#'   \item A categorical rule whose \code{values} are numeric but \code{input_col}
#'     is non-numeric raises a \strong{warning}, as implicit coercion may produce
#'     unexpected results.
#' }
#'
#' @examples
#' df <- data.frame(score = c(30, 55, 80, NA))
#'
#' rules <- list(
#'   list(op = "<",  threshold = 40, label = "Low"),
#'   list(op = "<=", threshold = 60, label = "Mid"),
#'   list(op = ">",  threshold = 60, label = "High")
#' )
#'
#' # Basic numeric categorisation
#' categorize(df, "score", "grade", conditions = rules, default = "Unknown")
#'
#' # With ordered factor output
#' categorize(df, "score", "grade", conditions = rules,
#'            lvl = c("Low", "Mid", "High"), default = "Unknown")
#'
#' @seealso \code{\link[dplyr]{case_when}}, \code{\link[dplyr]{mutate}}
#'
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#'
#' @export

categorize <- function(data, input_col, output_col, conditions, lvl = NULL, default = NA_character_) {
  
  # 1. Validate core inputs
  stopifnot(is.data.frame(data))
  stopifnot(is.character(input_col), length(input_col) == 1, input_col %in% names(data))
  stopifnot(is.character(output_col), length(output_col) == 1)
  
  # 2. Validate conditions
  if (is.null(conditions) || length(conditions) == 0)
    stop("'conditions' must be a non-empty list.", call. = FALSE)
  
  # 3. Validate rules against column data type
  col_data_type <- class(data[[input_col]])[1]
  
  for (rule in conditions) {
    if (!is.null(rule$op) && !is.null(rule$threshold)) {
      # Numeric rule applied to non-numeric column
      if (!col_data_type %in% c("numeric", "integer", "double")) {
        stop(
          paste0("Numeric rule provided for non-numeric column '", input_col,
                 "' (type: ", col_data_type, "). ",
                 "Numeric rules (with 'op' and 'threshold') can only be used with numeric columns."),
          call. = FALSE
        )
      }
    } else if (!is.null(rule$values)) {
      # Categorical rule with numeric values on non-numeric column
      if (all(is.numeric(rule$values)) && !col_data_type %in% c("numeric", "integer", "double")) {
        warning(
          paste0("Categorical rule with numeric 'values' applied to non-numeric column '",
                 input_col, "' (type: ", col_data_type, "). ",
                 "This might lead to unexpected behavior if values are not coerced correctly."),
          call. = FALSE
        )
      }
      # Numeric column with categorical (character) values rule
      if (all(is.character(rule$values)) && col_data_type %in% c("numeric", "integer", "double")) {
        warning(
          paste0("Categorical rule with character 'values' applied to numeric column '",
                 input_col, "' (type: ", col_data_type, "). ",
                 "Consider using a numeric rule with 'op' and 'threshold' instead."),
          call. = FALSE
        )
      }
    }
  }
  
  # 4. Validate lvl against condition labels
  if (!is.null(lvl)) {
    if (!is.character(lvl))
      stop("'lvl' must be a character vector.", call. = FALSE)
    labels <- sapply(conditions, `[[`, "label")
    unmatched <- setdiff(lvl, c(labels, default))
    if (length(unmatched) > 0)
      warning(
        "Some 'lvl' values don't match any condition label: ",
        paste(unmatched, collapse = ", "), call. = FALSE
      )
  }
  
  # 5. Build and apply case_when expression
  case_expr <- build_case_when(input_col, conditions, default)
  data <- data |> dplyr::mutate(!!rlang::sym(output_col) := !!case_expr)
  
  # 6. Optionally coerce output to factor
  if (!is.null(lvl)) {
    data <- data |>
      dplyr::mutate(!!rlang::sym(output_col) := factor(!!rlang::sym(output_col), levels = lvl))
  }
  
  return(data)
}
