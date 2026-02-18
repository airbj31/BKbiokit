#' categorize the input column
#'
#' The function categorize the column based on the colum value and. condition/rule.
#'
#' @param data data frame
#' @param input_col input column
#' @param output_col output column (categorized column)
#' @param lvl level
categorize <- function(data, input_col, output_col, conditions, lvl=NULL, default = NA_character_) {

  col_data_type <- class(data[[input_col]])[1] # Use [1] to get the primary class (e.g., "character", "numeric")

  # 2. Validate that the rules are consistent with the column's data type
  lapply(conditions, function(rule) {
    if (!is.null(rule$op) && !is.null(rule$threshold)) {
      # This is a numeric rule
      if (!col_data_type %in% c("numeric", "integer", "double")) {
        stop(
          paste0("Numeric rule provided for non-numeric column '", input_col, "' (type: ", col_data_type, "). ",
                 "Numeric rules (with 'op' and 'threshold') can only be used with numeric columns."),
          call. = FALSE # 'call. = FALSE' makes the error message cleaner
        )
      }
    } else if (!is.null(rule$values)) {
      if (all(is.numeric(rule$values)) && !col_data_type %in% c("numeric", "integer", "double")) {
        warning(paste0("Categorical rule with numeric 'values' applied to non-numeric column '", input_col, "' (type: ", col_data_type, "). ",
                       "This might lead to unexpected behavior if values are not coerced correctly."), call. = FALSE)
      }
    }
  })
  case_expr <- build_case_when(input_col, conditions, default)
  data <- data |> dplyr::mutate(!!rlang::sym(output_col) := !!case_expr)
  if(!is.null(lvl)) {
    data <- data |> mutate(!!rlang::sym(output_col) := factor(!!rlang::sym(output_col),levels=lvl))
  }
  return(data)
}
