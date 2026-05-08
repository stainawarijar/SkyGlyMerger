#' Version of `isTruthy` that works with non-existing reactive expressions
#'
#' When using `isTruthy` on a reactive expression but that does not yet exist,
#' an error occurs. To avoid that, this function first checks if the expression
#' x can be called before checking its "truthiness".
#'
#' @param x An expression for which to test the "truthiness".
#'
#' @return FALSE if calling the expression results in an error. Otherwise, the
#'  output of `isTruthy(x)` which is either TRUE or FALSE.
is_truthy <- function(x) {

  valid <- tryCatch(
    expr = {
      x
      TRUE
    },
    error = function(e) {
      FALSE
    }
  )

  if (valid) {
    return(isTruthy(x))
  }

  return(FALSE)
}
