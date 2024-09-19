#' @name round_right
#' @title Round right
#' @description Round to the nearest integer towards positive infinity
#'
#' @param x Numeric vector, the input vector
#' @return Numeric vector, rounded values
#' @examples
#' round_right(c(1.1, 2.3, 3.5, 4.7, 5.9))
round_right <- function(x) {
    return(trunc(x + 0.5))
}



