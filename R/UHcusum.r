#' @name UHcusum
#' @title CUSUM statistic
#' @description Calculate the cumulative sum of a series using the UH cusum method
#' This function calculates the cumulative sum of a series of numbers
#' using the UH cusum method.
#' @param X Numeric vector, series of numbers
#' @return Numeric vector, cumulative sum using the UH cusum method
UHcusum <- function(X) {
    n <- length(X)
    Y <- array(NA, (n - 1))
    for (i in 1:(n - 1)) {
        Y[i] <- sqrt((n - i) / (n * i)) * sum(X[1:i]) - sqrt(i / (n * (n - i))) * sum(X[(i + 1):n])
    }
    return(abs(Y))
}