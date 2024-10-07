#' @name partition
#' @title Partition
#' @description Partition the series of numbers into windows, and stack them vertically
#'
#' @param X Numeric vector, the series to partition
#' @param windowLen Numeric, the length of the window
#' @param overlap Numeric, the overlap between consecutive windows (default = 0)
#' @param padEnd Logical, Pad time series with zeros at end (if number of observations in time is not a
#' multiple of windowLen). Default is TRUE
#' @return Array, a 2D array where each column contains a partitioned window
#' @examples
#' series <- c(1:20)
#' partitioned <- partition(series, windowLen = 3, overlap = 1)
partition<-function (X, windowLen, overlap = 0, padEnd = TRUE) 
{
    nX <- length(X)
    winStart <- seq(1, nX - 1, by = windowLen - overlap)
    nWindow <- length(winStart)

    if (padEnd) {
        if ((nX - tail(winStart, 1)) != windowLen) {
            X <- c(X, rep(NA, (windowLen - nX + tail(winStart, 1)))) }
    } else { # If padEnd is FALSE, we throw the last part of the series that does not give a window of length windowLen
        if (tail(winStart, 1) > (nX - windowLen + overlap + 1)) {
            nWindow <- nWindow - 1
            winStart <- winStart[1:nWindow]
        }
    }


    Y <- array(0, c(windowLen, nWindow))
    for (iWindow in 1:nWindow){
        Y[1:windowLen, iWindow] <- X[winStart[iWindow]:(winStart[iWindow] + windowLen - 1)]
    }
    return(Y)
}
