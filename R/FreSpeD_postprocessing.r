#' @name FreSpeD_postprocessing
#' @title Optional postprocessing
#' @description Optional postprocessing to account for information of channel clusters 
#' 
#' @param X numeric array, Data 
#' @param cp list of numeric arrays, List of change points per channel/channel pair given by FreSpeD function 
#' @param singleTimeSeries logical, Check individual channels and channel pairs? This tests every change point of a
#'        given channel/channel pair versus its adjacent change-point neighbors. The test
#'        statistic is consistent with the initial thresholded CUSUM statistic.
#' @param allTimeSeries logical, Merge change points over all channels and channel pairs into global clusters. 
#'        A cluster is defined as a group of change points where each is at most clusterradius
#'        distant from the nearest cluster member. The resulting merged change point is
#'        the median of all cluster members
#' @param clusterradius list Maximum distance allowed between members of the clusters defined through
#'        allTimeSeries(default is empty list)
#' @param overlap numeric, the overlap between consecutive windows (default = 0)
#' @param normalize logical Indicates whether to normalize the spectral quantities.
#' @param logScale logical, Indicates whether to use a logarithmic scale to spectral quantities 
#' @param padEnd logical, logical, Pad time series with zeros at end (if number of observations in time is not a multiple of windowLen). Default is TRUE
#' @param transform character, transformation of the estimate. Fisher-z transform is default
#' @param thresFun function, Threshold function
#' 
#' @return If only individual tests are conducted, the output is a list similar to that of FreSpeD(). Otherwise,
#'         a single vector of merged change points is provided.
FreSpeD_postprocessing<-function (X, cp, windowLen, singleTimeSeries = TRUE, allTimeSeries = TRUE, 
    clusterradius = c(), deltaPer = 0.03, nFB = round_right(windowLen/20), 
    overlap = 0, normalize = FALSE, logScale = FALSE, padEnd = TRUE, 
    transform = "FZ", threshFun = function(Tv) {
        0.8 * log(Tv)^(1.1)
    }) 
{
    T <- nrow(X)
    D <- ncol(X)
    # cp is a list of matrices (the estimated spectrograms given by FreSpeD), we take the sum of the columns of each matrix, then takes the 
    # indices of the columns which have a positive sum
    cp <- lapply(cp, function(x) {
        which(colSums(x) > 0) # Now, cp is a list of vectors containing the indices of the columns that have a positive sum
    })
    delta <- round_right(deltaPer * T/windowLen)
    Tv <- ifelse(padEnd, ceiling(T/windowLen), floor(T/windowLen)) # number of the windows 
    delta <- round_right(deltaPer * T/windowLen)
    thresh <- threshFun(T/windowLen)
    check_par <- floor(log(Tv)^2.1/5)
    idlist <- define_identifier(D)
    if (singleTimeSeries) {
        cnt <- 1 # A counter used to enumerate the elements of idlist from 1 to length(idlist)
        for (id in idlist) {
            cpTmp <- cp[[cnt]] # temporary cp (depending on id)
            totest <- c(1, cpTmp, Tv) # A sub-window of a window
            s <- 1
            while (length(cpTmp) > 0 && totest[s + 1] != Tv) { # ???... check again for spurious change point ? (taking the frequencies together...)
                cp_s <- totest[s]
                cp_b <- totest[s + 1]
                cp_e <- totest[s + 2]
                S <- FreSpeD::computeTVSpectralQuant(X, id, windowLen = windowLen, 
                  normalize = normalize, overlap = overlap, nFB = nFB, 
                  logScale = logScale, padEnd = padEnd, transform = transform, 
                  plot = FALSE)
                ind_stop <- FreSpeD::FreSpeD_innerpost(S[, cp_s:cp_e], 
                  thresh = thresh, check_neighborhood = TRUE, 
                  check_par = check_par, b = cp_b - cp_s + 1, 
                  delta = delta/2)
                if (!ind_stop) { # if the point is spurious point
                  totest <- totest[-(s + 1)] # remove this point
                }
                else s <- s + 1
            }
            cp[[cnt]] <- totest[-c(1, length(totest))] # remove the first and the last element 
            cnt <- cnt + 1
        }
    }
    if (allTimeSeries) { 
        delta <- ifelse(length(clusterradius) == 0, round_right(delta/2), 
            clusterradius)
        lcp <- unlist(cp) # put the the change point windows index in the same vector (list change point)
        ucp <- sort(unique(lcp))
        # Computing the number of clusters, if two change point are close (have distance less than delta) then 
        # they are in the same cluster (unique change point)
        if (length(ucp) > 1) { 
            ucp <- cbind(ucp, 1)
            cl <- 1
            for (i in 2:nrow(ucp)) {
                if (ucp[i, 1] - ucp[i - 1, 1] > delta) {
                  cl <- cl + 1
                  ucp[i:nrow(ucp), 2] <- cl
                }
            }
            ncl <- cl # number of clusters 
            m <- c()
            for (cl in 1:ncl) {
                s <- head(ucp[ucp[, 2] == cl, 1], 1) # get the first element of the cluster cl 
                e <- tail(ucp[ucp[, 2] == cl, 1], 1) # get the last element of the cluster cl
                m <- c(m, median(lcp[intersect(which(lcp >= s), 
                  which(lcp <= e))])) # get all the change point points in the cluster cl, then take the median
            }
            print(paste("Returning change points merged over all ", 
                length(cp), " components and component pairs provided. Found ", 
                ncl, ifelse(ncl > 1, " clusters", " cluster"), 
                sep = ""))
            cp <- m
        }
    }
    return(cp) # return the medians of the clusters 
}