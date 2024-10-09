#' @name tv_estimations_cp
#' @title time-varying estimation of change points
#' @description time-varying estimation of change points. This function computes the localized spectral quantities; spectrums and cross-coherences,
#' then applies the cusum_stat function to these quantities.
#' 
#' @param X numeric array, rows: observations in time, columns: channels.
#' @param channelNames list of channel names; if not specified, channels are enumerated.
#' @param windowLen numeric, Number of observations used for the estimation of localized spectral quantities.
#' @param M_welch numeric, Welch regularisation parameter. Default is windowLen/20
#' @param overlap numeric, Number of observations by which localized spectral estimates can overlap. Default is zero.
#' @param normalize logical, Normalize spectrum (divide by total localized variance)
#' @param padEnd logical, Pad time series with zeros at end (if number of observations in time is not a multiple of windowLen). Default is TRUE
#' @param logScale logical, Indicates whether to use a logarithmic scale
#' @param transform logical, Transformation of coherence estimate. Fisher-z transform is default and provides consistent change-point estimates.
#' @param plot logical, Plot all autospectra and cross-coherences
#' @param check_neighborhood logical, Within the change-point detection, should the neighbourhood around a changepoint 
#' candidate be tested to reduce the risk of spurious detections
#' @param check_par numeric, \Delta_T the diameter of the interval around a candidate change-point used in check_neighborhood
#' @param threshFun function, Threshold function; see reference for theoretical requirements
#' @param deltaPer numeric, Percentage (of total observations in time) of minimum change-point distance.
#' @param nCores numeric, the number of cores used in the computation
#' @param bands_analysis logical, Should the bands analysis be performed
#' @param f_sampling numeric, sampling frequency of the time series
#' 
#' @return list of arrays, a list containing the the values of the CUSUM statistics.
#' 
#' @references 
#' Shroder & Ombao
#' 
tv_estimations_cp<-function (X, channelNames = 1:dim(X)[2], windowLen, M_welch = round_right(windowLen/20), 
    overlap = 0, normalize = FALSE, padEnd = TRUE, logScale = FALSE, 
    transform = "FZ", plot = FALSE, check_neighborhood = TRUE, check_par = NaN,
    threshFun = function(Tv) {
        0.8 * log(Tv)^(1.1)
    } ,deltaPer = 0.03, nCores = 1,bands_analysis = FALSE, f_sampling = NULL) 
{
    
    if (is.null(dim(X))) 
        stop("Requires dim(X)!=NULL")
    T <- nrow(X)
    D <- ncol(X)
    if (D > T) 
        warning("Number of channels D exceeds number of observations T")
    if (windowLen > (T/50)) 
        warning("Window length exceeds 2%T")
    if (windowLen < 100) 
        warning("Window length is small. Default number of frequency bands less than 5.")
    if (overlap < 0) 
        stop("Overlap must be positive")
    Tv <- ifelse(padEnd, ceiling(T/windowLen), floor(T/windowLen))
    delta <- round_right(deltaPer * T/windowLen)
    thresh <- threshFun(T/windowLen)
    
    if (is.nan(check_par)){check_par <- floor(log(Tv)^2.1/5)}
    
    idlist <- define_identifier(D)
    cp <- list() # an empty list, that will contain the change points for each pairs of channels
    if (nCores == 1) { 
        cnt <- 1
        for (id in (idlist)) {
            S <- computeTVSpectralQuant(X, id, windowLen = windowLen, 
                normalize = normalize, overlap = overlap, M_welch = M_welch, 
                logScale = logScale, padEnd = padEnd, transform = transform, 
                plot = plot, bands_analysis = bands_analysis, f_sampling = f_sampling)

            cp[[cnt]] <- cusum_stat(S, thresh = thresh, check_neighborhood = check_neighborhood, 
                check_par = check_par, delta = delta)
            cnt <- cnt + 1
        }
    }

    else { 
        # Check wether the packages "foreach" and "doParallel" are available
        t1 <- !requireNamespace("foreach", quietly = TRUE)
        t2 <- !requireNamespace("doParallel", quietly = TRUE)
        if (t1 || t2) 
            stop("Packages 'foreach' and 'doParallel' are required to use multiple cores. Set nCores=1 to continue or install packages and dependencies.")
        # registers a parallel backend, which allows parallel computing 
        registerDoParallel(cores = nCores)

        cp <- foreach(id = idlist, .inorder = TRUE) %dopar% {
            
            S <- FreSpeD:::computeTVSpectralQuant(X, id, windowLen = windowLen, 
                normalize = normalize, overlap = overlap, M_welch = M_welch, 
                logScale = logScale, padEnd = padEnd, transform = transform, 
                plot = plot, bands_analysis = bands_analysis, f_sampling = f_sampling)
            FreSpeD:::cusum_stat(S, thresh = thresh, check_neighborhood = check_neighborhood, 
                check_par = check_par, delta = delta)
        }
        # Close the parallel backend once the computation are finished
        stopImplicitCluster() 
    }

    names(cp) <- name_identifiers(channelNames) # Name each matrix that represent a pair of channels
    for (i in 1:length(cp)) colnames(cp[[i]]) <- paste("TimeInterval", 
        1:dim(cp[[i]])[2], sep = "")
    for (i in 1:length(cp)) rownames(cp[[i]]) <- paste("FreqBand", 
        1:dim(cp[[i]])[1], sep = "")
    return(cp)
}