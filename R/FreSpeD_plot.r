#' @name FreSpeD_plot
#' @title Plots the time-varying spectrum and cross-coherence
#' @description Plots the time-varying spectrum for a specified channel or coherence for a channel pair 
#' along with the estimated change points in time and frequency provided by FreSpeD.
#' 
#' @param X numeric array, Data
#' @param cp list of numeric arrays, List of change points per channel/channel pair given by FreSpeD function
#' @param id numeric, Which list entry should be plotted?
#' @param channelNames numeric vector, channel names. Default is 1:dim(X)[2]
#' @param windowLen numeric, Number of observations used for the estimation of localized spectral quantities
#' @param nFB numeric, Welch regularisation parameter.
#' @param overlap numeric, Number of observations by which localized spectral estimates can overlap. Default is zero
#' @param normalize logical, Normalize spectrum (dividing by total localized variance)?
#' @param padEnd logical, Pad time series with zeros at end (if number of observations in time is not a multiple of windowLen).
#' 
#' @return no output is returned 
#' 
FreSpeD_plot<-function (X, cp, id, channelNames = 1:dim(X)[2], windowLen, nFB = round_right(windowLen/20), 
    overlap = 0, normalize = FALSE, padEnd = TRUE) 
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
    logScale <- FALSE
    transform <- "FZ"
    plot <- FALSE
    check_neighborhood <- TRUE
    Tv <- ifelse(padEnd, ceiling(T/windowLen), floor(T/windowLen))
    idlist <- FreSpeD::define_identifier(D)
    if (!(id %in% idlist)) 
        stop("Id invalid")
    S <- FreSpeD::computeTVSpectralQuant(X, id, windowLen = windowLen, 
        normalize = normalize, overlap = overlap, nFB = nFB, 
        logScale = logScale, padEnd = padEnd, transform = transform, 
        plot = plot)
    freqSpec <- TRUE
    if (length(cp) == (D + D * (D - 1)/2)) {
        cph <- cp[[which(idlist == id)]]
        if (is.null(dim(cph))) {
            cph <- array(0, dim(S))
            cph[, cp[[which(idlist == id)]]] = 1
            freqSpec <- FALSE
        }
    }
    else {
        cph <- array(0, dim(S))
        cph[, round_right(cp)] = 1
        freqSpec <- FALSE
    }
    L <- dim(S)[1]
    try(dev.new(noRStudioGD = TRUE), silent = TRUE)
    layout(matrix(c(1, 2, 3, 0, 0, 4), nrow = 2, ncol = 1), widths = c(4), 
        heights = c(4, 1))
    pal.1 = colorRampPalette(c("yellow", "red", "black"), space = "rgb")
    channelN <- FreSpeD::name_identifiers(channelNames)[id]
    if (id <= D) 
        print(paste("Time-varying spectrum of channel", channelN))
    else print(paste("Time-varying coherence of channel pair", 
        channelN))
    print("Yellow (small) to black (large), x-axis: time, y-axis: frequency")
    if (freqSpec) 
        print("Circle indicates change-point location in time and frequency")
    else print("Circle indicates change-point location in time, frequency dimension is not carried over through postprocessing.")
    breaks <- seq(min(S), max(S), length.out = 100)
    par(mar = c(2, 2, 1, 1))
    image(seq(dim(S)[2]), seq(dim(S)[1])/L, t(S), col = pal.1(length(breaks) - 
        1), breaks = breaks, xaxs = "i", xaxt = "n")
    for (i in which(colSums(cph) > 0)) {
        v <- which(cph[, i] > 0)
        for (j in 1:length(v)) points(i, v[j]/L, pch = 1)
    }
    axis(1, at = seq(dim(S)[2]), labels = seq(round_right(windowLen/2), 
        T, by = windowLen))
    par(mar = c(1.9, 1, 1, 1))
    col_scale(S, col = pal.1(length(breaks) - 1), breaks = breaks, 
        horiz = TRUE)
}