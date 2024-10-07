#' @name computeTVSpectralQuant
#' @title Compute Time-Varying Spectral Quantities
#' 
#' @description This function calculates spectral quantities either for an individual channel or cross-coherence 
#' between a pair of channels.
#' 
#' @param X numeric array, Data.
#' @param id Numeric Identifier:
#'        - If id is less than or equal to the number of columns (channels) of X, it computes the 
#'          auto-spectrum of the id-th channel.
#'        - If id is greater than the number of columns, it computes the cross-coherence between 
#'          the two channels corresponding to the id.
#' @param windowLen numeric Number of observations used for the estimation of localized spectral quantities.
#' @param normalize logical Indicates whether to normalize the spectral quantities.
#' @param overlap numeric Amount of overlap between successive windows.
#' @param M_welch numeric, Welch regularisation parameter.
#' @param logScale logical Indicates whether to use a logarithmic scale for the spectral quantities.
#' @param padEnd logical Specifies whether to pad time series with zeros at the end if the number of 
#'        observations is not a multiple of windowLen. Default is TRUE.
#' @param transform character Transformation to apply to the data.
#' @param plot logical Indicates whether to plot the resulting spectral quantities.
#' @param bands_analysis logical Indicates whether to analyze by frequency bands.
#' @param f_sampling numeric Sampling frequency.
#' 
#' @return array Depending on the id:
#'        - If id corresponds to a single channel, it returns the auto-spectrum of that channel.
#'        - If id corresponds to a pair of channels, it returns the cross-coherence between them.
#' 
computeTVSpectralQuant<-function (X, id, windowLen = windowLen, normalize = FALSE, overlap = 0, 
                                  M_welch = round_right(windowLen/20), logScale = FALSE, padEnd = TRUE, 
                                  transform = "FZ", plot = FALSE, bands_analysis = FALSE, f_sampling = NULL) 
{
  D <- ncol(X)
  if (id <= D) {
    x <- X[, id]
    S <- tv_spectrum(x, windowLen = windowLen, overlap = overlap, 
                     normalize = normalize, M_welch = M_welch, padEnd = padEnd, 
                     logScale = logScale, plot = plot, bands_analysis = bands_analysis, f_sampling = f_sampling,
                     method = method, M_dwelch = M_dwelch, p = p, h = h, nFF = nFF)$S
  }
  else {
    if (id < 1000) {
      a <- as.numeric(substr(id, 1, 1))
      b <- as.numeric(substr(id, 2, 3))
    }
    else {
      a <- as.numeric(substr(id, 1, 2))
      b <- as.numeric(substr(id, 3, 4))
    }
    x <- cbind(X[, a], X[, b])
    S <- tv_coherence(x, windowLen = windowLen, overlap = overlap, 
                      M_welch = M_welch, logScale = logScale, padEnd = padEnd, 
                      transform = transform, plot = plot, bands_analysis = bands_analysis, f_sampling = f_sampling)$coh
  }
  return(S)
}
