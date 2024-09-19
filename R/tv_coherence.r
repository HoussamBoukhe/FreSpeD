#' @name tv_coherence
#' @title time varying coherence of two signals 
#' @description Estimate the cross coherence of two signals
#' 
#' @param X numeric two-dimensional data input, rows: observations in time, columns: channels.
#' @param windowLen numeric Number of observations used for the estimation of localized spectral quantities.
#' @param overlap numeric, the overlap between consecutive windows (default = 0)
#' @param M_welch numeric, Welch regularisation parameter.
#' @param logScale logical, Indicates whether to use a logarithmic scale to the estimated coherence
#' @param padEnd logical, Pad time series with zeros at end (if number of observations in time is not a multiple of windowLen). Default is TRUE
#' @param tranform character, Apply a tranformation to the coherence (default = 'FZ')
#'  'FZ' Fisher-z transformation, 'ARC' arcsin the inverse of sin function, '' no transformation is applied
#' @param plot logical Indicates whether to plot the estimated cross-coherence
#' @param bands_analysis logical, if TRUE, the function will make the analysis by bands
#' @param f_sampling numeric, sampling frequency
#' 
#' @return list, Output a list containing the estimated cross-coherence, the set of windows used in the estimation, and the applied transform.
tv_coherence<-function (X, windowLen, overlap = 0, M_welch = windowLen/2 * 0.1, 
                        logScale = FALSE, padEnd = TRUE, transform = "FZ", plot = FALSE,bands_analysis = FALSE, f_sampling = NULL) 
{
  if (is.null(ncol(X))) 
    stop("Two-column array of input data required")
  else if (ncol(X) != 2) 
    stop("Two-column array of input data required")
  if (!is.numeric(X)) 
    stop("Data not numeric")
  if (windowLen < 1 || windowLen%%1 != 0) 
    stop("window length has to be a positive interger")
  if (overlap < 0 || overlap%%1 != 0) 
    stop("Overlap has to be a non-negative interger")
  if (windowLen > nrow(X)) 
    stop("Window length has to be smaller than number of observations")
  #if ((windowLen/2/M_welch) < 4 || (windowLen/2/M_welch) > (windowLen/2)) # ???... 
  if ((windowLen/2/M_welch) < 4 || (M_welch) > (windowLen/2))
    stop("Number of frequency blocks badly specified: frequencies per band have to exceed 4 and cannot exceed 0.5 window length")
  
  fWindowLen <- floor(windowLen/M_welch)

  smallno <- 1e-10
  
  # We do similar steps as in the function tv_spectrum 
  Y1 <- partition(X[, 1], windowLen, overlap, padEnd)
  Y2 <- partition(X[, 2], windowLen, overlap, padEnd)
  Y1[is.na(Y1)] <- 0
  Y2[is.na(Y2)] <- 0
  Y1 <- sweep(Y1, 2, apply(Y1, 2, mean), "-")
  Y2 <- sweep(Y2, 2, apply(Y2, 2, mean), "-")

  for (iWin in 1:ncol(Y1)) {
          varY1 <- var(Y1[, iWin])
          if (varY1 < smallno) varY1 <- smallno
          Y1[, iWin] <- Y1[, iWin] / sqrt(varY1)
      }

  for (iWin in 1:ncol(Y2)) {
        varY2 <- var(Y2[, iWin])
        if (varY2 < smallno) varY2 <- smallno
        Y2[, iWin] <- Y2[, iWin] / sqrt(varY2)
    }
  
  if (windowLen%%2 == 1) {
    Y1 <- rbind(Y1, 0)
    Y2 <- rbind(Y2, 0)
  }
  
  windowLen <- nrow(Y1)
  # nF <- floor(windowLen/2)  

  coh <- array(NA, c((fWindowLen/2), ncol(Y1)))
  
  for (iInterval in 1:ncol(Y1)) {
    xfft <- array(NA, c(M_welch, fWindowLen, 2))
    for (iFB in 1:M_welch) {
      Y1int <- Y1[((iFB * fWindowLen) - fWindowLen + 1):(iFB * fWindowLen), iInterval] 
      Y1int <- Y1int - mean(Y1int)

      Y2int <- Y2[((iFB * fWindowLen) - fWindowLen + 1):(iFB * fWindowLen), iInterval] 
      Y2int <- Y2int - mean(Y2int)

      # Combine the two signals, then compute the FFT of each signal, and return the result by combining the transforms as columns
      xfft[iFB, , ] <- mvfft(cbind(Y1int, Y2int)) 
      }
      
    for (iF in 1:(fWindowLen/2)) { # ... Computation of cross coherence between the two signals 
      x1x2 <- xfft[, iF, 1] * Conj(xfft[, iF, 2])
      denom1 <- mean(abs(xfft[, iF, 1])^2)
      denom2 <- mean(abs(xfft[, iF, 2])^2)
            
      if (denom1 < smallno) denom1 <- smallno
      if (denom2 < smallno) denom2 <- smallno
            
      coh[iF, iInterval] <- Mod(mean(x1x2))^2 / (denom1 * denom2)
      
    }
  }
  coh[coh > 1] <- 1

  if (transform == "FZ"){
    # Before applying Fisher transform, we change the values equal to 1 to 1 - epsilon
    coh[coh > 1 - smallno] <- 1 - smallno
    coh <- 1/2 * log((1 + sqrt(coh))/(1 - sqrt(coh)))}
  if (transform == "ARC") {
    coh[coh > 1 - smallno] <- 1 - smallno
    coh <- asin(abs(coh))}

  if (bands_analysis){
      
      coh_bands <- array(NA, dim = c(5, ncol(coh)))

      if  (nrow(coh) <= f_sampling/4){stop("Choose a smaller number of Welch regulariser or a larger window length to be able to analyse the bands")}

      indices <- get_bands_indices(2*nrow(coh), f_sampling)
      
      coh_bands[1, ] <- apply(coh[indices$delta_indices, ], 2, mean)
      coh_bands[2, ] <- apply(coh[indices$theta_indices, ], 2, mean)
      coh_bands[3, ] <- apply(coh[indices$alpha_indices, ], 2, mean)
      coh_bands[4, ] <- apply(coh[indices$beta_indices, ], 2, mean)
      coh_bands[5, ] <- apply(coh[indices$gamma_indices, ], 2, mean)
      

      coh <- coh_bands
    }

  if (plot) {
  df <- reshape2::melt(coh)
  
  # Define custom labels for x-axis
  custom_labels <- seq(round_right(windowLen/2), length(Y1), by = windowLen)
  
  # Define custom labels for y-axis
  custom_labels_y <- seq(dim(coh)[1])
  
  # Plot heatmap with custom x-axis and y-axis labels
  p <- ggplot(df, aes(x = Var2, y = factor(Var1, levels = seq(dim(coh)[1]), labels = custom_labels_y), fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black") +  # Define colors
    labs(x = "time", y = "Fundamental frequencies", title = "Cross-coherence") +  # Add labels
    theme_minimal() +                                    # Choose a theme
    guides(fill = guide_colorbar(title = "Value")) +    # Add color scale legend
    scale_x_continuous(breaks = seq(dim(coh)[2]), labels = custom_labels) # Custom x-axis labels
  
  print(p)
}

  out <- list(coh = coh, data = list(Y1, Y2), "transform" == transform)
}