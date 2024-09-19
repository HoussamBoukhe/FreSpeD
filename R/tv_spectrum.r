#' @name tv_spectrum
#' @title time-varying spectrum of a univariate time series
#' @description Estimate the auto-spectrum of a univarite time series 
#' 
#' @param X numeric vector, a signal 
#' @param windowLen numeric, the length of the windows used in the partitin of time domain
#' @param overlap numeric, wether the windows overlap or not 
#' @param M_welch numeric, Welch regularisation parameter. Default is 20% of the windowLen
#' @param normalize logical Indicates whether to normalize the spectral quantities.
#' @param padEnd logical, Pad time series with zeros at end (if number of observations in time is not a multiple of windowLen). Default is TRUE
#' @param logScale logical, Indicates whether to use a logarithmic scale to the estimated periodogram
#' @param bands_analysis logical, if TRUE, the analysis is performed by bands
#' @param f_sampling numeric, sampling frequency
#' @param plot logical Indicates whether to plot the estimated spectrum.
#' 
#' 
#' @return list, output a list of length 2, the first element is the estimated auto-spectrum, and 
#' the second, the the set of windows used in the estimation
tv_spectrum<-function (X, windowLen, overlap = 0, M_welch = windowLen/2 * 0.1,
                       normalize = FALSE, padEnd = TRUE, logScale = FALSE, bands_analysis = FALSE, f_sampling = NULL, plot = FALSE, ...) 
{
  if (!is.numeric(X)) 
    stop("Data not numeric")
  if (windowLen < 1 || windowLen%%1 != 0) 
    stop("window length has to be positive and interger")
  if (overlap < 0 || overlap%%1 != 0) 
    stop("Overlap has to be non-negative and interger")
  if (!is.null(dim(X))) 
    if (ncol(X) > 1) 
      stop("Data has to be univariate")
  if (windowLen > length(X)) 
    stop("Window length has to be smaller than number of observations")   
  if ((windowLen/2/M_welch) < 4 || (M_welch) > (windowLen/2) || M_welch < 1) 
    stop("Welch regularisation number is badly specified: frequencies per band have to exceed 4 and cannot exceed 0.5 window length")

  if (bands_analysis){
    if (is.null(f_sampling)){stop("Please specify a sampling frequency")}
    }  

   
  
  smallno <- 1e-10 
  
  Y <- partition(X, windowLen, overlap, padEnd) # Partition the time domain into windows, then stuck them verticaly into a matrix Y
  Y[is.na(Y)] <- 0
  
  Y <- sweep(Y, 2, apply(Y, 2, mean), "-") # Center each column of the matrix Y
  for (iWin in 1:ncol(Y)){
    var_i <- var(Y[, iWin])
    if (var_i == 0){ var_i <- smallno} 
    
    Y[, iWin] <- Y[, iWin]/sqrt(var_i) # Normalise each column of the matrix Y

    }

  
  # If the length of the window is odd
    if (windowLen %% 2 == 1) {
      # Add a row of zeros at the end of the matrix Y
      Y <- rbind(Y, 0)
    }

  # The length of each window is now even
  windowLen <- nrow(Y)

  fWindowLen <- floor(windowLen/M_welch) # Number of possible frequencies in each window 
 
  S <- array(NA, c((fWindowLen/2), ncol(Y))) # Create an empty array with fWindowLen/2 rows and ncol(Y) columns

  for (iInterval in 1:ncol(Y)) {
      # For each window
      fTmp <- array(NA, c(M_welch, fWindowLen / 2))  # Create a temporary array

      for (m in 1:M_welch) {
          # Divide the window into M_welch sub-windows
          Yint <- Y[((m * fWindowLen) - fWindowLen + 1):(m * fWindowLen), iInterval]

          # Center the sub-window
          Yint <- Yint - mean(Yint)

          # Compute the Fourier transform of the sub-window
          rawfft <- abs(fft(Yint))^2

          # Store the raw Periodogram of the window iInterval in the temporary array
          fTmp[m, ] <- c(rawfft[1], 2 * rawfft[2:(fWindowLen / 2)]) / (fWindowLen^2)
      }

      # Optionally apply log scale to the periodogram
      if (logScale) {
          fTmp <- log(fTmp)
          S[, iInterval] <- apply(fTmp, 2, mean)
      }

      # Calculate the Welch Periodogram by taking the mean along each column
      S[, iInterval] <- apply(fTmp, 2, mean)

      # Adjust the first element if log scale is applied
      if (logScale) {
          S[1, iInterval] <- 0
      }

      # Normalize each column if required
      if (normalize) {
          S[, iInterval] <- S[, iInterval] / sum(S[, iInterval])
      }
  } 
  

    if (bands_analysis){
      
      S_bands <- array(NA, dim = c(5, ncol(S)))

      if  (nrow(S) <= f_sampling/4){stop("Choose a smaller number of Welch regulariser or a larger window length to be able to analyse the bands.")}

      indices <- FreSpeD::get_bands_indices(2*nrow(S), f_sampling)
      

      S_bands[1, ] <- apply(S[indices$delta_indices, ], 2, mean)
      # print(indices$delta_indices)
      S_bands[2, ] <- apply(S[indices$theta_indices, ], 2, mean)
      # print(indices$theta_indices)
      S_bands[3, ] <- apply(S[indices$alpha_indices, ], 2, mean)
      # print(indices$alpha_indices)
      S_bands[4, ] <- apply(S[indices$beta_indices, ], 2, mean)
      # print(indices$beta_indices)
      S_bands[5, ] <- apply(S[indices$gamma_indices, ], 2, mean)
      # print(indices$gamma_indices)

      S <- S_bands
    }

  if (plot) {
    df <- reshape2::melt(S)
    
    # Define custom labels for x-axis
    custom_labels <- seq(round_right(windowLen/2), length(Y), by = windowLen)
    
    # Define custom labels for y-axis
    custom_labels_y <- seq(dim(S)[1])
    
    # Select every 10th label for the x-axis
    selected_indices <- seq(1, length(custom_labels), by = ceiling(length(custom_labels) / 10))
    selected_custom_labels <- custom_labels[selected_indices]
    
    # Select corresponding labels for breaks
    selected_breaks <- seq(dim(S)[2])[selected_indices]
    
    # Plot heatmap with custom x-axis and y-axis labels
    p <- ggplot(df, aes(x = Var2, y = factor(Var1, levels = seq(dim(S)[1]), labels = custom_labels_y), fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black") +  # Define colors
      labs(x = "time", y = "Fundamental frequencies", title = "Periodogram") +  # Add labels
      theme_minimal() +                                    # Choose a theme
      guides(fill = guide_colorbar(title = "Value")) +    # Add color scale legend
      scale_x_continuous(breaks = selected_breaks, labels = selected_custom_labels) # Custom x-axis labels
    
    print(p)
}


  out <- list(S = S, data = Y)
}