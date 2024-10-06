#' @name cusum_stat
#' @title Apply CUSUM statistics
#' @description Applyi the CUSUM statistics to the estimated spectral quantities to get the change points 
#' 
#' @param S numeric array, The estimated spectral quantities for each window (lines) and for the corresponding frequencies (columns)
#' @param thresh numeric, The threshold to decide that a point is a change point, using the statistics eq (5), the 
#'               (default is the theoretical threshold \tau_T Theorem 1 Shorder & Ombao 2019)
#' @param check_neighborhood logical, Within the change-point detection, should the neighbourhood around 
#'  a changepoint candidate be tested to reduce the risk of spurious detections?
#' @param check_par numeric, \Delta_T the diameter of the interval around a change-point used in check_neighborhood
#' @param cp numeric array, containing only seros and has the same dimension as S
#' @param s numeric, The first index of the interval for which FreSpeD is applied
#' @param delta numeric, \delta_T' in the algorithm: the minimal distance between two change points 
#' 
#' @return numric array, the values of UHcusum on the detected change points for each frquencies, and 0 otherwise
#' 
cusum_stat <- function(S, thresh = 0.8 * log(dim(S)[2])^(1.1), check_neighborhood = TRUE, 
                          check_par = floor(log(dim(S)[2])^2.1/5), cp = array(0, dim(S)), 
                          s = 1, delta = 1) {
  D <- dim(S)[1]
  n <- dim(S)[2]
  
  # The creterion on the length of the interval to perform FreSped, (step 2 in Figure 2, Shroder & Ombao)
  if (n > (2 * delta + 1)) {
    XM <- array(NA, dim(S))
    
    # Loop through each row of the input matrix (skip the first raw)
    for (d in 2:D) {
      X <- S[d, ]
      sigM <- mad((X[2:n] - X[1:(n - 1)])/sqrt(2)) # computation of the median absolute deviation of X
      YM <- UHcusum(X)
      YM[1:delta] <- 0 
      YM[(n - delta):(n - 1)] <- 0
      if (sigM != 0){
      XM[d, which((YM/sigM) > thresh)] <- YM[which((YM/sigM) > thresh)]/sigM # fill the matrix XM with the UHcusum values that exceeds the threshold
      } else {XM[d, which((YM/sigM) > thresh)] <- YM[which((YM/sigM) > thresh)]}
    }
    
    # Compute the sum of each column, ignoring NA values (equation (5) Shroder & Ombao)
    thrCUSUM <- colSums(XM, na.rm = TRUE)
    
    # Check for spurious change points in the neighborhood ???... 
    if (check_neighborhood) {
      b <- c()
      cnt <- sum(thrCUSUM > 0)
      i <- 1
      while (i <= cnt) {
        b <- which(thrCUSUM == thrCUSUM[order(thrCUSUM, decreasing = TRUE)[i]])
        if (length(b) == 1) {
          if (any(thrCUSUM[max(b - check_par, 1):min(b + check_par, n)] == 0)) {
            i <- i + 1
          } else {
            i <- cnt + 1
          }
        } else { 
          for (j in 1:length(b)) {
            if (all(thrCUSUM[max(b[j] - check_par, 1):min(b[j] + check_par, n)] > 0)) {
              i <- cnt + 1
            } else if (j == length(b)) {
              i <- i + 1
            }
          }
        }
      }
      exceeded <- !(is.null(b)) 
    } else {
      exceeded <- any(!is.na(XM))
    }
    
    # Update change point matrix if change points are detected
    if (exceeded) {
      if (!check_neighborhood) {
        # b is the index of the maximum element of thrCUSUM
        b <- which(thrCUSUM == max(thrCUSUM))
        if (length(b)>1){
          # if the there are multiple maximums of thrCUSUM pick the smallest index 
          b <- b[1]
        } 
      }
      
      db <- which(!is.na(XM[, b])) 
      
      if (all(db <= nrow(cp))) {
        cp[db, b + s - 1] <- cp[db, b + s - 1] + XM[db, b]
        cp <- cusum_stat(S[, 1:b], thresh, check_neighborhood, check_par, cp, s = s, delta)
        cp <- cusum_stat(S[, (b + 1):n], thresh, check_neighborhood, check_par, cp, s = (s + b), delta)
      } else {
        warning("Index out of bounds, skipping update.") 
      }
    }
  
  }
  
  return(cp)
}