#' @name FreSpeD_innerpost
#' @title Check change point
#' @description Check wether a given point b is a true change point or spurious change point 
#' 
#' @param S numeric array, the estimated spectral quantities 
#' @param thresh numeric, numeric, The threshold to decide that a point is a change point
#' @param check_neighborhood logical, logical, Within the change-point detection, should the neighbourhood around 
#'  a changepoint candidate be tested
#' @param b numeric, the candidate to be a change point 
#' @param delta numeric, the minimal distance between two change points 
#' 
#' @return logical, if b is a spurious change point return FALSE, and TRUE otherwise.
#' 
FreSpeD_innerpost <- function(S, thresh, check_neighborhood = TRUE, check_par = floor(log(dim(S)[2])^2.1/5), 
                              b, delta = 1) {
  L <- dim(S)[1]
  n <- dim(S)[2]
  
  # Check if the length of the sequence is suitable for FreSped
  if (n > (2 * delta + 1)) {
    XM <- array(NA, dim(S))
    
    # Iterate through rows (skipping the first one)
    for (l in 2:L) {
      X <- S[l, ]
      sigM <- mad((X[2:n] - X[1:(n - 1)])/sqrt(2))
      YM <- FreSpeD::UHcusum(X)
      YM[1:delta] <- 0
      YM[(n - delta):(n - 1)] <- 0
      XM[l, which((YM/sigM) > thresh)] <- YM[which((YM/sigM) > thresh)]/sigM
    }
    
    # Compute the sum of each column, ignoring NA values
    thrCUSUM <- colSums(XM, na.rm = TRUE)
    
    # Check for spurious change points in the neighborhood
    if (check_neighborhood) {
      if (any(thrCUSUM[max(b - check_par, 1):min(b + check_par, n)] == 0)) 
        return(FALSE)
      else return(TRUE)
    } else {
      if (thrCUSUM[b] != 0) 
        return(TRUE)
      else 
        return(FALSE)
    }
  } else {
    print("warning - cp dist too small in postprocessing") 
  }
}

