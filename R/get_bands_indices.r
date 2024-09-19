#' @name get_bands_indices
#' @title Get the indices of the frequency bands
#' @description Given the length of the window and the sampling frequency, return the indices of the bands delta, theta, alpha, 
#' beta and gamma.
#' 
#' @param length_window Numeric vector, representing the length of the number of equispaced fundamental frequencies in [0, 0.5].
#' @param f_sampling Numeric, representing the sampling rate frequency.
#' 
#' @return List of vectors, representing the indices of the bands delta, theta, alpha, beta and gamma.
get_bands_indices <- function(length_window, f_sampling){
  m_delta <- 0.5
  M_delta <- 4
  m_theta <- 4 
  M_theta <- 8
  m_alpha <- 8
  M_alpha <- 12
  m_beta <- 12
  M_beta <- 30
  m_gamma <- 30
  M_gamma <- 50

  nF <- floor(length_window/2)

  delta_indices <- c(floor(length_window*m_delta/f_sampling):floor(length_window*M_delta/f_sampling))
  theta_indices <- c(ceiling(length_window*m_theta/f_sampling):floor(length_window*M_theta/f_sampling))
  alpha_indices <- c(ceiling(length_window*m_alpha/f_sampling):floor(length_window*M_alpha/f_sampling))
  beta_indices <- c(ceiling(length_window*m_beta/f_sampling):floor(length_window*M_beta/f_sampling))
  gamma_indices <- c(ceiling(length_window*m_gamma/f_sampling):nF)

  out <- list(delta_indices = delta_indices, theta_indices = theta_indices, alpha_indices = alpha_indices, beta_indices = beta_indices, gamma_indices = gamma_indices)
  
}