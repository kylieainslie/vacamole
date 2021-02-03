#' Calculate force of infection 
#' @param dat data frame of states
#' @param beta transmission parameter
#' @param contact_matrix contact matrix
#' @param N vector of total group sizes
#' @return matrix of force of infection in each age group (columns) at each
#' time point (rows)
#' @keywords vacamole
#' @export
get_foi <- function(dat, beta, contact_matrix, N){

  # sum over different I states for each time step and age group
  I_all <- dat$I + dat$Iv_1d + dat$Iv_2d
  
  # calculate force of infection for each time point
  foi <- t(apply(I_all, 1, function(x){beta * (contact_matrix %*% x/N)}))

  return(foi)
}