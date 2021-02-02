#' Calculate force of infection 
#' @param I data frame of infectious (unvaccinated)
#' @param Iv_1d data frame of infectious (vaccinated with dose 1)
#' @param Iv_2d data frame of infectious (vaccinated with dose 2)
#' @param beta transmission parameter
#' @param contact_matrix contact matrix
#' @param N vector of total group sizes
#' @return List of summary results
#' @keywords vacamole
#' @export
get_foi <- function(I, Iv_1d, Iv_2d, beta, contact_matrix, N){
  foi <- 0
  #beta * C %*% rowSums(seir_out[,c("I1", "Iv1", "Iv21")])/N[1]

  return(foi)
}