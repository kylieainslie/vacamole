#' Calculate transmission parameter beta 
#' @param R0 basic reproduction number
#' @param contact_matrix contact matrix
#' @param N vector of total group sizes
#' @param sigma 1/latent period
#' @param gamma 1/infectious period
#' @return matrix of force of infection in each age group (columns) at each
#' time point (rows)
#' @keywords vacamole
#' @export
get_beta <- function(R0, contact_matrix, N, sigma, gamma){
  
  n_groups <- length(N)
  Nj <- matrix(rep(N, n_groups),nrow = n_groups)
  Ni <- t(Nj)
  Deff <- contact_matrix %*% (Ni / Nj)
  
  ones <- rep(1,n_groups)
  F_mat <- matrix(rep(0,(2*n_groups)^2),nrow = 2*n_groups)
  F_mat[1:n_groups,(n_groups+1):dim(F_mat)[2]] <- Deff
  
  v_vec <- c(sigma*ones, gamma*ones);
  V <- diag(v_vec)
  V[(n_groups+1):(2*n_groups),1:n_groups] <- diag(-sigma*ones)
  
  GD <- F_mat %*% solve(V)    # next generation matrix
  d <- as.numeric(eigs(GD,1)$values)
  beta <- R0/d
  
  return(beta)
  
}