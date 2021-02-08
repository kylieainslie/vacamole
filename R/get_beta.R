#' Calculate transmission parameter beta 
#' @param R0 basic reproduction number
#' @param contact_matrix initial contact matrix
#' @param N vector of total group sizes
#' @param sigma 1/latent period
#' @param gamma 1/infectious period
#' @param contact_matrix2 current contact matrix
#' @param Reff effective reproduction number
#' @param init_s susceptibles used to calculate Reff
#' @return matrix of force of infection in each age group (columns) at each
#' time point (rows)
#' @keywords vacamole
#' @export
get_beta <- function(R0, contact_matrix, N, sigma, gamma, 
                     Reff, contact_matrix2 = NULL, init_s = NULL){
  
  n_groups <- length(N)
  Nj <- matrix(rep(N, n_groups),nrow = n_groups)
  Ni <- t(Nj)
  Deff <- contact_matrix * (Ni / Nj)
  
  ones <- rep(1,n_groups)
  F_mat <- matrix(rep(0,(2*n_groups)^2),nrow = 2*n_groups)
  F_mat[1:n_groups,(n_groups+1):dim(F_mat)[2]] <- Deff
  
  v_vec <- c(-sigma*ones, -gamma*ones);
  V <- diag(v_vec)
  V[(n_groups+1):(2*n_groups),1:n_groups] <- diag(sigma*ones)
  
  GD <- -F_mat %*% solve(V)    # next generation matrix
  d <- as.numeric(eigs(GD,1)$values)
  beta <- R0/d
  
  if(!is.null(contact_matrix2) & !is.null(init_s)){
    Sj <- matrix(rep(init_s, n_groups), nrow = n_groups)
    Si <- t(Sj)
    Deff2 <- contact_matrix2 * (Si / Nj)
    F_mat2 <- matrix(rep(0,(6*n_groups)^2),nrow = 6*n_groups)
    F_mat2[1:n_groups,(3*n_groups+1):(4*n_groups)] <- Deff2 # E -> I
    F_mat2[(2*n_groups+1):(3*n_groups),(4*n_groups+1):(5*n_groups)] <- Deff2 # Ev_1d -> Iv_1d
    F_mat2[(3*n_groups+1):(4*n_groups),(5*n_groups+1):(6*n_groups)] <- Deff2 # Ev_2d -> Iv_2d
    
    v_vec2 <- c(-sigma*ones, -sigma*ones, -sigma*ones,
                -gamma*ones, -gamma*ones, -gamma*ones)
    V2 <- diag(v_vec2)
    F_mat2[(3*n_groups+1):(4*n_groups),1:n_groups] <- diag(sigma*ones) # I
    F_mat2[(4*n_groups+1):(5*n_groups), (2*n_groups+1):(3*n_groups)] <- diag(sigma*ones) # Iv_1d
    F_mat2[(5*n_groups+1):(6*n_groups), (3*n_groups+1):(4*n_groups)] <- diag(sigma*ones) # Iv_2d
    
    GD2 <- -F_mat2 %*% solve(V2)
    d2 <- as.numeric(eigs(GD2,1)$values)
    beta2 <- Reff/d2
  }
  
  rtn <- list(beta = beta,
              beta2 = beta2)
  return(rtn)
  
}