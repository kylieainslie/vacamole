#' Calculate transmission parameter beta 
#' @param R0 basic reproduction number
#' @param contact_matrix initial contact matrix
#' @param N vector of total group sizes
#' @param sigma 1/latent period
#' @param gamma 1/infectious period
#' @param contact_matrix2 current contact matrix
#' @param Reff effective reproduction number
#' @param init_s susceptibles used to calculate Reff
#' @param eta 1 - VE after dose 1
#' @param eta2 1- VE after dose 2
#' @return matrix of force of infection in each age group (columns) at each
#' time point (rows)
#' @keywords vacamole
#' @export
get_beta <- function(R0, contact_matrix, N, sigma, gamma, 
                     Reff, contact_matrix2 = NULL, init_s = NULL,
                     eta = NULL, eta2 = NULL, vac_started = FALSE){
  
  n_groups <- length(N)
  Ni <- matrix(rep(N, n_groups),nrow = n_groups)
  Nj <- t(Ni)
  Deff <- contact_matrix * (Ni / Nj)
  
  ones <- rep(1,n_groups)
  F_mat <- matrix(rep(0,(2*n_groups)^2),nrow = 2*n_groups)
  F_mat[1:n_groups,(n_groups+1):dim(F_mat)[2]] <- Deff
  
  v_vec <- c(sigma*ones, gamma*ones);
  V <- diag(v_vec)
  V[(n_groups+1):(2*n_groups),1:n_groups] <- diag(-sigma*ones)
  
  GD <- F_mat %*% solve(V)    # next generation matrix
  d <- as.numeric(eigs(GD,1)$values)
  beta <- R0/d
  
  rtn <- list(beta = beta)
  
  if(!is.null(contact_matrix2) & !is.null(init_s) &
     !is.null(eta) & !is.null(eta2)){
    
    Si <- matrix(rep(init_s, n_groups), nrow = n_groups)
    Deff2 <- contact_matrix2 * (Si / Nj) # effective contact matrix
    
    if(vac_started){
    # create empty matrix for F
    F_mat2 <- matrix(rep(0,(6*n_groups)^2),nrow = 6*n_groups)
    # fill in F matrix with copies of Deff
    # block row 1: E -> I
    F_mat2[1:n_groups,(3*n_groups+1):(4*n_groups)] <- Deff2 # E -> I
    F_mat2[1:n_groups,(4*n_groups+1):(5*n_groups)] <- Deff2 # E -> I
    F_mat2[1:n_groups,(5*n_groups+1):(6*n_groups)] <- Deff2 # E -> I
    # block row 2: Ev_1d -> Iv_1d
    F_mat2[(2*n_groups+1):(3*n_groups),(3*n_groups+1):(4*n_groups)] <- eta * Deff2 # Ev_1d -> Iv_1d
    F_mat2[(2*n_groups+1):(3*n_groups),(4*n_groups+1):(5*n_groups)] <- eta * Deff2 # Ev_1d -> Iv_1d
    F_mat2[(2*n_groups+1):(3*n_groups),(5*n_groups+1):(6*n_groups)] <- eta * Deff2 # Ev_1d -> Iv_1d
    # block row 3: Ev_2d -> Iv_2d
    F_mat2[(3*n_groups+1):(4*n_groups),(3*n_groups+1):(4*n_groups)] <- eta2 * Deff2 # Ev_2d -> Iv_2d
    F_mat2[(3*n_groups+1):(4*n_groups),(4*n_groups+1):(5*n_groups)] <- eta2 * Deff2 # Ev_2d -> Iv_2d
    F_mat2[(3*n_groups+1):(4*n_groups),(5*n_groups+1):(6*n_groups)] <- eta2 * Deff2 # Ev_2d -> Iv_2d
    
    # V matrix
    v_vec2 <- c(-sigma*ones, -sigma*ones, -sigma*ones,
                -gamma*ones, -gamma*ones, -gamma*ones)
    V2 <- diag(v_vec2)
    V2[(3*n_groups+1):(4*n_groups),1:n_groups] <- diag(sigma*ones) # I
    V2[(4*n_groups+1):(5*n_groups), (2*n_groups+1):(3*n_groups)] <- diag(sigma*ones) # Iv_1d
    V2[(5*n_groups+1):(6*n_groups), (3*n_groups+1):(4*n_groups)] <- diag(sigma*ones) # Iv_2d
    
    GD2 <- -F_mat2 %*% solve(V2)
    } else {
      
      F_mat2 <- matrix(rep(0,(2*n_groups)^2),nrow = 2*n_groups)
      F_mat2[1:n_groups,(n_groups+1):dim(F_mat)[2]] <- Deff2
      GD2 <- -F_mat2 %*% solve(V)
    }
    
    d2 <- as.numeric(eigs(GD2,1)$values)
    beta2 <- Reff/d2
    
    rtn <- list(beta = beta,
                beta2 = beta2)
  }
  
  
  return(rtn)
  
}