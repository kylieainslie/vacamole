#' Calculate force of infection 
#' @param times vector of times
#' @param params list of parameter values
#' @return list of percentage of each group to be vaccinated at time point
#' and the composite ve for each group (based on how much of each vaccine is
#' used and the ve of each vaccine)
#' @keywords vacamole
#' @export
get_lambda <- function(params, I, I_vd, Iv_2d,
                       S, Shold_1d, Sv_1d, Shold_2d, Sv_2d,
                       eta, eta2){
  with(as.list(c(params,I, I_vd, Iv_2d,
                 S, Shold_1d, Sv_1d, Shold_2d, Sv_2d,
                 eta, eta2)),{
    
    l <- beta * (C %*% ((I + Iv_1d + Iv_2d)/N))
    
    if (constant_contact_matrix == FALSE){
      incidence <- (S + Shold_1d + (eta * (Sv_1d + Shold_2d)) + (eta2 * Sv_2d)) * l
      if(sum(incidence) < 7/100000) {
        flag <- 0
      } else if (sum(incidence) > 21/100000){
        flag <- 1
      }
      C <- ifelse(flag == 1, c_lockdown, c_relaxed)
      l <- beta * (C %*% ((I + Iv_1d + Iv_2d)/N))
    }
    
    list(l)
    
  }) 
}