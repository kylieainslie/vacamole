#' Age-structured SEIR ODE model of vaccination with 2 doses and delay to protection
#' @param times vector of times
#' @param init list of inititial states
#' @param params list of parameter values
#' @return List of summary results
#' @keywords vacamole
#' @export
# Define model -----------------------------------------------------
age_struct_seir_ode <- function(times,init,params){
  with(as.list(c(params,init)), {
    # determine vaccination rate
    if (t>tv && S/N>(1-uptake)){
      alpha <- vac_per_day
    } else {
      alpha <- 0
    }
    
    if (t>tv2 && S/N>(1-uptake)){ 
      alpha2 <- vac_per_day2
    } else {
      alpha2 <- 0
    }
    
    # define initial state vectors from input ----------------------
    S = c(S1, S2)
    Shold_1d = c(Shold_1d1, Shold_1d2)
    Sv_1d = c(Sv_1d1, Sv_1d2)
    Shold_2d = c(Shold_2d1, Shold_2d2)
    Sv_2d = c(Sv_2d1, Sv_2d2)
    E = c(E1, E2)
    Ev_1d = c(Ev_1d1, Ev_1d2)
    Ev_2d = c(Ev_2d1, Ev_2d2)
    I = c(I1, I2)
    Iv_1d = c(Iv_1d1, Iv_1d2)
    Iv_2d = c(Iv_2d1, Iv_2d2)
    H = c(H1, H2)
    Hv_1d = c(Hv_1d1, Hv_1d2)
    Hv_2d = c(Hv_2d1, Hv_2d2)
    D = c(D1, D2)
    R = c(R1, R2)
    Rv_1d = c(Rv_1d1, Rv_1d2)
    Rv_2d = c(Rv_2d1, Rv_2d2)
    ################################################################
    # ODEs:
    lambda <- beta * (C%*%((I + Iv_1d + Iv_2d)/N))
    
    dS <- -lambda * S - alpha * S/N
    dShold_1d <- alpha * S/N - (1/delay) * Shold_1d - lambda * Shold_1d
    dSv_1d <- (1/delay) * Shold_1d - eta * lambda * Sv_1d - ifelse(Sv_1d>0,alpha2 * Sv_1d/(Sv_1d+Ev_1d+Iv_1d+Rv_1d), alpha2*0)
    dShold_2d <- ifelse(Sv_1d>0,alpha2 * Sv_1d/(Sv_1d+Ev_1d+Iv_1d+Rv_1d), alpha2*0) - (1/delay2) * Shold_2d - eta * lambda * Shold_2d
    dSv_2d <- (1/delay2) * Shold_2d - eta2 * lambda * Sv_2d
    dE <- lambda * (S + Shold_1d) - sigma * E
    dEv_1d <- eta * lambda * (Sv_1d + Shold_2d) - sigma * Ev_1d
    dEv_2d <- eta * lambda * Sv_2d - sigma * Ev_2d 
    dI <- sigma * E - (gamma + h) * I 
    dIv_1d <- sigma * Ev_1d - (gamma + h) * Iv_1d  
    dIv_2d <- sigma * Ev_2d - (gamma + h) * Iv_2d
    dH <- h * I - (d + r) * H
    dHv_1d <- h * Iv_1d - (d + r) * Hv_1d
    dHv_2d <- h * Iv_2d - (d + r) * Hv_2d
    dD <- d * (H + Hv_1d + Hv_2d) 
    dR <- gamma * I + r * H 
    dRv_1d <- gamma * Iv_1d + r * Hv_1d
    dRv_2d <- gamma * Iv_2d + r * Hv_2d
    
    ################################################################
    
    dt <- 1
    list(c(dt,dS,dShold_1d,dSv_1d,dShold_2d,dSv_2d,dE,dEv_1d,dEv_2d,
           dI,dIv_1d,dIv_2d,dH, dHv_1d,dHv_2d,dD,dR,dRv_1d,dRv_2d))
  })
}
