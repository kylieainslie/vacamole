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
    
    # define initial state vectors from input ----------------------
    S = c(S1, S2, S3, S4, S5, S6, S7, S8, S9)
    Shold_1d = c(Shold_1d1, Shold_1d2, Shold_1d3, Shold_1d4, Shold_1d5, Shold_1d6, 
                 Shold_1d7, Shold_1d8, Shold_1d9)
    Sv_1d = c(Sv_1d1, Sv_1d2, Sv_1d3, Sv_1d4, Sv_1d5, Sv_1d6, Sv_1d7, Sv_1d8, Sv_1d9)
    Shold_2d = c(Shold_2d1, Shold_2d2, Shold_2d3, Shold_2d4, Shold_2d5, Shold_2d6, 
                 Shold_2d7, Shold_2d8, Shold_2d9)
    Sv_2d = c(Sv_2d1, Sv_2d2, Sv_2d3, Sv_2d4, Sv_2d5, Sv_2d6, Sv_2d7, Sv_2d8, Sv_2d9)
    E = c(E1, E2, E3, E4, E5, E6, E7, E8, E9)
    Ev_1d = c(Ev_1d1, Ev_1d2, Ev_1d3, Ev_1d4, Ev_1d5, Ev_1d6, Ev_1d7, Ev_1d8, Ev_1d9)
    Ev_2d = c(Ev_2d1, Ev_2d2, Ev_2d3, Ev_2d4, Ev_2d5, Ev_2d6, Ev_2d7, Ev_2d8, Ev_2d9)
    I = c(I1, I2, I3, I4, I5, I6, I7, I8, I9)
    Iv_1d = c(Iv_1d1, Iv_1d2, Iv_1d3, Iv_1d4, Iv_1d5, Iv_1d6, Iv_1d7, Iv_1d8, Iv_1d9)
    Iv_2d = c(Iv_2d1, Iv_2d2, Iv_2d3, Iv_2d4, Iv_2d5, Iv_2d6, Iv_2d7, Iv_2d8, Iv_2d9)
    H = c(H1, H2, H3, H4, H5, H6, H7, H8, H9)
    Hv_1d = c(Hv_1d1, Hv_1d2, Hv_1d3, Hv_1d4, Hv_1d5, Hv_1d6, Hv_1d7, Hv_1d8, Hv_1d9)
    Hv_2d = c(Hv_2d1, Hv_2d2, Hv_2d3, Hv_2d4, Hv_2d5, Hv_2d6, Hv_2d7, Hv_2d8, Hv_2d9)
    D = c(D1, D2, D3, D4, D5, D6, D7, D8, D9)
    R = c(R1, R2, R3, R4, R5, R6, R7, R8, R9)
    Rv_1d = c(Rv_1d1, Rv_1d2, Rv_1d3, Rv_1d4, Rv_1d5, Rv_1d6, Rv_1d7, Rv_1d8, Rv_1d9)
    Rv_2d = c(Rv_2d1, Rv_2d2, Rv_2d3, Rv_2d4, Rv_2d5, Rv_2d6, Rv_2d7, Rv_2d8, Rv_2d9)
    
    # determine vaccination rate -----------------------------------
    alpha <- ifelse(t>tv && S/N > (1 - uptake), vac_per_day, 0)
    alpha2 <- ifelse(t>tv2 && S/N > (1 - uptake), vac_per_day2, 0)
    
    ################################################################
    # ODEs:
    lambda <- beta * (C %*% ((I + Iv_1d + Iv_2d)/N))
    
    dS <- -lambda * S - alpha * (S/N)
    dShold_1d <- alpha * (S/N) - (1/delay) * Shold_1d - lambda * Shold_1d
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
