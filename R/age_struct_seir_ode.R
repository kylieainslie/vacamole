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
    Shold = c(Shold1, Shold2)
    Sv = c(Sv1, Sv2)
    Shold2 = c(Shold21, Shold22)
    Sv2 = c(Sv21, Sv22)
    E = c(E1, E2)
    Ev = c(Ev1, Ev2)
    Ev2 = c(Ev21, Ev22)
    I = c(I1, I2)
    Iv = c(Iv1, Iv2)
    Iv2 = c(Iv21, Iv22)
    H = c(H1, H2)
    D = c(D1, D2)
    R = c(R1, R2)
    Rv = c(Rv1, Rv2)
    Rv2 = c(Rv21, Rv22)
    ################################################################
    # ODEs:
    lambda <- beta * (C%*%((I + Iv + Iv2)/N))
    
    dS <- -lambda * S - alpha * S/N
    dShold <- alpha * S/N - (1/delay) * Shold - lambda * Shold
    dSv <- (1/delay) * Shold - eta * lambda * Sv - ifelse(Sv>0,alpha2 * Sv/(Sv+Ev+Iv+Rv), alpha2*0)
    dShold2 <- ifelse(Sv>0,alpha2 * Sv/(Sv+Ev+Iv+Rv), alpha2*0) - (1/delay2) * Shold2 - eta * lambda * Shold2
    dSv2 <- (1/delay2) * Shold2 - eta2 * lambda * Sv2
    dE <- lambda * (S + Shold) - sigma * E
    dEv <- eta * lambda * (Sv + Shold2) - sigma * Ev
    dEv2 <- eta * lambda * Sv2 - sigma * Ev2 
    dI <- sigma * E - (gamma + h) * I 
    dIv <- sigma * Ev - (gamma + h) * Iv  
    dIv2 <- sigma * Ev2 - (gamma + h) * Iv2
    dH <- h * (I + Iv + Iv2) - (d + r) * H
    dD <- d * H 
    dR <- gamma * I + r * H 
    dRv <- gamma * Iv + r * H
    dRv2 <- gamma * Iv2 + r * H
    
    ################################################################
    
    dt <- 1
    list(c(dt,dS,dShold,dSv,dShold2,dSv2,dE,dEv,dEv2,
           dI,dIv,dIv2,dH,dD,dR,dRv,dRv2))
  })
}
