#' SEIR ODE model of vaccination with 2 doses and delay to protection
#' @param times vector of times
#' @param init list of inititial states
#' @param params list of parameter values
#' @return List of summary results
#' @keywords vacamole
#' @export
# Define model -----------------------------------------------------
seir_ode_2 <- function(times,init,params){
  with(as.list(c(params,init)), {
    # determine vaccination rate
    if(is.null(tv)){
      alpha <- vac_per_day[t]
    } else {
      if (t>tv && S/N>(1-uptake)){
        alpha <- vac_per_day
      } else {
        alpha <- 0
      }
    }
    
    if(is.null(tv2)){
      alpha2 <- vac_per_day2[t]
    } else{
      if (t>tv2 && S/N>(1-uptake)){ 
        alpha2 <- vac_per_day2
      } else {
        alpha2 <- 0
      }
    }
    
    if(vac_input_perc){
      alpha <- alpha * N
      alpha2 <- alpha2 * (Sv+Ev+Iv+Rv)
    }
    ################################################################
    # ODEs:
    if (constant_foi){
      lambda <- beta*init_inf/N
    } else{lambda <- beta * (I + Iv + Iv2)/N}
    
      dS <- -lambda * S - alpha * S/N
      dShold <- alpha * S/N - (1/delay) * Shold - lambda * Shold
      dSv <- (1/delay) * Shold - eta * lambda * Sv - ifelse(Sv>0,alpha2 * Sv/(Sv+Ev+Iv+Rv), alpha2*0)
      dShold2 <- ifelse(Sv>0,alpha2 * Sv/(Sv+Ev+Iv+Rv), alpha2*0) - (1/delay2) * Shold2 - eta * lambda * Shold2
      dSv2 <- (1/delay2) * Shold2 - eta2 * lambda * Sv2
      dE <- lambda * (S + Shold) - sigma * E
      dEv <- eta * lambda * (Sv + Shold2) - sigma * Ev
      dEv2 <- eta * lambda * Sv2 - sigma * Ev2 
      dI <- sigma * E - (gamma * I) # (gamma + h) * I
      dIv <- sigma * Ev - (gamma * Iv) # (gamma + h) * Iv
      dIv2 <- sigma * Ev2 - (gamma * Iv2) # (gamma + h) * Iv2
      # dH <- h * (I + Iv + Iv2) - (d + r) * H
      # dD <- d * H 
      dR <- gamma * I #+ r * H 
      dRv <- gamma * Iv #+ r * H
      dRv2 <- gamma * Iv2 #+ r * H
      
    ################################################################
    
    dt <- 1
    list(c(dt,dS,dShold,dSv,dShold2,dSv2,dE,dEv,dEv2,
           dI,dIv,dIv2,dR,dRv,dRv2))
  })
}
