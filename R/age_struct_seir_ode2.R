#' Age-structured SEIR ODE model of vaccination with 2 dose primary series with 3 booster doses (3rd, 4th, and 5th doses)
#' @param times vector of times
#' @param init list of initial states
#' @param params list of parameter values
#' @return List of summary results
#' @keywords vacamole
#' @export
# Define model -----------------------------------------------------
age_struct_seir_ode2 <- function(times, init, params) {
  with(as.list(c(params, init)), {
    #print(t)
    # define initial state vectors from input ----------------------
    # susceptible
    S <- c(S1, S2, S3, S4, S5, S6, S7, S8, S9)     
    Shold_1d <- c(
      Shold_1d1, Shold_1d2, Shold_1d3, Shold_1d4, Shold_1d5, Shold_1d6,
      Shold_1d7, Shold_1d8, Shold_1d9
    )
    Sv_1d <- c(Sv_1d1, Sv_1d2, Sv_1d3, Sv_1d4, Sv_1d5, Sv_1d6, Sv_1d7, Sv_1d8, Sv_1d9)
    Shold_2d <- c(
      Shold_2d1, Shold_2d2, Shold_2d3, Shold_2d4, Shold_2d5, Shold_2d6,
      Shold_2d7, Shold_2d8, Shold_2d9
    )
    Sv_2d <- c(Sv_2d1, Sv_2d2, Sv_2d3, Sv_2d4, Sv_2d5, Sv_2d6, Sv_2d7, Sv_2d8, Sv_2d9)
    Shold_3d <- c(
      Shold_3d1, Shold_3d2, Shold_3d3, Shold_3d4, Shold_3d5, Shold_3d6,
      Shold_3d7, Shold_3d8, Shold_3d9
    )
    Sv_3d <- c(Sv_3d1, Sv_3d2, Sv_3d3, Sv_3d4, Sv_3d5, Sv_3d6, Sv_3d7, Sv_3d8, Sv_3d9)
    Shold_4d <- c(
      Shold_4d1, Shold_4d2, Shold_4d3, Shold_4d4, Shold_4d5, Shold_4d6,
      Shold_4d7, Shold_4d8, Shold_4d9
    )
    Sv_4d <- c(Sv_4d1, Sv_4d2, Sv_4d3, Sv_4d4, Sv_4d5, Sv_4d6, Sv_4d7, Sv_4d8, Sv_4d9)
    Shold_5d <- c(
      Shold_5d1, Shold_5d2, Shold_5d3, Shold_5d4, Shold_5d5, Shold_5d6,
      Shold_5d7, Shold_5d8, Shold_5d9
    )
    Sv_5d <- c(Sv_5d1, Sv_5d2, Sv_5d3, Sv_5d4, Sv_5d5, Sv_5d6, Sv_5d7, Sv_5d8, Sv_5d9)
    # exposed
    E <- c(E1, E2, E3, E4, E5, E6, E7, E8, E9)     
    Ev_1d <- c(Ev_1d1, Ev_1d2, Ev_1d3, Ev_1d4, Ev_1d5, Ev_1d6, Ev_1d7, Ev_1d8, Ev_1d9)
    Ev_2d <- c(Ev_2d1, Ev_2d2, Ev_2d3, Ev_2d4, Ev_2d5, Ev_2d6, Ev_2d7, Ev_2d8, Ev_2d9)
    Ev_3d <- c(Ev_3d1, Ev_3d2, Ev_3d3, Ev_3d4, Ev_3d5, Ev_3d6, Ev_3d7, Ev_3d8, Ev_3d9)
    Ev_4d <- c(Ev_4d1, Ev_4d2, Ev_4d3, Ev_4d4, Ev_4d5, Ev_4d6, Ev_4d7, Ev_4d8, Ev_4d9)
    Ev_5d <- c(Ev_5d1, Ev_5d2, Ev_5d3, Ev_5d4, Ev_5d5, Ev_5d6, Ev_5d7, Ev_5d8, Ev_5d9)
    # infectious
    I <- c(I1, I2, I3, I4, I5, I6, I7, I8, I9)     
    Iv_1d <- c(Iv_1d1, Iv_1d2, Iv_1d3, Iv_1d4, Iv_1d5, Iv_1d6, Iv_1d7, Iv_1d8, Iv_1d9)
    Iv_2d <- c(Iv_2d1, Iv_2d2, Iv_2d3, Iv_2d4, Iv_2d5, Iv_2d6, Iv_2d7, Iv_2d8, Iv_2d9)
    Iv_3d <- c(Iv_3d1, Iv_3d2, Iv_3d3, Iv_3d4, Iv_3d5, Iv_3d6, Iv_3d7, Iv_3d8, Iv_3d9)
    Iv_4d <- c(Iv_4d1, Iv_4d2, Iv_4d3, Iv_4d4, Iv_4d5, Iv_4d6, Iv_4d7, Iv_4d8, Iv_4d9)
    Iv_5d <- c(Iv_5d1, Iv_5d2, Iv_5d3, Iv_5d4, Iv_5d5, Iv_5d6, Iv_5d7, Iv_5d8, Iv_5d9)
    # hospitalized
    H <- c(H1, H2, H3, H4, H5, H6, H7, H8, H9)     
    Hv_1d <- c(Hv_1d1, Hv_1d2, Hv_1d3, Hv_1d4, Hv_1d5, Hv_1d6, Hv_1d7, Hv_1d8, Hv_1d9)
    Hv_2d <- c(Hv_2d1, Hv_2d2, Hv_2d3, Hv_2d4, Hv_2d5, Hv_2d6, Hv_2d7, Hv_2d8, Hv_2d9)
    Hv_3d <- c(Hv_3d1, Hv_3d2, Hv_3d3, Hv_3d4, Hv_3d5, Hv_3d6, Hv_3d7, Hv_3d8, Hv_3d9)
    Hv_4d <- c(Hv_4d1, Hv_4d2, Hv_4d3, Hv_4d4, Hv_4d5, Hv_4d6, Hv_4d7, Hv_4d8, Hv_4d9)
    Hv_5d <- c(Hv_5d1, Hv_5d2, Hv_5d3, Hv_5d4, Hv_5d5, Hv_5d6, Hv_5d7, Hv_5d8, Hv_5d9)
    # ICU
    IC <- c(IC1, IC2, IC3, IC4, IC5, IC6, IC7, IC8, IC9)
    ICv_1d <- c(ICv_1d1, ICv_1d2, ICv_1d3, ICv_1d4, ICv_1d5, ICv_1d6, ICv_1d7, ICv_1d8, ICv_1d9)
    ICv_2d <- c(ICv_2d1, ICv_2d2, ICv_2d3, ICv_2d4, ICv_2d5, ICv_2d6, ICv_2d7, ICv_2d8, ICv_2d9)
    ICv_3d <- c(ICv_3d1, ICv_3d2, ICv_3d3, ICv_3d4, ICv_3d5, ICv_3d6, ICv_3d7, ICv_3d8, ICv_3d9)
    ICv_4d <- c(ICv_4d1, ICv_4d2, ICv_4d3, ICv_4d4, ICv_4d5, ICv_4d6, ICv_4d7, ICv_4d8, ICv_4d9)
    ICv_5d <- c(ICv_5d1, ICv_5d2, ICv_5d3, ICv_5d4, ICv_5d5, ICv_5d6, ICv_5d7, ICv_5d8, ICv_5d9)
    # return to hospital ward after ICU
    H_IC <- c(H_IC1, H_IC2, H_IC3, H_IC4, H_IC5, H_IC6, H_IC7, H_IC8, H_IC9)
    H_ICv_1d <- c(H_ICv_1d1, H_ICv_1d2, H_ICv_1d3, H_ICv_1d4, H_ICv_1d5, H_ICv_1d6, H_ICv_1d7, H_ICv_1d8, H_ICv_1d9)
    H_ICv_2d <- c(H_ICv_2d1, H_ICv_2d2, H_ICv_2d3, H_ICv_2d4, H_ICv_2d5, H_ICv_2d6, H_ICv_2d7, H_ICv_2d8, H_ICv_2d9)
    H_ICv_3d <- c(H_ICv_3d1, H_ICv_3d2, H_ICv_3d3, H_ICv_3d4, H_ICv_3d5, H_ICv_3d6, H_ICv_3d7, H_ICv_3d8, H_ICv_3d9)
    H_ICv_4d <- c(H_ICv_4d1, H_ICv_4d2, H_ICv_4d3, H_ICv_4d4, H_ICv_4d5, H_ICv_4d6, H_ICv_4d7, H_ICv_4d8, H_ICv_4d9)
    H_ICv_5d <- c(H_ICv_5d1, H_ICv_5d2, H_ICv_5d3, H_ICv_5d4, H_ICv_5d5, H_ICv_5d6, H_ICv_5d7, H_ICv_5d8, H_ICv_5d9)
    # death
    D <- c(D1, D2, D3, D4, D5, D6, D7, D8, D9) 
    # recovered
    R <- c(R1, R2, R3, R4, R5, R6, R7, R8, R9)     
    Rv_1d <- c(Rv_1d1, Rv_1d2, Rv_1d3, Rv_1d4, Rv_1d5, Rv_1d6, Rv_1d7, Rv_1d8, Rv_1d9)
    Rv_2d <- c(Rv_2d1, Rv_2d2, Rv_2d3, Rv_2d4, Rv_2d5, Rv_2d6, Rv_2d7, Rv_2d8, Rv_2d9)
    Rv_3d <- c(Rv_3d1, Rv_3d2, Rv_3d3, Rv_3d4, Rv_3d5, Rv_3d6, Rv_3d7, Rv_3d8, Rv_3d9)
    Rv_4d <- c(Rv_4d1, Rv_4d2, Rv_4d3, Rv_4d4, Rv_4d5, Rv_4d6, Rv_4d7, Rv_4d8, Rv_4d9)
    Rv_5d <- c(Rv_5d1, Rv_5d2, Rv_5d3, Rv_5d4, Rv_5d5, Rv_5d6, Rv_5d7, Rv_5d8, Rv_5d9)
    # extra recovered compartments so that waning immunity is not exponential
    R_1w <- c(R_1w1, R_1w2, R_1w3, R_1w4, R_1w5, R_1w6, R_1w7, R_1w8, R_1w9)
    Rv_1d_1w <- c(Rv_1d_1w1, Rv_1d_1w2, Rv_1d_1w3, Rv_1d_1w4, Rv_1d_1w5, Rv_1d_1w6, Rv_1d_1w7, Rv_1d_1w8, Rv_1d_1w9)
    Rv_2d_1w <- c(Rv_2d_1w1, Rv_2d_1w2, Rv_2d_1w3, Rv_2d_1w4, Rv_2d_1w5, Rv_2d_1w6, Rv_2d_1w7, Rv_2d_1w8, Rv_2d_1w9)
    Rv_3d_1w <- c(Rv_3d_1w1, Rv_3d_1w2, Rv_3d_1w3, Rv_3d_1w4, Rv_3d_1w5, Rv_3d_1w6, Rv_3d_1w7, Rv_3d_1w8, Rv_3d_1w9)
    Rv_4d_1w <- c(Rv_4d_1w1, Rv_4d_1w2, Rv_4d_1w3, Rv_4d_1w4, Rv_4d_1w5, Rv_4d_1w6, Rv_4d_1w7, Rv_4d_1w8, Rv_4d_1w9)
    Rv_5d_1w <- c(Rv_5d_1w1, Rv_5d_1w2, Rv_5d_1w3, Rv_5d_1w4, Rv_5d_1w5, Rv_5d_1w6, Rv_5d_1w7, Rv_5d_1w8, Rv_5d_1w9)
    
    R_2w <- c(R_2w1, R_2w2, R_2w3, R_2w4, R_2w5, R_2w6, R_2w7, R_2w8, R_2w9)
    Rv_1d_2w <- c(Rv_1d_2w1, Rv_1d_2w2, Rv_1d_2w3, Rv_1d_2w4, Rv_1d_2w5, Rv_1d_2w6, Rv_1d_2w7, Rv_1d_2w8, Rv_1d_2w9)
    Rv_2d_2w <- c(Rv_2d_2w1, Rv_2d_2w2, Rv_2d_2w3, Rv_2d_2w4, Rv_2d_2w5, Rv_2d_2w6, Rv_2d_2w7, Rv_2d_2w8, Rv_2d_2w9)
    Rv_3d_2w <- c(Rv_3d_2w1, Rv_3d_2w2, Rv_3d_2w3, Rv_3d_2w4, Rv_3d_2w5, Rv_3d_2w6, Rv_3d_2w7, Rv_3d_2w8, Rv_3d_2w9)
    Rv_4d_2w <- c(Rv_4d_2w1, Rv_4d_2w2, Rv_4d_2w3, Rv_4d_2w4, Rv_4d_2w5, Rv_4d_2w6, Rv_4d_2w7, Rv_4d_2w8, Rv_4d_2w9)
    Rv_5d_2w <- c(Rv_5d_2w1, Rv_5d_2w2, Rv_5d_2w3, Rv_5d_2w4, Rv_5d_2w5, Rv_5d_2w6, Rv_5d_2w7, Rv_5d_2w8, Rv_5d_2w9)
    
    R_3w <- c(R_3w1, R_3w2, R_3w3, R_3w4, R_3w5, R_3w6, R_3w7, R_3w8, R_3w9)
    Rv_1d_3w <- c(Rv_1d_3w1, Rv_1d_3w2, Rv_1d_3w3, Rv_1d_3w4, Rv_1d_3w5, Rv_1d_3w6, Rv_1d_3w7, Rv_1d_3w8, Rv_1d_3w9)
    Rv_2d_3w <- c(Rv_2d_3w1, Rv_2d_3w2, Rv_2d_3w3, Rv_2d_3w4, Rv_2d_3w5, Rv_2d_3w6, Rv_2d_3w7, Rv_2d_3w8, Rv_2d_3w9)
    Rv_3d_3w <- c(Rv_3d_3w1, Rv_3d_3w2, Rv_3d_3w3, Rv_3d_3w4, Rv_3d_3w5, Rv_3d_3w6, Rv_3d_3w7, Rv_3d_3w8, Rv_3d_3w9)
    Rv_4d_3w <- c(Rv_4d_3w1, Rv_4d_3w2, Rv_4d_3w3, Rv_4d_3w4, Rv_4d_3w5, Rv_4d_3w6, Rv_4d_3w7, Rv_4d_3w8, Rv_4d_3w9)
    Rv_5d_3w <- c(Rv_5d_3w1, Rv_5d_3w2, Rv_5d_3w3, Rv_5d_3w4, Rv_5d_3w5, Rv_5d_3w6, Rv_5d_3w7, Rv_5d_3w8, Rv_5d_3w9)
    
    # define vaccination parameters ---------------------------------
    # don't index parameters when there's no vaccination, it's faster!
    index <- floor(times) + 1              # use floor of time point + 1 to index df
    # daily vac rate
    alpha1 <- params$alpha1[index, -1] # remove date column
    alpha2 <- params$alpha2[index, -1]
    alpha3 <- params$alpha3[index, -1]
    alpha4 <- params$alpha4[index, -1]
    alpha5 <- params$alpha5[index, -1]
    
    # delay to protection
    delay1 <- params$delay1[index, -1]
    delay2 <- params$delay2[index, -1]
    delay3 <- params$delay3[index, -1]
    delay4 <- params$delay4[index, -1]
    delay5 <- params$delay5[index, -1]
    
    # protection against infection (1 - VE_inf)
    eta1   <- params$eta1[index, -1]
    eta2   <- params$eta2[index, -1]
    eta3   <- params$eta3[index, -1]
    eta4   <- params$eta4[index, -1]
    eta5   <- params$eta5[index, -1]
    
    # protection against hospitalisation 1 - (1 - VE_hosp) / (1 - VE_inf)
    eta_hosp1   <- params$eta_hosp1[index, -1]
    eta_hosp2   <- params$eta_hosp2[index, -1]
    eta_hosp3   <- params$eta_hosp3[index, -1]
    eta_hosp4   <- params$eta_hosp4[index, -1]
    eta_hosp5   <- params$eta_hosp5[index, -1]
    
    # protection against transmission (1 - VE_trans)
    eta_trans1   <- as.numeric(params$eta_trans1[index, -1])
    eta_trans2   <- as.numeric(params$eta_trans2[index, -1])
    eta_trans3   <- as.numeric(params$eta_trans3[index, -1])
    eta_trans4   <- as.numeric(params$eta_trans4[index, -1])
    eta_trans5   <- as.numeric(params$eta_trans5[index, -1])
    # ---------------------------------------------------------------
    
    # determine force of infection ----------------------------------
    # seasonality 
    calendar_day <- lubridate::yday(as.Date(times, origin = calendar_start_date))
    # emergence of new variants
    var1 <- (sigma1 + sigma2*tanh((t-t_var1)/l))
    var2 <- (sigma1 + sigma2*tanh((t-t_var2)/l))
    # determine transmission rate with seasonal/variant effects
    beta_t <- beta * (1 + beta1 * cos(2 * pi * calendar_day / 365.24)) * 
      (var_emerg) * (1 + var1 + var2) + (!var_emerg) * 1
    # calculate force of infection (lambda)
    lambda <- beta_t * (contact_mat %*% (I + (eta_trans1 * Iv_1d) + (eta_trans2 * Iv_2d) + (eta_trans3 * Iv_3d)
                                         + (eta_trans4 * Iv_4d) + (eta_trans5 * Iv_5d)
    )) 
    # ---------------------------------------------------------------
    
    #################################################################
    # ODEs:
    dS <- -lambda * S - alpha1 * S + (omega*4) * R_3w
    dShold_1d <- alpha1 * S - (1/delay1) * Shold_1d - lambda * Shold_1d
    dSv_1d <- (1/delay1) * Shold_1d - eta1 * lambda * Sv_1d - alpha2 * Sv_1d + (omega*4) * Rv_1d_3w
    dShold_2d <- alpha2 * Sv_1d - (1/delay2) * Shold_2d - eta1 * lambda * Shold_2d
    dSv_2d <- (1/delay2) * Shold_2d - eta2 * lambda * Sv_2d - alpha3 * Sv_2d + (omega*4) * Rv_2d_3w
    dShold_3d <- alpha3 * Sv_2d - (1/delay3) * Shold_3d - eta2 * lambda * Shold_3d
    dSv_3d <- (1/delay3) * Shold_3d - eta3 * lambda * Sv_3d - alpha4 * Sv_3d + (omega*4) * Rv_3d_3w     
    dShold_4d <- alpha4 * Sv_3d - (1/delay4) * Shold_4d - eta3 * lambda * Shold_4d
    dSv_4d <- (1/delay4) * Shold_4d - eta4 * lambda * Sv_4d - alpha5 * Sv_4d + (omega*4) * Rv_4d_3w  
    dShold_5d <- alpha5 * Sv_4d - (1/delay5) * Shold_5d - eta4 * lambda * Shold_5d
    dSv_5d <- (1/delay5) * Shold_5d - eta5 * lambda * Sv_5d + (omega*4) * Rv_5d_3w
    
    dE     <- lambda * S + lambda * Shold_1d - sigma * E + epsilon
    dEv_1d <- eta1 * lambda * Sv_1d + eta1 * lambda * Shold_2d - sigma * Ev_1d
    dEv_2d <- eta2 * lambda * Sv_2d + eta2 * lambda * Shold_3d - sigma * Ev_2d
    dEv_3d <- eta3 * lambda * Sv_3d + eta3 * lambda * Shold_4d - sigma * Ev_3d    
    dEv_4d <- eta4 * lambda * Sv_4d + eta4 * lambda * Shold_5d - sigma * Ev_4d  
    dEv_5d <- eta5 * lambda * Sv_5d - sigma * Ev_5d
    
    dI     <- sigma * E - (gamma + h) * I
    dIv_1d <- sigma * Ev_1d - (gamma + eta_hosp1 * h) * Iv_1d
    dIv_2d <- sigma * Ev_2d - (gamma + eta_hosp2 * h) * Iv_2d
    dIv_3d <- sigma * Ev_3d - (gamma + eta_hosp3 * h) * Iv_3d
    dIv_4d <- sigma * Ev_4d - (gamma + eta_hosp4 * h) * Iv_4d
    dIv_5d <- sigma * Ev_5d - (gamma + eta_hosp5 * h) * Iv_5d
    
    dH     <- h * I - (i1 + d + r) * H
    dHv_1d <- eta_hosp1 * h * Iv_1d - (i1 + d + r) * Hv_1d
    dHv_2d <- eta_hosp2 * h * Iv_2d - (i1 + d + r) * Hv_2d
    dHv_3d <- eta_hosp3 * h * Iv_3d - (i1 + d + r) * Hv_3d
    dHv_4d <- eta_hosp4 * h * Iv_4d - (i1 + d + r) * Hv_4d
    dHv_5d <- eta_hosp5 * h * Iv_5d - (i1 + d + r) * Hv_5d
    
    dIC     <- i1 * H - (i2 + d_ic) * IC
    dICv_1d <- i1 * Hv_1d - (i2 + d_ic) * ICv_1d
    dICv_2d <- i1 * Hv_2d - (i2 + d_ic) * ICv_2d
    dICv_3d <- i1 * Hv_3d - (i2 + d_ic) * ICv_3d
    dICv_4d <- i1 * Hv_4d - (i2 + d_ic) * ICv_4d
    dICv_5d <- i1 * Hv_5d - (i2 + d_ic) * ICv_5d
    
    dH_IC     <- i2 * IC - (r_ic + d_hic) * H_IC
    dH_ICv_1d <- i2 * ICv_1d - (r_ic + d_hic) * H_ICv_1d
    dH_ICv_2d <- i2 * ICv_2d - (r_ic + d_hic) * H_ICv_2d
    dH_ICv_3d <- i2 * ICv_3d - (r_ic + d_hic) * H_ICv_3d
    dH_ICv_4d <- i2 * ICv_4d - (r_ic + d_hic) * H_ICv_4d
    dH_ICv_5d <- i2 * ICv_5d - (r_ic + d_hic) * H_ICv_5d
    
    dD <- d * (H + Hv_1d + Hv_2d + Hv_3d + Hv_4d + Hv_5d) +               # 
      d_ic * (IC + ICv_1d + ICv_2d + ICv_3d + ICv_4d + ICv_5d) +           # 
      d_hic * (H_IC + H_ICv_1d + H_ICv_2d + H_ICv_3d + H_ICv_4d + H_ICv_5d) # 
    
    dR     <- (gamma * I) + (r * H) + (r_ic * H_IC) - ((omega*4) * R)
    dRv_1d <- (gamma * Iv_1d) + (r * Hv_1d) + (r_ic * H_ICv_1d) - ((omega*4) * Rv_1d)
    dRv_2d <- (gamma * Iv_2d) + (r * Hv_2d) + (r_ic * H_ICv_2d) - ((omega*4) * Rv_2d)
    dRv_3d <- (gamma * Iv_3d) + (r * Hv_3d) + (r_ic * H_ICv_3d) - ((omega*4) * Rv_3d)
    dRv_4d <- (gamma * Iv_4d) + (r * Hv_4d) + (r_ic * H_ICv_4d) - ((omega*4) * Rv_4d)
    dRv_5d <- (gamma * Iv_5d) + (r * Hv_5d) + (r_ic * H_ICv_5d) - ((omega*4) * Rv_5d)
    
    dR_1w     <- (omega*4) * R - (omega*4) * R_1w
    dRv_1d_1w <- (omega*4) * Rv_1d - (omega*4) * Rv_1d_1w
    dRv_2d_1w <- (omega*4) * Rv_2d - (omega*4) * Rv_2d_1w
    dRv_3d_1w <- (omega*4) * Rv_3d - (omega*4) * Rv_3d_1w
    dRv_4d_1w <- (omega*4) * Rv_4d - (omega*4) * Rv_4d_1w
    dRv_5d_1w <- (omega*4) * Rv_5d - (omega*4) * Rv_5d_1w
    
    dR_2w     <- (omega*4) * R_1w - (omega*4) * R_2w
    dRv_1d_2w <- (omega*4) * Rv_1d_1w - (omega*4) * Rv_1d_2w
    dRv_2d_2w <- (omega*4) * Rv_2d_1w - (omega*4) * Rv_2d_2w
    dRv_3d_2w <- (omega*4) * Rv_3d_1w - (omega*4) * Rv_3d_2w
    dRv_4d_2w <- (omega*4) * Rv_4d_1w - (omega*4) * Rv_4d_2w
    dRv_5d_2w <- (omega*4) * Rv_5d_1w - (omega*4) * Rv_5d_2w
    
    dR_3w     <- (omega*4) * R_2w - (omega*4) * R_3w
    dRv_1d_3w <- (omega*4) * Rv_1d_2w - (omega*4) * Rv_1d_3w
    dRv_2d_3w <- (omega*4) * Rv_2d_2w - (omega*4) * Rv_2d_3w
    dRv_3d_3w <- (omega*4) * Rv_3d_2w - (omega*4) * Rv_3d_3w
    dRv_4d_3w <- (omega*4) * Rv_4d_2w - (omega*4) * Rv_4d_3w
    dRv_5d_3w <- (omega*4) * Rv_5d_2w - (omega*4) * Rv_5d_3w
    #################################################################
    dt <- 1
    # output --------------------------------------------------------
    list(c(dt, dS, dSv_1d, dSv_2d,dSv_3d, dSv_4d, dSv_5d, 
           dShold_1d, dShold_2d,  dShold_3d, dShold_4d, dShold_5d, 
           dE, dEv_1d, dEv_2d, dEv_3d, dEv_4d, dEv_5d,
           dI, dIv_1d, dIv_2d, dIv_3d, dIv_4d, dIv_5d,
           dH, dHv_1d, dHv_2d, dHv_3d, dHv_4d, dHv_5d,
           dIC, dICv_1d, dICv_2d, dICv_3d, dICv_4d, dICv_5d,
           dH_IC, dH_ICv_1d, dH_ICv_2d, dH_ICv_3d, dH_ICv_4d, dH_ICv_5d,
           dD, 
           dR, dRv_1d, dRv_2d, dRv_3d, dRv_4d, dRv_5d,
           dR_1w, dRv_1d_1w, dRv_2d_1w, dRv_3d_1w, dRv_4d_1w, dRv_5d_1w,
           dR_2w, dRv_1d_2w, dRv_2d_2w, dRv_3d_2w, dRv_4d_2w, dRv_5d_2w,
           dR_3w, dRv_1d_3w, dRv_2d_3w, dRv_3d_3w, dRv_4d_3w, dRv_5d_3w
    ))
  })
}
