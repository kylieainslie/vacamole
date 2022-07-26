# Minimal working example of SEIR model with hosptitalised and IC compartments
# each compartment is divided into 9 age groups
# Problem:
# With specified initial conditions, the ODE solver produces negative values
# in the IC compartment

# Load required packages -------------------------------------------
library(deSolve)
library(lubridate)
library(readxl)
library(ggplot2)

source("R/convert_vac_schedule2.R")
source("R/na_to_zero.R")
source("R/calc_waning.R")
# Define model -----------------------------------------------------
age_struct_seir_ode_test <- function(times, init, params) {
  with(as.list(c(params, init)), {
    # print(t)
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
    # incorporate seasonality in transmission rate 
    calendar_day <- lubridate::yday(as.Date(times, origin = calendar_start_date))
    beta_t <- beta * (1 + beta1 * cos(2 * pi * calendar_day / 365.24)) 
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
# -------------------------------------------------------------------
times <- seq(697, 719, by = 1)

# Specify initial conditions ----------------------------------------
empty_state <- c(rep(0,9))
n_vec <- 17407585 * c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 0.14514332, 0.12092904, 0.08807406, 0.04622194) # Dutch population size
s_vec <- c(1757166, 878617.4, 776667.4, 747091.8, 724414.8, 762523.6, 664041.6, 247366.5, 169643)
shold1_vec <- c(0, 20962.29,  14818.1,  13510.63,  7585.556,  4515.087,  1304.304,  1535.003)
sv1_vec <- c(0, 639545.1, 719967.7, 710214.3, 727050.4, 801066, 539380.6, 307956.5, 195745.2, 715.6142)
shold2_vec <- c(0, 5562.067, 4672.31, 4615.642, 3337.849, 2365.948, 676.7691, 1100.721, 389.9541)
sv2_vec <- c(0, 287791.2, 405376.3, 544364.5, 705437.3, 829093.4, 784944.8, 905975.1, 336207.8)
shold3_vec <- c(0, 28.01181, 1322.589, 2059.032, 3306.735, 4440.677, 2726.498, 1598.139, 12506.28)
sv3_vec <- c(33.65339, 1588.34, 2464.115, 3959.912, 5325.896, 3333.597, 2071.44, 16812.62)
shold4_vec <- empty_state
sv4_vec <- empty_state
shold5_vec <- empty_state
sv5_vec <- empty_state
e_vec <- c(2386.128, 7805.516, 10981.63, 3771.885, 3769.941, 4241.147, 3921.587, 1520.052, 2004.645)
ev1_vec <- c(0, 2533.785, 4544.659, 1627.346, 1727.1, 2026.257, 1795.951, 924.0714, 1187.933)
ev2_vec <- c(0, 358.4641, 884.8641, 523.5487, 710.9326, 934.8989, 1063.881, 1735.172, 1376.026)
ev3_vec <- c(0, 0.01708742, 1.286057, 0.7167056, 1.198273, 1.732874, 1.166555, 0.7673278, 12.10313)
ev4_vec <- empty_state
ev5_vec <- empty_state
i_vec <- c(2393.964, 7920.415, 11180.77, 3813.357, 3814.714, 4305.289, 3987.789, 1564.516, 2060.12)
iv1_vec <- c(0, 2546.347, 4583.856, 1634.677, 1739.419, 2048.863, 1820.342, 947.9065, 1213.283)
iv2_vec <- c(0, 357.1614, 886.2606, 524.1614, 714.3546, 943.9066, 1076.205, 1775.135, 1417.425)
iv3_vec <- c(0, 0.01060861, 0.7982189, 0.4436163, 0.742976, 1.078147, 0.7371704, 0.5077865, 8.15221)
iv4_vec <- empty_state
iv5_vec <- empty_state
h_vec <- c(27.18996, 3.716007, 13.60442, 16.03164, 27.18963, 62.6223, 84.38554, 77.09082, 147.9043)
hv1_vec <- c(0, 1.034427, 5.01926, 6.210725, 11.40022, 27.75593, 36.61349, 43.12757, 81.65341)
hv2_vec <- c(0, 0.1555604, 0.9425662, 1.977678, 4.445326, 12.09392, 22.06691, 85.16457, 99.29292)
hv3_vec <- c(0, 9.360707e-07, 0.0001901773, 0.0003795419, 0.001168993, 0.00365982, 0.004046525, 
             0.007012951, 0.1429213)
hv4_vec <- empty_state
hv5_vec <- empty_state
ic_vec <- c(0, 1.244071, 7.19287, 9.327737, 23.51088, 66.66829, 116.7801, 79.60983, 26.97763)
icv1_vec <- c(0, 0.326927, 2.515957, 3.477652, 9.576357, 28.83535, 49.522, 43.35214, 14.3759)
icv2_vec <- c(0, 0.04753805, 0.4588652, 1.093194, 3.690472, 12.43055, 29.61921, 85.37708, 17.49398)
icv3_vec <- c(0, 3.507978e-08, 1.100045e-05, 2.451984e-05, 0.0001134246, 0.0004429778, 0.0006808471, 
              0.001195892, 0.004727641)
icv4_vec <- empty_state
icv5_vec <- empty_state
hic_vec <- c(0, 0.5427055, 3.130428, 3.977067, 9.995462, 27.00499, 6.750479, 21.76603, 7.087477)
hicv1_vec <- c(0, 0.1371999,  1.060901,  1.448596,  4.000353,  11.50886,  2.826694, 11.66692, 3.700065)
hicv2_vec <- c(0, 0.01940543, 0.1896714, 0.4507069,  1.528743, 4.926877,   1.68333, 22.94502, 4.496207)
hicv3_vec <- c(0, 2.215896e-09, 6.81384e-07, 1.44918e-06, 6.776026e-06, 2.589946e-05, 7.141737e-06,
               7.790817e-05, 0.0003174733)
hicv4_vec <- empty_state
hicv5_vec <- empty_state
d_vec <- c(1.020095, 2.568894, 16.95339, 36.2586, 84.95562, 317.2533, 2170.362, 2172.682, 1893.923)
r_vec <- c(12965.45, 45395.17, 64970.73, 21163.58, 20838.08, 23145.87, 20872.23, 7944.972, 10673.15)
rv1_vec <- c(0, 13853.25, 25327.69, 8763.372, 9265.962, 10788.17, 9370.418, 4709.47, 6088.05)
rv2_vec <- c(0, 1877.575, 4747.051, 2769.703, 3756.09, 4910.486, 5481.361, 8757.952, 7164.384)
rv3_vec <- c(0, 0.00911274, 0.6822761, 0.3734436, 0.6233937, 0.9009441, 0.6505618, 0.5101846, 8.822465)
rv4_vec <- empty_state
rv5_vec <- empty_state
r_vec1 <- c(9858.296, 36074.89, 51848.93, 16378.41, 15986.89, 17615.19, 15649.99, 5898.077, 7984.347)
rv1_vec1 <- c(0,10549.48, 19453.3, 6602.408, 6972.201, 8082.101, 6939.408,  3440.285,  4450.496)
rv2_vec1 <- c(0, 1387.75, 3556.726, 2060.699, 2794.943, 3641.184, 4023.562,  6361.184 , 5216.777)
rv3_vec1 <- c(0, 0.0009841821, 0.0729049, 0.03891708, 0.06523884, 0.09528294, 0.07852455, 0.07795941,  1.460834)
rv4_vec1 <- empty_state
rv5_vec1 <- empty_state
r_vec2 <- c(6894.898, 26380.91, 37413.62, 11562.24, 11216.42, 12262.62, 10721.59, 3993.396, 5414.175)
rv1_vec2 <- c(0, 7381.78, 13584.93, 4553.088, 4810.678, 5550.896, 4705.327, 2297.737, 2961.422)
rv2_vec2 <- c(0, 935.1921, 2420.256, 1398.805,  1903.739, 2474.218,  2705.909, 4226.689, 3452.321)
rv3_vec2 <- c(0, 9.195519e-05, 0.00666582, 0.003381631, 0.005725258, 0.008557778, 0.008865433, 0.01145796, 0.2290288)
rv4_vec2 <- empty_state
rv5_vec2 <- empty_state
r_vec3 <- c(4755.703, 19319.74, 26158.76, 7973.445, 7707.248, 8350.497, 7141.973, 2613.521, 3517.111)
rv1_vec3 <- c(0, 5078.007,  9178.228,  3063.103, 3252.445,  3731.313,  3105.754,   1484.29,  1890.89)
rv2_vec3 <- c(0, 605.4878,  1575.553,  914.6719,   1258.31,  1636.313,  1767.933,  2714.133, 2191.516)
rv3_vec3 <- c(0, 8.160795e-06, 0.0005707839, 0.0002642459, 0.0004561424, 0.0007115875, 0.001003167, 
              0.001605036, 0.03335941)
rv4_vec3 <- empty_state
rv5_vec3 <- empty_state
# r_vec3 <- n_vec - s_vec - e_vec - i_vec - h_vec - hic_vec - ic_vec - d_vec - r_vec - r_vec1 - r_vec2

init_cond <- c(t = times[1],
               S = s_vec, Sv_1d = sv1_vec, Sv_2d = sv2_vec, Sv_3d = sv3_vec, Sv_4d = sv4_vec, Sv_5d = sv5_vec,
               Shold_1d = shold1_vec, Shold_2d = shold2_vec, Shold_3d = shold3_vec, Shold_4d = shold4_vec, Shold_5d = shold5_vec, 
               E = e_vec, Ev_1d = ev1_vec, Ev_2d = ev2_vec, Ev_3d = ev3_vec, Ev_4d = ev4_vec, Ev_5d = ev5_vec,
               I = i_vec, Iv_1d = iv1_vec, Iv_2d = iv2_vec, Iv_3d = iv3_vec, Iv_4d = iv4_vec, Iv_5d = iv5_vec,
               H = h_vec, Hv_1d = hv1_vec, Hv_2d = hv2_vec, Hv_3d = hv3_vec, Hv_4d = hv4_vec, Hv_5d = hv5_vec,
               IC = ic_vec, ICv_1d = icv1_vec, ICv_2d = icv2_vec, ICv_3d = icv3_vec, ICv_4d = icv4_vec, ICv_5d = icv5_vec,
               H_IC = hic_vec, H_ICv_1d = hicv1_vec, H_ICv_2d = hicv2_vec, H_ICv_3d = hicv3_vec, H_ICv_4d = hicv4_vec, H_ICv_5d = hicv5_vec,
               D = d_vec,
               R = r_vec, Rv_1d = rv1_vec, Rv_2d = rv2_vec, Rv_3d = rv3_vec, Rv_4d = rv4_vec, Rv_5d = rv5_vec,
               R_1w = r_vec1, Rv_1d_1w = rv1_vec1, Rv_2d_1w = rv2_vec1, Rv_3d_1w = rv3_vec1, Rv_4d_1w = rv4_vec1, Rv_5d_1w = rv5_vec1,
               R_2w = r_vec2, Rv_1d_2w = rv1_vec2, Rv_2d_2w = rv2_vec2, Rv_3d_2w = rv3_vec2, Rv_4d_2w = rv4_vec2, Rv_5d_2w = rv5_vec2,
               R_3w = r_vec3, Rv_1d_3w = rv1_vec3, Rv_2d_3w = rv2_vec3, Rv_3d_3w = rv3_vec3, Rv_4d_3w = rv4_vec3, Rv_5d_3w = rv5_vec3
)

# Specify model parameters ------------------------------------------
# define contact/transmission matrix
path <- "inst/extdata/inputs/contact_matrices/converted/"
# path <- "/rivm/s/ainsliek/data/contact_matrices/converted/" 
november_2021  <- readRDS(paste0(path,"transmission_matrix_november_2021.rds"))

# probabilities -----------------------------------------------------
p_infection2admission <- c(0.00347, 0.000377, 0.000949, 0.00388, 0.00842, 0.0165, 0.0251, 0.0494, 0.0463)
p_admission2death     <- c(0.00191, 0.00433, 0.00976, 0.0219, 0.025, 0.0401, 0.106, 0.229, 0.311)
p_admission2IC        <- c(0, 0.0618, 0.0962, 0.11, 0.164, 0.202, 0.244, 0.196, 0.035)
p_IC2hospital         <- c(0.866, 0.866, 0.866, 0.866, 0.866, 0.829, 0.120, 0.573, 0.555)
p_hospital2death      <- c(rep(0, 5), 0.01, 0.04, 0.12, 0.29) # (after ICU)
p_reported_by_age     <- c(0.29, 0.363, 0.381, 0.545, 0.645, 0.564, 0.365, 0.33, 0.409) # from Jantien

# delays ------------------------------------------------------------
time_symptom2admission   <- c(2.29, 5.51, 5.05, 5.66, 6.55, 5.88, 5.69, 5.09, 4.33) # assume same as infectious2admission
time_admission2discharge <- 7.9
time_admission2IC        <- 2.28
time_IC2hospital         <- 15.6
time_hospital2discharge  <- 10.1 # (after ICU)
time_admission2death     <- 7
time_IC2death            <- 19
time_hospital2death      <- 10 # (after ICU)

# define transition rates -------------------------------------------
i2r    <- (1-p_infection2admission) / 2                   # I -> R
i2h    <- p_infection2admission / time_symptom2admission  # I -> H

h2ic   <- p_admission2IC / time_admission2IC              # H -> IC
h2d    <- p_admission2death / time_admission2death        # H -> D
h2r    <- (1 - (p_admission2IC + p_admission2death)) / time_admission2discharge
                                                          # H -> R

ic2hic <- p_IC2hospital / time_IC2hospital                # IC -> H_IC
ic2d   <- (1 - p_IC2hospital) / time_IC2death             # IC -> D

hic2d  <- p_hospital2death / time_hospital2death          # H_IC -> D
hic2r  <- (1 - p_hospital2death) / time_hospital2discharge# H_IC -> R

# vaccination schedule ----------------------------------------------
# read in vaccination schedule
raw_vac_schedule <- read_csv("inst/extdata/inputs/vac_schedule_real_w_4th_and_5th_dose.csv") %>%
  select(-X1)

# subset for only pfizer and 2 doses
# pf_schedule <- vac_schedule %>%
#   select(date:pf_d1_9, pf_d2_1:pf_d2_9) 

# read in xlsx file with VEs (there is 1 sheet for each variant)
# we'll only use wildtype values for now
# wt_ve <- read_excel("inst/extdata/inputs/ve_dat.xlsx", sheet = "wildtype") 
# alpha_ve <- read_excel("inst/extdata/inputs/ve_dat.xlsx", sheet = "alpha")
# delta_ve <- read_excel("inst/extdata/inputs/ve_dat.xlsx", sheet = "delta")
omicron_ve <- read_excel("inst/extdata/inputs/ve_estimates/ve_dat.xlsx", sheet = "omicron")

# convert vaccination schedule for input into model
vac_rates_wt <- convert_vac_schedule_debug(
  vac_schedule = raw_vac_schedule,
  ve_pars = wt_ve,
  wane = TRUE)
                       
# data wrangle for model input
vac_rates <- convert_vac_schedule2(
  vac_schedule = vac_schedule_extra,
  ve_pars = ve_params,
  wane = TRUE,
  k_inf = 0.012,
  k_sev = 0.006)

# data wrangle for model input
df_input <- pivot_wider(vac_rates %>% 
                          filter(param != "comp_ve") %>%
                          mutate(param = ifelse(param == "comp_delay", "delay", param)), 
                        names_from = c("param", "age_group"), 
                        names_sep = "", values_from = "value")

# parameters must be in a named list
params <- list(N = n_vec,
               beta = 0.0004,
               beta1 = 0.14,
               sigma = 0.5,
               gamma = i2r,
               h = i2h,
               i1 = h2ic,
               d = empty_state, #h2d,
               r = h2r,
               i2 = ic2hic,
               d_ic = empty_state, #ic2d,
               d_hic = hic2d,
               r_ic = hic2r,
               epsilon = 0.00,
               omega = wane_8months,
               # daily vaccination rate
               alpha1 = df_input %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha2 = df_input %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha3 = df_input %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha4 = df_input %>%
                 filter(dose == "d4", outcome == "infection") %>%
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha5 = df_input %>%
                 filter(dose == "d5", outcome == "infection") %>%
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               # delay to protection
               delay1 = df_input %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay2 = df_input %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay3 = df_input %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay4 = df_input %>%
                 filter(dose == "d4", outcome == "infection") %>%
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay5 = df_input %>%
                 filter(dose == "d5", outcome == "infection") %>%
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               # protection against infection
               eta1 = df_input %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta2 = df_input %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta3 = df_input %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta4 = df_input %>%
                 filter(dose == "d4", outcome == "infection") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta5 = df_input %>%
                 filter(dose == "d5", outcome == "infection") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               # protection from hospitalisation
               eta_hosp1 = df_input %>% 
                 filter(dose == "d1", outcome == "hospitalisation") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp2 = df_input %>% 
                 filter(dose == "d2", outcome == "hospitalisation") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp3 = df_input %>% 
                 filter(dose == "d3", outcome == "hospitalisation") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp4 = df_input %>%
                 filter(dose == "d4", outcome == "hospitalisation") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp5 = df_input %>%
                 filter(dose == "d5", outcome == "hospitalisation") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               # protection from transmission
               eta_trans1 = df_input %>% 
                 filter(dose == "d1", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans2 = df_input %>% 
                 filter(dose == "d2", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans3 = df_input %>% 
                 filter(dose == "d3", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans4 = df_input %>%
                 filter(dose == "d4", outcome == "transmission") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans5 = df_input %>%
                 filter(dose == "d5", outcome == "transmission") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               p_report = p_reported_by_age,
               contact_mat = november_2021,#$mean,  # change contact matrix
               calendar_start_date = as.Date("2020-01-01")#,
               #no_vac = nv
)
# Run model -----------------------------------------------------------
rk45 <- rkMethod("rk45dp7")
seir_out <- ode(init_cond, times, age_struct_seir_ode_test, params, method = rk45) #, rtol = 1e-08, hmax = 0.02
out <- as.data.frame(seir_out) 
# ---------------------------------------------------------------------

# check for negative values
any(tail(seir_out, 1) < 0)

# check population size
all.equal(sum(tail(seir_out, 1)[-c(1:2)]),sum(params$N))

# Plot output ---------------------------------------------------------
out_long <- out %>% 
  select(-t) %>%
  pivot_longer(-c(time), names_to = "state", values_to = "value") %>%
  separate(state, into = c("state", "age_group"), sep = -1) %>%
  group_by(time, state) %>%
  summarise(total = sum(value)) %>%
  ungroup()

p <- ggplot(out_long %>% filter(grepl('Ev_1d', state)), 
            aes(x = time, y = total, color = state)) +
  geom_line()
p
# --------------------------------------------------------------------




