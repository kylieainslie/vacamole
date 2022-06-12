# Minimal working example of SEIR model with hosptitalised and IC compartments
# each compartment is divided into 9 age groups
# Problem:
# With specified initial conditions, the ODE solver produces negative values
# in the IC compartment

# Load required packages -------------------------------------------
library(deSolve)
library(lubridate)
library(readxl)

source("R/convert_vac_schedule_debug.R")
source("R/na_to_zero.R")
source("R/calc_waning.R")
# Define model -----------------------------------------------------
age_struct_seir_ode_test <- function(times, init, params) {
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
    # exposed
    E <- c(E1, E2, E3, E4, E5, E6, E7, E8, E9)     
    Ev_1d <- c(Ev_1d1, Ev_1d2, Ev_1d3, Ev_1d4, Ev_1d5, Ev_1d6, Ev_1d7, Ev_1d8, Ev_1d9)
    Ev_2d <- c(Ev_2d1, Ev_2d2, Ev_2d3, Ev_2d4, Ev_2d5, Ev_2d6, Ev_2d7, Ev_2d8, Ev_2d9)
    Ev_3d <- c(Ev_3d1, Ev_3d2, Ev_3d3, Ev_3d4, Ev_3d5, Ev_3d6, Ev_3d7, Ev_3d8, Ev_3d9)
    # infectious
    I <- c(I1, I2, I3, I4, I5, I6, I7, I8, I9)     
    Iv_1d <- c(Iv_1d1, Iv_1d2, Iv_1d3, Iv_1d4, Iv_1d5, Iv_1d6, Iv_1d7, Iv_1d8, Iv_1d9)
    Iv_2d <- c(Iv_2d1, Iv_2d2, Iv_2d3, Iv_2d4, Iv_2d5, Iv_2d6, Iv_2d7, Iv_2d8, Iv_2d9)
    Iv_3d <- c(Iv_3d1, Iv_3d2, Iv_3d3, Iv_3d4, Iv_3d5, Iv_3d6, Iv_3d7, Iv_3d8, Iv_3d9)
    # hospitalized
    H <- c(H1, H2, H3, H4, H5, H6, H7, H8, H9)     
    Hv_1d <- c(Hv_1d1, Hv_1d2, Hv_1d3, Hv_1d4, Hv_1d5, Hv_1d6, Hv_1d7, Hv_1d8, Hv_1d9)
    Hv_2d <- c(Hv_2d1, Hv_2d2, Hv_2d3, Hv_2d4, Hv_2d5, Hv_2d6, Hv_2d7, Hv_2d8, Hv_2d9)
    Hv_3d <- c(Hv_3d1, Hv_3d2, Hv_3d3, Hv_3d4, Hv_3d5, Hv_3d6, Hv_3d7, Hv_3d8, Hv_3d9)
    # ICU
    IC <- c(IC1, IC2, IC3, IC4, IC5, IC6, IC7, IC8, IC9)
    ICv_1d <- c(ICv_1d1, ICv_1d2, ICv_1d3, ICv_1d4, ICv_1d5, ICv_1d6, ICv_1d7, ICv_1d8, ICv_1d9)
    ICv_2d <- c(ICv_2d1, ICv_2d2, ICv_2d3, ICv_2d4, ICv_2d5, ICv_2d6, ICv_2d7, ICv_2d8, ICv_2d9)
    ICv_3d <- c(ICv_3d1, ICv_3d2, ICv_3d3, ICv_3d4, ICv_3d5, ICv_3d6, ICv_3d7, ICv_3d8, ICv_3d9)
    # return to hospital ward after ICU
    H_IC <- c(H_IC1, H_IC2, H_IC3, H_IC4, H_IC5, H_IC6, H_IC7, H_IC8, H_IC9)
    H_ICv_1d <- c(H_ICv_1d1, H_ICv_1d2, H_ICv_1d3, H_ICv_1d4, H_ICv_1d5, H_ICv_1d6, H_ICv_1d7, H_ICv_1d8, H_ICv_1d9)
    H_ICv_2d <- c(H_ICv_2d1, H_ICv_2d2, H_ICv_2d3, H_ICv_2d4, H_ICv_2d5, H_ICv_2d6, H_ICv_2d7, H_ICv_2d8, H_ICv_2d9)
    H_ICv_3d <- c(H_ICv_3d1, H_ICv_3d2, H_ICv_3d3, H_ICv_3d4, H_ICv_3d5, H_ICv_3d6, H_ICv_3d7, H_ICv_3d8, H_ICv_3d9)
    # death
    D <- c(D1, D2, D3, D4, D5, D6, D7, D8, D9) 
    # recovered
    R <- c(R1, R2, R3, R4, R5, R6, R7, R8, R9)     
    Rv_1d <- c(Rv_1d1, Rv_1d2, Rv_1d3, Rv_1d4, Rv_1d5, Rv_1d6, Rv_1d7, Rv_1d8, Rv_1d9)
    Rv_2d <- c(Rv_2d1, Rv_2d2, Rv_2d3, Rv_2d4, Rv_2d5, Rv_2d6, Rv_2d7, Rv_2d8, Rv_2d9)
    Rv_3d <- c(Rv_3d1, Rv_3d2, Rv_3d3, Rv_3d4, Rv_3d5, Rv_3d6, Rv_3d7, Rv_3d8, Rv_3d9)
    # extra recovered compartments so that waning immunity is not exponential
    R_1w <- c(R_1w1, R_1w2, R_1w3, R_1w4, R_1w5, R_1w6, R_1w7, R_1w8, R_1w9)
    Rv_1d_1w <- c(Rv_1d_1w1, Rv_1d_1w2, Rv_1d_1w3, Rv_1d_1w4, Rv_1d_1w5, Rv_1d_1w6, Rv_1d_1w7, Rv_1d_1w8, Rv_1d_1w9)
    Rv_2d_1w <- c(Rv_2d_1w1, Rv_2d_1w2, Rv_2d_1w3, Rv_2d_1w4, Rv_2d_1w5, Rv_2d_1w6, Rv_2d_1w7, Rv_2d_1w8, Rv_2d_1w9)
    Rv_3d_1w <- c(Rv_3d_1w1, Rv_3d_1w2, Rv_3d_1w3, Rv_3d_1w4, Rv_3d_1w5, Rv_3d_1w6, Rv_3d_1w7, Rv_3d_1w8, Rv_3d_1w9)
    
    R_2w <- c(R_2w1, R_2w2, R_2w3, R_2w4, R_2w5, R_2w6, R_2w7, R_2w8, R_2w9)
    Rv_1d_2w <- c(Rv_1d_2w1, Rv_1d_2w2, Rv_1d_2w3, Rv_1d_2w4, Rv_1d_2w5, Rv_1d_2w6, Rv_1d_2w7, Rv_1d_2w8, Rv_1d_2w9)
    Rv_2d_2w <- c(Rv_2d_2w1, Rv_2d_2w2, Rv_2d_2w3, Rv_2d_2w4, Rv_2d_2w5, Rv_2d_2w6, Rv_2d_2w7, Rv_2d_2w8, Rv_2d_2w9)
    Rv_3d_2w <- c(Rv_3d_2w1, Rv_3d_2w2, Rv_3d_2w3, Rv_3d_2w4, Rv_3d_2w5, Rv_3d_2w6, Rv_3d_2w7, Rv_3d_2w8, Rv_3d_2w9)
    
    R_3w <- c(R_3w1, R_3w2, R_3w3, R_3w4, R_3w5, R_3w6, R_3w7, R_3w8, R_3w9)
    Rv_1d_3w <- c(Rv_1d_3w1, Rv_1d_3w2, Rv_1d_3w3, Rv_1d_3w4, Rv_1d_3w5, Rv_1d_3w6, Rv_1d_3w7, Rv_1d_3w8, Rv_1d_3w9)
    Rv_2d_3w <- c(Rv_2d_3w1, Rv_2d_3w2, Rv_2d_3w3, Rv_2d_3w4, Rv_2d_3w5, Rv_2d_3w6, Rv_2d_3w7, Rv_2d_3w8, Rv_2d_3w9)
    Rv_3d_3w <- c(Rv_3d_3w1, Rv_3d_3w2, Rv_3d_3w3, Rv_3d_3w4, Rv_3d_3w5, Rv_3d_3w6, Rv_3d_3w7, Rv_3d_3w8, Rv_3d_3w9)
    
    # define vaccination parameters ---------------------------------
    index <- floor(t) + 1              # use floor of time point + 1 to index df
    alpha1 <- params$alpha1[index, -1] # remove date column
    alpha2 <- params$alpha2[index, -1]
    alpha3 <- params$alpha3[index, -1]
    delay1 <- params$delay1[index, -1]
    delay2 <- params$delay2[index, -1]
    delay3 <- params$delay3[index, -1]
    eta1   <- params$eta1[index, -1]
    eta2   <- params$eta2[index, -1]
    eta3   <- params$eta3[index, -1]
    eta_hosp1   <- params$eta1[index, -1]
    eta_hosp2   <- params$eta2[index, -1]
    eta_hosp3   <- params$eta3[index, -1]
    eta_trans1   <- as.numeric(params$eta_trans1[index, -1])
    eta_trans2   <- as.numeric(params$eta_trans2[index, -1])
    eta_trans3   <- as.numeric(params$eta_trans3[index, -1])
    # ---------------------------------------------------------------
    
    # determine force of infection ----------------------------------
    # incorporate seasonality in transmission rate 
    calendar_day <- lubridate::yday(as.Date(times, origin = calendar_start_date))
    beta_t <- beta * (1 + beta1 * cos(2 * pi * calendar_day / 365.24)) 
    lambda <- beta_t * (contact_mat %*% (I + (eta_trans1 * Iv_1d) + (eta_trans2 * Iv_2d) + (eta_trans3 * Iv_3d)))
    # ---------------------------------------------------------------
    
    #################################################################
    # ODEs:
    dS <- -lambda * S - alpha1 * S + (omega*4) * R_3w
    dShold_1d <- alpha1 * S - delay1 * Shold_1d - lambda * Shold_1d
    dSv_1d <- delay1 * Shold_1d - eta1 * lambda * Sv_1d - alpha2 * Sv_1d + (omega*4) * Rv_1d_3w
    dShold_2d <- alpha2 * Sv_1d - delay2 * Shold_2d - eta1 * lambda * Shold_2d
    dSv_2d <- delay2 * Shold_2d - eta2 * lambda * Sv_2d - alpha3 * Sv_2d + (omega*4) * Rv_2d_3w
    dShold_3d <- alpha3 * Sv_2d - delay3 * Shold_3d - eta2 * lambda * Shold_3d
    dSv_3d <- delay3 * Shold_3d - eta3 * lambda * Sv_3d + (omega*4) * Rv_3d_3w
    
    dE     <- lambda * S + lambda * Shold_1d - sigma * E + epsilon
    dEv_1d <- eta1 * lambda * Sv_1d + eta1 * lambda * Shold_2d - sigma * Ev_1d
    dEv_2d <- eta2 * lambda * Sv_2d + eta2 * lambda * Shold_3d - sigma * Ev_2d
    dEv_2d <- eta3 * lambda * Sv_3d - sigma * Ev_3d
    
    dI     <- sigma * E - (gamma + h) * I
    dIv_1d <- sigma * Ev_1d - (gamma + eta_hosp1 * h) * Iv_1d
    dIv_2d <- sigma * Ev_2d - (gamma + eta_hosp2 * h) * Iv_2d
    dIv_3d <- sigma * Ev_3d - (gamma + eta_hosp3 * h) * Iv_3d
    
    dH     <- h * I - (i1 + d + r) * H
    dHv_1d <- eta_hosp1 * h * Iv_1d - (i1 + d + r) * Hv_1d
    dHv_2d <- eta_hosp2 * h * Iv_2d - (i1 + d + r) * Hv_2d
    dHv_3d <- eta_hosp3 * h * Iv_3d - (i1 + d + r) * Hv_3d
    
    dIC     <- i1 * H - (i2 + d_ic) * IC
    dICv_1d <- i1 * Hv_1d - (i2 + d_ic) * ICv_1d
    dICv_2d <- i1 * Hv_2d - (i2 + d_ic) * ICv_2d
    dICv_3d <- i1 * Hv_3d - (i2 + d_ic) * ICv_3d
    
    dH_IC     <- i2 * IC - (r_ic + d_hic) * H_IC
    dH_ICv_1d <- i2 * ICv_1d - (r_ic + d_hic) * H_ICv_1d
    dH_ICv_2d <- i2 * ICv_2d - (r_ic + d_hic) * H_ICv_2d
    dH_ICv_3d <- i2 * ICv_3d - (r_ic + d_hic) * H_ICv_3d
    
    dD <- d * (H + Hv_1d + Hv_2d + Hv_3d) + 
      d_ic * (IC + ICv_1d + ICv_2d + ICv_3d) + 
      d_hic * (H_IC + H_ICv_1d + H_ICv_2d + H_ICv_3d)
    
    dR     <- (gamma * I) + (r * H) + (r_ic * H_IC) - ((omega*4) * R)
    dRv_1d <- (gamma * Iv_1d) + (r * Hv_1d) + (r_ic * H_ICv_1d) - ((omega*4) * Rv_1d)
    dRv_2d <- (gamma * Iv_2d) + (r * Hv_2d) + (r_ic * H_ICv_2d) - ((omega*4) * Rv_2d)
    dRv_3d <- (gamma * Iv_3d) + (r * Hv_3d) + (r_ic * H_ICv_3d) - ((omega*4) * Rv_3d)
    
    dR_1w     <- (omega*4) * R - (omega*4) * R_1w
    dRv_1d_1w <- (omega*4) * Rv_1d - (omega*4) * Rv_1d_1w
    dRv_2d_1w <- (omega*4) * Rv_2d - (omega*4) * Rv_2d_1w
    dRv_3d_1w <- (omega*4) * Rv_3d - (omega*4) * Rv_3d_1w
    
    dR_2w     <- (omega*4) * R_1w - (omega*4) * R_2w
    dRv_1d_2w <- (omega*4) * Rv_1d_1w - (omega*4) * Rv_1d_2w
    dRv_2d_2w <- (omega*4) * Rv_2d_1w - (omega*4) * Rv_2d_2w
    dRv_3d_2w <- (omega*4) * Rv_3d_1w - (omega*4) * Rv_3d_2w
    
    dR_3w     <- (omega*4) * R_2w - (omega*4) * R_3w
    dRv_1d_3w <- (omega*4) * Rv_1d_2w - (omega*4) * Rv_1d_3w
    dRv_2d_3w <- (omega*4) * Rv_2d_2w - (omega*4) * Rv_2d_3w
    dRv_3d_3w <- (omega*4) * Rv_3d_2w - (omega*4) * Rv_3d_3w
    #################################################################
    dt <- 1
    # output --------------------------------------------------------
    list(c(dt, dS, dShold_1d, dSv_1d, dShold_2d, dSv_2d, dShold_3d, dSv_3d, 
           dE, dEv_1d, dEv_2d, dEv_3d,
           dI, dIv_1d, dIv_2d, dIv_3d,
           dH, dHv_1d, dHv_2d, dHv_3d,
           dIC, dICv_1d, dICv_2d, dICv_3d,
           dH_IC, dH_ICv_1d, dH_ICv_2d, dH_ICv_3d,
           dD, 
           dR, dRv_1d, dRv_2d, dRv_3d,
           dR_1w, dRv_1d_1w, dRv_2d_1w, dRv_3d_1w,
           dR_2w, dRv_1d_2w, dRv_2d_2w, dRv_3d_2w,
           dR_3w, dRv_1d_3w, dRv_2d_3w, dRv_3d_3w))
  })
}
# -------------------------------------------------------------------
times <- seq(369, 400, by = 1)

# Specify initial conditions ----------------------------------------
empty_state <- c(rep(0,9))
n_vec <- 17407585 * c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 0.14514332, 0.12092904, 0.08807406, 0.04622194) # Dutch population size
s_vec <- c(1758200.1, 1795008.2, 1822164.1, 1939971.7, 2081005.6, 2305117.0, 1909464.3, 1401363.8, 696014.1)
shold1_vec <- empty_state
sv1_vec <- empty_state
shold2_vec <- empty_state
sv2_vec <- empty_state
shold3_vec <- empty_state
sv3_vec <- empty_state
e_vec <- c(629.0361, 3594.4890, 5915.3246, 2910.9371, 3130.5281, 3523.5570, 3116.9297, 2113.7468, 1700.1333)
ev1_vec <- empty_state
ev2_vec <- empty_state
ev3_vec <- empty_state
i_vec <- c(671.8259, 3846.0577, 6342.2980, 3119.7126, 3366.3184, 3809.1308, 3388.1608, 2331.4422, 1867.0974)
iv1_vec <- empty_state
iv2_vec <- empty_state
iv3_vec <- empty_state
h_vec <- c(9.747787, 2.181694, 9.196823, 15.838136, 28.679582, 65.825066, 84.881661, 138.027275, 169.278401)
hv1_vec <- empty_state
hv2_vec <- empty_state
hv3_vec <- empty_state
ic_vec <- c(0, 0.9713496,  6.5248700,  12.6613143,  34.6104265,  99.2475311, 175.8306823, 209.6044387, 44.6070593)
icv1_vec <- empty_state
icv2_vec <- empty_state
icv3_vec <- empty_state
hic_vec <- c(0, 0.5002954, 3.3952939, 6.5203019, 17.8682707, 49.1401382, 12.5871195, 71.4254831, 14.6514841)
hicv1_vec <- empty_state
hicv2_vec <- empty_state
hicv3_vec <- empty_state
d_vec <- c(0.2322485, 0.6282215, 4.8540845, 11.4236062, 28.2727435, 108.5131182, 755.6668547, 870.7532525, 824.5838766)
r_vec <- c(19070.29, 113274.53, 193344.72, 90546.75, 96975.90, 108942.40, 95565.06, 64079.59, 52759.36)
rv1_vec <- empty_state
rv2_vec <- empty_state
rv3_vec <- empty_state
r_vec1 <- c(10587.70, 63859.80, 111614.59,  51041.91,  54561.39,  61339.48,  53778.07, 35909.68, 29848.27)
rv1_vec1 <- empty_state
rv2_vec1 <- empty_state
rv3_vec1 <- empty_state
r_vec2 <- c(5073.434, 30457.485, 54493.620, 24804.111, 26513.001, 29958.204, 26456.079, 17704.939, 14638.817)
rv1_vec2 <- empty_state
rv2_vec2 <- empty_state
rv3_vec2 <- empty_state
r_vec3 <- c(2206.518, 12865.510, 23865.844, 11058.733, 11852.807, 13582.198, 12285.014, 8363.689, 6731.436)
rv1_vec3 <- empty_state
rv2_vec3 <- empty_state
# r_vec3 <- n_vec - s_vec - e_vec - i_vec - h_vec - hic_vec - ic_vec - d_vec - r_vec - r_vec1 - r_vec2

init_cond <- c(t = times[1],
               S = s_vec, Shold_1d = shold1_vec, Sv_1d = sv1_vec, Shold_2d = shold2_vec, Sv_2d = sv2_vec, Shold_3d = shold3_vec, Sv_3d = sv3_vec,
               E = e_vec, Ev_1d = ev1_vec, Ev_2d = ev2_vec, Ev_3d = ev3_vec,
               I = i_vec, Iv_1d = iv1_vec, Iv_2d = iv2_vec, Iv_3d = Iv3_vec,
               H = h_vec, Hv_1d = hv1_vec, Hv_2d = hv2_vec, Hv_3d = hv3_vec,
               IC = ic_vec, ICv_1d = icv1_vec, ICv_2d = icv2_vec, ICv_3d = icv3_vec,
               H_IC = hic_vec, H_ICv_1d = hicv1_vec, H_ICv_2d = hicv2_vec, H_ICv_3d = hicv3_vec,
               D = d_vec,
               R = r_vec, Rv_1d = rv1_vec, Rv_2d = rv2_vec, Rv_3d = rv3_vec,
               R_1w = r_vec1, Rv_1d_1w = rv1_vec1, Rv_2d_1w = rv2_vec1, Rv_3d_1w = rv3_vec1,
               R_2w = r_vec2, Rv_1d_2w = rv1_vec2, Rv_2d_2w = rv2_vec2, Rv_3d_2w = rv3_vec2,
               R_3w = r_vec3, Rv_1d_3w = rv1_vec3, Rv_2d_3w = rv2_vec3, Rv_3d_3w = rv3_vec3
)

# Specify model parameters ------------------------------------------
# define contact/transmission matrix
cm <- matrix(c(0.00008, 0.00004, 0.00004, 0.00003, 0.00003, 0.00001, 0.00001, 0.00001, 0.00001,
               0.00003, 0.00045, 0.00016, 0.00006, 0.00007, 0.00005, 0.00003, 0.00002, 0.00001,
               0.00003, 0.00015, 0.00049, 0.00014, 0.00011, 0.00010, 0.00007, 0.00004, 0.00006,
               0.00003, 0.00006, 0.00014, 0.00015, 0.00009, 0.00007, 0.00005, 0.00003, 0.00004,
               0.00002, 0.00006, 0.00011, 0.00009, 0.00014, 0.00008, 0.00006, 0.00004, 0.00004,
               0.00001, 0.00004, 0.00009, 0.00006, 0.00007, 0.00010, 0.00007, 0.00005, 0.00007,
               0.00001, 0.00003, 0.00007, 0.00006, 0.00006, 0.00009, 0.00019, 0.00011, 0.00012,
               0.00001, 0.00002, 0.00006, 0.00005, 0.00006, 0.00008, 0.00016, 0.00029, 0.00023,
               0.00001, 0.00004, 0.00017, 0.00010, 0.00013, 0.00023, 0.00032, 0.00043, 0.00120), nrow = 9)

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
vac_schedule <- read_csv("inst/extdata/inputs/vac_schedule_real_w_4th_and_5th_dose.csv") %>%
  select(-X1)

# subset for only pfizer and 2 doses
pf_schedule <- vac_schedule %>%
  select(date:pf_d1_9, pf_d2_1:pf_d2_9) 

# read in xlsx file with VEs (there is 1 sheet for each variant)
# we'll only use wildtype values for now
wt_ve <- read_excel("inst/extdata/inputs/ve_dat.xlsx", sheet = "wildtype") 
# alpha_ve <- read_excel("inst/extdata/inputs/ve_dat.xlsx", sheet = "alpha") 
# delta_ve <- read_excel("inst/extdata/inputs/ve_dat.xlsx", sheet = "delta") 
# omicron_ve <- read_excel("inst/extdata/inputs/ve_dat.xlsx", sheet = "omicron") 

# convert vaccination schedule for input into model
vac_rates_wt <- convert_vac_schedule_debug(
  vac_schedule = pf_schedule,
  ve_pars = wt_ve,
  wane = TRUE)
                       
# data wrangle for model input
df_input <- pivot_wider(vac_rates_wt %>% 
                          filter(param != "comp_ve") %>%
                          mutate(param = ifelse(param == "comp_delay", "delay", param)), 
                        names_from = c("param", "age_group"), 
                        names_sep = "", values_from = "value")

# parameters must be in a named list
params <- list(N = n_vec,
               beta = 0.0002531514,
               beta1 = 0.14,
               sigma = 0.5,
               gamma = i2r,
               h = i2h,
               i1 = h2ic,
               d = h2d,
               r = h2r,
               i2 = ic2hic,
               d_ic = ic2d,
               d_hic = hic2d,
               r_ic = hic2r,
               epsilon = 0.00,
               omega = 0.0038,
               alpha1 = df_input %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha2 = df_input %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha3 = df_input %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               delay1 = df_input %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay2 = df_input %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay3 = df_input %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               eta1 = df_input %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta2 = df_input %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta3 = df_input %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp1 = df_input %>% 
                 filter(dose == "d1", outcome == "hospitalisation") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp2 = df_input %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp3 = df_input %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans1 = df_input %>% 
                 filter(dose == "d1", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans2 = df_input %>% 
                 filter(dose == "d2", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans3 = df_input %>% 
                 filter(dose == "d3", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               p_report = p_reported_by_age,
               contact_mat = cm,
               calendar_start_date = as.Date("2020-01-01")
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
out_for_plot <- out %>% 
  select(-t) %>%
  pivot_longer(-c(time), names_to = "state", values_to = "value") %>%
  separate(state, into = c("state", "age_group"), sep = -1) %>%
  group_by(time, state) %>%
  summarise(total = sum(value)) %>%
  ungroup()

p <- ggplot(out_for_plot %>% filter(grepl('Rv_2d', state)), 
            aes(x = time, y = total, color = state)) +
  geom_line()
p
# --------------------------------------------------------------------




