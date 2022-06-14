#' Age-structured SEIR ODE model test code (for debugging)
#' @param times vector of times
#' @param init list of initial states
#' @param params list of parameter values
#' @return List of summary results
#' @keywords vacamole
#' @export
# Define model -----------------------------------------------------
age_struct_seir_ode_test <- function(times, init, params) {
  with(as.list(c(params, init)), {
    # print(t)
    # define initial state vectors from input ----------------------
    # susceptible
    S <- c(S1, S2, S3, S4, S5, S6, S7, S8, S9)
    # exposed
    E <- c(E1, E2, E3, E4, E5, E6, E7, E8, E9)
    # infectious
    I <- c(I1, I2, I3, I4, I5, I6, I7, I8, I9)
    # hospitalized
    H <- c(H1, H2, H3, H4, H5, H6, H7, H8, H9)
    # ICU
    IC <- c(IC1, IC2, IC3, IC4, IC5, IC6, IC7, IC8, IC9)
    # return to hospital ward after ICU
    H_IC <- c(H_IC1, H_IC2, H_IC3, H_IC4, H_IC5, H_IC6, H_IC7, H_IC8, H_IC9)
    # death
    D <- c(D1, D2, D3, D4, D5, D6, D7, D8, D9)
    # recovered
    R <- c(R1, R2, R3, R4, R5, R6, R7, R8, R9)
    # extra recovered compartments so that waning immunity is not exponential
    R_1w <- c(R_1w1, R_1w2, R_1w3, R_1w4, R_1w5, R_1w6, R_1w7, R_1w8, R_1w9)
    R_2w <- c(R_2w1, R_2w2, R_2w3, R_2w4, R_2w5, R_2w6, R_2w7, R_2w8, R_2w9)
    R_3w <- c(R_3w1, R_3w2, R_3w3, R_3w4, R_3w5, R_3w6, R_3w7, R_3w8, R_3w9)

    
    # determine force of infection ----------------------------------
    calendar_day <- lubridate::yday(as.Date(times, origin = calendar_start_date))
    beta_t <- beta * (1 + beta1 * cos(2 * pi * calendar_day / 365.24)) # incorporate seasonality in transmission rate
    lambda <- beta_t * (contact_mat %*% I)
    # lambda <- ifelse(lambda < 0, 0, lambda)
    # print(beta_t)
    # ---------------------------------------------------------------
    
    ################################################################
    # ODEs:
    dS    <- -lambda * S + (omega * 4 * R_3w)
    dE    <- lambda * S - sigma * E + epsilon 
    dI    <- sigma * E - gamma * I - h * I
    dH    <- (h * I) - (i1 * H) - (d * H) - (r * H)
    dIC   <- (i1 * H) - (i2 * IC) - (d_ic * IC)
    dH_IC <- (i2 * IC) - (r_ic * H_IC) - (d_hic * H_IC)
    dD    <- (d * H) + (d_ic * IC) + (d_hic * H_IC)
    dR    <- (gamma * I) + (r * H) + (r_ic * H_IC) - (omega * 4 * R)
    dR_1w     <- (omega*4) * R - (omega*4) * R_1w
    dR_2w     <- (omega*4) * R_1w - (omega*4) * R_2w
    dR_3w     <- (omega*4) * R_2w - (omega*4) * R_3w
    ################################################################
    # dt <- 1
    list(c(dt, dS, dE, dI, dH, dIC, dH_IC, dD, dR, dR_1w, dR_2w, dR_3w))
  })
}
