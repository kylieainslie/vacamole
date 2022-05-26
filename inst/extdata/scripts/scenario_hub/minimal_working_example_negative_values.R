# Minimal working example of SEIR model with hosptitalised and IC compartments
# each compartment is divided into 9 age groups
# Problem:
# With specified initial conditions, the ODE solver produces negative values
# in the IC compartment

# Load required packages -------------------------------------------
library(deSolve)
library(lubridate)

# Define model -----------------------------------------------------
age_struct_seir_ode_test <- function(times, init, params) {
  with(as.list(c(params, init)), {
    # define initial state vectors from input ----------------------
    S <- c(S1, S2, S3, S4, S5, S6, S7, S8, S9)     # susceptible
    E <- c(E1, E2, E3, E4, E5, E6, E7, E8, E9)     # exposed
    I <- c(I1, I2, I3, I4, I5, I6, I7, I8, I9)     # infectious
    H <- c(H1, H2, H3, H4, H5, H6, H7, H8, H9)     # hospitalized
    IC <- c(IC1, IC2, IC3, IC4, IC5, IC6, IC7, 
            IC8, IC9)                              # ICU
    H_IC <- c(H_IC1, H_IC2, H_IC3, H_IC4, H_IC5, 
              H_IC6, H_IC7, H_IC8, H_IC9)          # return to hospital ward after ICU
    D <- c(D1, D2, D3, D4, D5, D6, D7, D8, D9)     # death
    R <- c(R1, R2, R3, R4, R5, R6, R7, R8, R9)     # recovered
    
    # determine force of infection ----------------------------------
    # incorporate seasonality in transmission rate 
    calendar_day <- lubridate::yday(as.Date(times, origin = calendar_start_date))
    beta_t <- beta * (1 + beta1 * cos(2 * pi * calendar_day / 365.24)) 
    lambda <- beta_t * (contact_mat %*% I)
    # ---------------------------------------------------------------
    
    #################################################################
    # ODEs:
    dS    <- -lambda * S + (omega * R)
    dE    <- lambda * S - sigma * E + epsilon 
    dI    <- sigma * E - gamma * I - h * I
    dH    <- (h * I) - (i1 * H) - (d * H) - (r * H)
    dIC   <- (i1 * H) - (i2 * IC) - (d_ic * IC)
    dH_IC <- (i2 * IC) - (r_ic * H_IC) - (d_hic * H_IC)
    dD    <- (d * H) + (d_ic * IC) + (d_hic * H_IC)
    dR    <- (gamma * I) - (omega * R) + (r * H) + (r_ic * H_IC)
    #################################################################
    
    # output --------------------------------------------------------
    list(c(dt, dS, dE, dI, dH, dIC, dH_IC, dD, dR))
  })
}
# -------------------------------------------------------------------

# Specify initial conditions ----------------------------------------
init_cond <- c(t = 76,
               S = c(1.796233e+06, 2.021641e+06, 2.214668e+06, 2.121977e+06, 2.275925e+06, 2.524675e+06, 2.103065e+06, 1.531171e+06, 8.030364e+05),
               E = c(3.707938e+01, 2.175486e+02, 5.307151e+02, 2.610720e+02, 2.724314e+02, 3.290284e+02, 3.457338e+02, 3.402400e+02, 2.699806e+02), 
               I = c(3.105244e+01, 1.826302e+02, 4.454958e+02, 2.189514e+02, 2.282486e+02, 2.749749e+02, 2.881797e+02, 2.811478e+02, 2.227498e+02), 
               H = c(2.091348e-01, 4.961192e-02, 3.137189e-01, 5.496294e-01, 9.891533e-01, 2.461984e+00, 3.820313e+00, 8.696155e+00, 9.720779e+00), 
               H_IC = c(0, 2.382904e-03, 2.348453e-02, 4.706867e-02, 1.265463e-01, 3.778414e-01, 1.658210e-01, 9.850088e-01, 1.745855e-01), 
               IC = c(0, 8.179626e-03, 8.058121e-02, 1.614645e-01, 4.336877e-01, 1.333769e+00, 2.629808e+00, 4.654115e+00, 9.273612e-01), 
               D = c(5.505607e-04, 8.362371e-04, 9.543505e-03, 2.727756e-02, 6.302218e-02, 2.526043e-01, 1.707657e+00, 3.840142e+00, 4.409166e+00), 
               R = c(1.482346e+02, 8.692593e+02, 2.120647e+03, 1.042637e+03, 1.088236e+03, 1.311541e+03, 1.375330e+03, 1.346433e+03, 1.068415e+03)
)

# Specify model parameters ------------------------------------------
# define contact/transmission matrix
april_2017 <- matrix(c(0.00005, 0.00003, 0.00004, 0.00004, 0.00003, 0.00002, 0.00002, 0.00001, 0.00001,
                       0.00003, 0.00043, 0.00020, 0.00008, 0.00010, 0.00007, 0.00005, 0.00003, 0.00002,
                       0.00003, 0.00018, 0.00063, 0.00021, 0.00017, 0.00016, 0.00011, 0.00006, 0.00010,
                       0.00003, 0.00008, 0.00022, 0.00020, 0.00013, 0.00010, 0.00008, 0.00005, 0.00007,
                       0.00002, 0.00009, 0.00017, 0.00012, 0.00018, 0.00011, 0.00008, 0.00005, 0.00007,
                       0.00001, 0.00006, 0.00014, 0.00008, 0.00010, 0.00014, 0.00011, 0.00007, 0.00009,
                       0.00001, 0.00005, 0.00011, 0.00008, 0.00009, 0.00013, 0.00028, 0.00018, 0.00019,
                       0.00001, 0.00003, 0.00009, 0.00006, 0.00008, 0.00011, 0.00025, 0.00056, 0.00055,
                       0.00001, 0.00006, 0.00026, 0.00018, 0.00018, 0.00030, 0.00050, 0.00106, 0.00211), nrow = 9)

# parameters must be in a named list
dt <- 1
params <- list(dt = dt,                                     # time step
               beta = 4.848224e-04/dt,                      # transmission probability
               beta1 = 0.14/dt,                             # seasonal forcing
               sigma = 0.5/dt,                              # rate from S -> E (1/sigma = latent period)
               gamma = 0.5/dt,                              # rate from I -> R (1/gamma = infectious period)
               epsilon = 0.005/dt,                          # case importation rate
               omega = 0.0038/dt,                           # waning rate (R -> S) 
               N = c(1796448.8, 2022910.4, 2217764.5,  2123500.3, 2277514.9, 2526594.7, 2105082.5, 1533156.7, 804612.3), # size of each age group (total population size = sum(N))
               h = c(0.0015,0.00007,0.00019,0.00069, 0.00129, 0.00281,0.00441, 0.0097, 0.010693)/dt,                     # rate from I -> H
               i1 = c(0, 0.02711, 0.04219, 0.04825, 0.07193, 0.08860, 0.10701, 0.08596,0.01535)/dt,                      # rate from H -> IC
               i2 = c(0.05551, 0.05551, 0.05551, 0.05551, 0.05551, 0.05314, 0.00769, 0.03673, 0.03558)/dt,               # rate from IC -> H_IC
               d = c(0.00027, 0.00062, 0.00139, 0.00313, 0.00357, 0.00573, 0.01514, 0.03271, 0.04443)/dt,                # rate from H -> D
               d_ic = c(0.00705, 0.00705, 0.00705, 0.00705, 0.00705, 0.009, 0.04632, 0.02247, 0.02342)/dt,               # rate from IC -> D 
               d_hic = c(0, 0, 0, 0, 0, 0.001, 0.004, 0.012, 0.029)/dt,                                                  # rate from H_IC -> D
               r = c(0.12634, 0.12603, 0.12535, 0.12381, 0.12342, 0.12151, 0.11316, 0.09759, 0.08722)/dt,                # rate from H -> R
               r_ic = c(0.0857, 0.0857, 0.0857, 0.0857, 0.0857, 0.0821, 0.0119, 0.0567, 0.0550)/dt,                      # rate from H_IC -> R
               contact_mat = april_2017,                    # transmission matrix
               calendar_start_date = as.Date("2020-01-01")) # calendar date of t=0 (this is used to calculate seasonality in transmission)

times <- seq(76, 119, by = 1)

# Run model -----------------------------------------------------------
seir_out <- lsoda(init_cond, times, age_struct_seir_ode_test, params)
out <- as.data.frame(seir_out)
# ---------------------------------------------------------------------

# Plot output ---------------------------------------------------------
# get number of people in each compartment
susceptibles <- rowSums(out[,c(paste0("S",1:9))])
exposed <- rowSums(out[,c(paste0("E",1:9))])
infected <- rowSums(out[,c(paste0("I",1:9))])
hospitalised <- rowSums(out[,c(paste0("H",1:9))])
ic <- rowSums(out[,c(paste0("IC",1:9))])
hosp_after_ic <- rowSums(out[,c(paste0("H_IC",1:9))])
deaths <- rowSums(out[,c(paste0("D",1:9))])
recovered <- rowSums(out[,c(paste0("R",1:9))]) 

# plot SEIR compartments
plot(susceptibles ~ times, type = "l", ylim = c(0, sum(params$N)))
abline(h = sum(params$N), lty = "dashed")
lines(recovered ~ times, type = "l", col = "blue") #, ylim = c(0,max(recovered))
lines(exposed ~ times, col = "green")
lines(infected ~ times, col = "red")
# plot severe disease compartments
plot(hospitalised ~ times, type = "l", col = "orange", ylim = c(min(ic),max(hospitalised)))
lines(ic ~ times, col = "pink", type = "l")
lines(hosp_after_ic ~ times, col = "purple")
lines(deaths ~ times, col = "grey")
# --------------------------------------------------------------------





