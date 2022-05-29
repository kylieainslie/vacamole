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
    R_1w <- c(R_1w1, R_1w2, R_1w3, R_1w4, R_1w5, R_1w6, R_1w7, R_1w8, R_1w9)
    R_2w <- c(R_2w1, R_2w2, R_2w3, R_2w4, R_2w5, R_2w6, R_2w7, R_2w8, R_2w9)
    R_3w <- c(R_3w1, R_3w2, R_3w3, R_3w4, R_3w5, R_3w6, R_3w7, R_3w8, R_3w9)
    
    # determine force of infection ----------------------------------
    # incorporate seasonality in transmission rate 
    calendar_day <- lubridate::yday(as.Date(times, origin = calendar_start_date))
    beta_t <- beta * (1 + beta1 * cos(2 * pi * calendar_day / 365.24)) 
    lambda <- beta_t * (contact_mat %*% I)
    # ---------------------------------------------------------------
    
    #################################################################
    # ODEs:
    dS    <- -lambda * S + (omega * 4 * R)
    dE    <- lambda * S - sigma * E + epsilon 
    dI    <- sigma * E - gamma * I - h * I
    dH    <- (h * I) - (i1 * H) - (d * H) - (r * H)
    dIC   <- (i1 * H) - (i2 * IC) - (d_ic * IC)
    dH_IC <- (i2 * IC) - (r_ic * H_IC) - (d_hic * H_IC)
    dD    <- (d * H) + (d_ic * IC) + (d_hic * H_IC)
    dR    <- (gamma * I) + (r * H) + (r_ic * H_IC) - (omega*4 * R) 
    dR_1w     <- (omega*4) * R - (omega*4) * R_1w
    dR_2w     <- (omega*4) * R_1w - (omega*4) * R_2w
    dR_3w     <- (omega*4) * R_2w - (omega*4) * R_3w
    #################################################################
    
    # output --------------------------------------------------------
    list(c(dt, dS, dE, dI, dH, dIC, dH_IC, dD, dR, dR_1w, dR_2w, dR_3w))
  })
}
# -------------------------------------------------------------------
times <- seq(83, 300, by = 1)

# Specify initial conditions ----------------------------------------
n_vec <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 0.14514332, 0.12092904, 0.08807406, 0.04622194) * 17407585 # Dutch population size
s_vec <- c(1796059.7, 2020649.2, 2212198.9, 2120736.6, 2274621.2, 2523059.8, 2101300.6, 1529360.5, 801594.8)
e_vec <- c(75.86727,  440.70015, 1084.04776,  538.51798,  563.87157,  688.51887, 736.66532, 739.16223, 587.27497)
i_vec <- c(60.18361, 349.56693, 860.32740, 427.98032, 449.32464, 550.90366, 592.06253, 600.81061, 475.48768)
h_vec <- c(0.36178738, 0.08791759,  0.57106796,  1.01938830,  1.89264935,  4.86809927, 7.86665679, 18.30326062, 18.95899090)
ic_vec <- c(0, 0.005122535, 0.051756889, 0.105509990, 0.291407610, 0.912051533, 1.475548634, 3.098444148, 0.573477983)
hic_vec <- c(0, 0.01020844, 0.10313389, 0.21039705, 0.58155193, 1.83324077, 3.25508254, 6.48969267, 1.20301073)
d_vec <- c(0.0008745288, 0.0011873604, 0.0141734602, 0.0427083027, 0.0997896145, 0.4174739580, 2.7315543074, 7.0906894849, 7.8916046775 )
r_vec <- c(222.5397, 1294.9318, 3187.3854, 1581.0716, 1653.2164, 2013.4445, 2146.5559, 2132.0079, 1695.9995)
r_vec1 <- c(26.63577, 155.03183, 381.65575, 189.25263, 197.88501, 241.30536, 256.80275, 255.00317, 202.96027)
r_vec2 <- c(3.206318, 18.669604, 45.948723, 22.776022, 23.815222, 29.211642, 30.857473, 30.600270, 24.357005)
r_vec3 <- n_vec - s_vec - e_vec - i_vec - h_vec - hic_vec - ic_vec - d_vec - r_vec - r_vec1 - r_vec2

init_cond <- c(t = times[1],
               S = s_vec,
               E = e_vec,
               I = i_vec,
               H = h_vec,
               IC = ic_vec,
               H_IC = hic_vec,
               D = d_vec,
               R = r_vec,
               R_1w = r_vec1,
               R_2w = r_vec2,
               R_3w = r_vec3
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

# probabilities -----------------------------------------------------
dons_probs <- read_xlsx("inst/extdata/inputs/ProbabilitiesDelays_20210107.xlsx")
p_infection2admission <- dons_probs$P_infection2admission
p_admission2death <- dons_probs$P_admission2death
p_admission2IC <- dons_probs$P_admission2IC
p_IC2hospital <- dons_probs$P_IC2hospital
p_hospital2death <- c(rep(0, 5), 0.01, 0.04, 0.12, 0.29) # (after ICU)
p_reported_by_age <- c(0.29, 0.363, 0.381, 0.545, 0.645, 0.564, 0.365, 0.33, 0.409) # from Jantien

# delays ------------------------------------------------------------
time_symptom2admission <- c(2.29, 5.51, 5.05, 5.66, 6.55, 5.88, 5.69, 5.09, 4.33) # assume same as infectious2admission
time_admission2discharge <- 7.9
time_admission2IC <- 2.28
time_IC2hospital <- 15.6
time_hospital2discharge <- 10.1 # (after ICU)
time_admission2death <- 7
time_IC2death <- 19
time_hospital2death <- 10 # (after ICU)

# define transition rates -------------------------------------------
i2r    <- (1-p_infection2admission) / 2                   # I -> R
i2h    <- p_infection2admission / time_symptom2admission      # I -> H

h2ic   <- p_admission2IC / time_admission2IC                 # H -> IC
h2d    <- p_admission2death / time_admission2death            # H -> D
h2r    <- (1 - (p_admission2IC + p_admission2death)) / time_admission2discharge
# H -> R

ic2hic <- p_IC2hospital / time_IC2hospital                   # IC -> H_IC
ic2d   <- (1 - p_IC2hospital) / time_IC2death                       # IC -> D

hic2d  <- p_hospital2death / time_hospital2death          # H_IC -> D
hic2r  <- (1 - p_hospital2death) / time_hospital2discharge # H_IC -> R

dt <- 1
# parameters must be in a named list
params <- list(dt = dt,
               N = n_vec,
               beta = 0.0004734092/dt,
               beta1 = 0.14/dt,
               sigma = 0.5/dt,
               gamma = i2r/dt,
               h = i2h/dt,
               i1 = h2ic/dt,
               d = h2d/dt,
               r = h2r/dt,
               i2 = ic2hic/dt,
               d_ic = ic2d/dt,
               d_hic = hic2d/dt,
               r_ic = hic2r/dt,
               epsilon = 0.00/dt,
               omega = 0.0038/dt,
               p_report = p_reported_by_age,
               contact_mat = april_2017,
               calendar_start_date = as.Date("2020-01-01")
)

# Run model -----------------------------------------------------------
rk45 <- rkMethod("rk45dp7")
seir_out <- ode(init_cond, times, age_struct_seir_ode_test, params, method = rk45) #, rtol = 1e-08, hmax = 0.02
out <- as.data.frame(seir_out) 
# ---------------------------------------------------------------------

# Plot output ---------------------------------------------------------
# get number of people in each compartment
susceptibles  <- rowSums(out[,c(paste0("S",1:9))])
exposed       <- rowSums(out[,c(paste0("E",1:9))])
infected      <- rowSums(out[,c(paste0("I",1:9))])
hospitalised  <- rowSums(out[,c(paste0("H",1:9))])
ic            <- rowSums(out[,c(paste0("IC",1:9))])
hosp_after_ic <- rowSums(out[,c(paste0("H_IC",1:9))])
deaths        <- rowSums(out[,c(paste0("D",1:9))])
recovered     <- rowSums(out[,c(paste0("R",1:9))]) 
recovered1    <- rowSums(out[,c(paste0("R_1w",1:9))]) 
recovered2    <- rowSums(out[,c(paste0("R_2w",1:9))]) 
recovered3    <- rowSums(out[,c(paste0("R_3w",1:9))]) 

# plot SEIR compartments
plot(susceptibles ~ times, type = "l", ylim = c(0, sum(params$N)))
abline(h = sum(params$N), lty = "dashed")
lines(recovered ~ times, type = "l", col = "blue") #, ylim = c(0,max(recovered))
lines(recovered1 ~ times, col = "blue", lty = "dashed")
lines(recovered2 ~ times, col = "blue", lty = "dotted")
lines(recovered3 ~ times, col = "blue", lty = "twodash")
lines(exposed ~ times, col = "green")
lines(infected ~ times, col = "red")
# plot severe disease compartments
plot(hospitalised ~ times, type = "l", col = "orange", ylim = c(min(ic),max(hospitalised)))
lines(ic ~ times, col = "pink", type = "l")
lines(hosp_after_ic ~ times, col = "purple")
lines(deaths ~ times, col = "grey")
# --------------------------------------------------------------------




