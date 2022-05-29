# ------------------------------------------------------------------
# Fit model to data script
# ------------------------------------------------------------------
# Getting negative compartment values, so re-writing the entire fit 
# script. The negative compartment values are not reproduced in 
# minimal working example script with the same initial conditions/
# parameter values. Therefore, problem is likely somewhere in the 
# original script.
# ------------------------------------------------------------------

# Load required packages -------------------------------------------
library(deSolve)
library(lubridate)
library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(lubridate)
# ------------------------------------------------------------------
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
# -------------------------------------------------------------------

# Specify initial conditions ----------------------------------------
# define population size (by age group)
n_vec <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463,
           0.14514332, 0.12092904, 0.08807406, 0.04622194) * 17407585 # Dutch population size
empty_state <- c(rep(0, 9)) # vector of zeros
seed_age_group <- sample(1:9,1)
inf_seed_vec <- empty_state
inf_seed_vec[seed_age_group] <- 1

s_vec   <- n_vec - inf_seed_vec
e_vec   <- empty_state
i_vec   <- inf_seed_vec
h_vec   <- empty_state
ic_vec  <- empty_state
hic_vec <- empty_state
d_vec   <- empty_state
r_vec   <- empty_state
r_vec1  <- empty_state
r_vec2  <- empty_state
r_vec3  <- n_vec - s_vec - e_vec - i_vec - h_vec - ic_vec - hic_vec - d_vec - r_vec - r_vec1 - r_vec2

init <- c(t    = 0,
          S    = s_vec,
          E    = e_vec,
          I    = i_vec,
          H    = h_vec,
          IC   = ic_vec,
          H_IC = hic_vec,
          D    = d_vec,
          R    = r_vec,
          R_1w = r_vec1,
          R_2w = r_vec2,
          R_3w = r_vec3
)

# Specify model parameters ------------------------------------------
# define contact/transmission matrix --------------------------------
path <- "inst/extdata/inputs/contact_matrices/converted/"
april_2017     <- readRDS(paste0(path,"transmission_matrix_april_2017.rds"))
april_2020     <- readRDS(paste0(path,"transmission_matrix_april_2020.rds"))
june_2020      <- readRDS(paste0(path,"transmission_matrix_june_2020.rds"))
september_2020 <- readRDS(paste0(path,"transmission_matrix_september_2020.rds"))
february_2021  <- readRDS(paste0(path,"transmission_matrix_february_2021.rds"))
june_2021      <- readRDS(paste0(path,"transmission_matrix_june_2021.rds"))
november_2021  <- readRDS(paste0(path,"transmission_matrix_november_2021.rds"))

# put contact matrices into a list for input into fit_to_data_func()
contact_matrices <- list(
  april_2017 = april_2017,
  april_2020 = april_2020,
  june_2020 = june_2020,
  september_2020 = september_2020,
  february_2021 = february_2021,
  june_2021 = june_2021,
  november_2021 = november_2021
)

# probabilities -------------------------------------------------------
dons_probs <- read_xlsx("inst/extdata/inputs/ProbabilitiesDelays_20210107.xlsx")
p_infection2admission <- dons_probs$P_infection2admission
p_admission2death <- dons_probs$P_admission2death
p_admission2IC <- dons_probs$P_admission2IC
p_IC2hospital <- dons_probs$P_IC2hospital
p_hospital2death <- c(rep(0, 5), 0.01, 0.04, 0.12, 0.29) # (after ICU)
p_reported_by_age <- c(0.29, 0.363, 0.381, 0.545, 0.645, 0.564, 0.365, 0.33, 0.409) # from Jantien

# delays --------------------------------------------------------------
time_symptom2admission <- c(2.29, 5.51, 5.05, 5.66, 6.55, 5.88, 5.69, 5.09, 4.33) # assume same as infectious2admission
time_admission2discharge <- 7.9
time_admission2IC <- 2.28
time_IC2hospital <- 15.6
time_hospital2discharge <- 10.1 # (after ICU)
time_admission2death <- 7
time_IC2death <- 19
time_hospital2death <- 10 # (after ICU)

# define transition rates ---------------------------------------------
i2r    <- (1-p_infection2admission) / 2                    # I -> R
i2h    <- p_infection2admission / time_symptom2admission   # I -> H

h2ic   <- p_admission2IC / time_admission2IC               # H -> IC
h2d    <- p_admission2death / time_admission2death         # H -> D
h2r    <- (1 - (p_admission2IC + p_admission2death)) / time_admission2discharge
                                                           # H -> R

ic2hic <- p_IC2hospital / time_IC2hospital                 # IC -> H_IC
ic2d   <- (1 - p_IC2hospital) / time_IC2death              # IC -> D

hic2d  <- p_hospital2death / time_hospital2death           # H_IC -> D
hic2r  <- (1 - p_hospital2death) / time_hospital2discharge # H_IC -> R

dt <- 1

# model input parameters ---------------------------------------------
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

# -------------------------------------------------------------------
# Define likelihood function ----------------------------------------
# we're fitting the transmission probability (beta) and an
# over-dispersion parameter (alpha)
likelihood_func_test <- function(x, t, data, params, init) {
  # parameters to be estimated
  beta <- x[1]
  alpha <- x[2]
  
  # observed daily cases
  inc_obs <- data$inc
  
  # run model with current parameter values
  params$beta <- x[1]/10000
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init, t, age_struct_seir_ode_test, params, method = rk45, rtol = 1e-08, hmax = 0.02) # , rtol = 1e-08, hmax = 0.02
  out <- as.data.frame(seir_out)
  
  # modeled cases
  daily_cases <- rowSums(params$sigma * out[,c(paste0("E",1:9))] * params$p_report)
  daily_cases <- ifelse(daily_cases == 0, 0.0001, daily_cases) # prevent likelihood function function from being Inf
  
  # log-likelihood function
  # lik <- sum(dpois(x = inc_obs,lambda = incidence,log=TRUE))
  lik <- -sum(stats::dnbinom(x = inc_obs, mu = daily_cases, size = alpha, log = TRUE))
  
  print(x)
  print(lik)
  lik
}
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# define fit windows ------------------------------------------------
df_breakpoints <- read_csv2("inst/extdata/inputs/breakpoints_for_model_fit_v3.csv") %>%
  mutate(date = as.Date(date, format = "%d-%m-%Y"),
         time = as.numeric(date - date[1])) %>%
  select(date, time, variant, contact_matrix)

bp_for_fit <- df_breakpoints[1:5,]
n_bp <- dim(bp_for_fit)[1] - 1

# specify initial values and bounds for fitted parameters -----------
fit_params <- list(
  init_value = c(4, 1),
  lower_bound = c(0.5, 0.0001),
  upper_bound = c(Inf, Inf)
)

# create empty containers to store outputs from fit -----------------
init_cond <- list()     # store initial conditions for each window
out <- list()           # model output for each time point
times <- list()         # time points
cases <- list()         # daily cases to plot against real data
mles <- list()          # MLEs for each time window
# load case data ----------------------------------------------------
data_date <- "2022-03-12"
case_data <- readRDS(paste0("inst/extdata/data/case_data_upto_", data_date, ".rds"))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Fit model to data -------------------------------------------------
init_cond[[1]] <- init

# loop over time windows --------------------------------------------
for (j in 1:n_bp) {
  
  print(paste(paste0(j,")"),"Fitting from", bp_for_fit$date[j], "to", bp_for_fit$date[j+1]))
  
  # set contact matrix for time window ------------------------------
  if (bp_for_fit$contact_matrix[j+1] == "april_2017"){contact_matrix <- contact_matrices$april_2017 #; print("april_2017")
  } else if (bp_for_fit$contact_matrix[j+1] == "april_2020"){contact_matrix <- contact_matrices$april_2020 #; print("april_2020")
  } else if (bp_for_fit$contact_matrix[j+1] == "june_2020"){contact_matrix <- contact_matrices$june_2020 #; print("june_2020")
  } else if (bp_for_fit$contact_matrix[j+1] == "september_2020"){contact_matrix <- contact_matrices$september_2020 #; print("septemeber_2020")
  } else if (bp_for_fit$contact_matrix[j+1] == "february_2021"){contact_matrix <- contact_matrices$february_2021 #; print("february_2021")
  } else if (bp_for_fit$contact_matrix[j+1] == "june_2021"){contact_matrix <- contact_matrices$june_2021 #; print("june_2021")
  } else {contact_matrix <- contact_matrices$november_2021} 
  
  # change contact matrix in params list ----------------------------
  params$contact_mat <- contact_matrix
  
  # set time sequence -----------------------------------------------
  times[[j]] <- seq(bp_for_fit$time[j], bp_for_fit$time[j+1], by = 1)
  
  # subset data for time window -------------------------------------
  case_data_sub <- case_data[times[[j]] + 1, ]
  
  # run optimization procedure --------------------------------------
  res <- optim(par = fit_params$init_value, 
               fn = likelihood_func_test,
               method = "L-BFGS-B",
               lower = fit_params$lower_bound,
               upper = fit_params$upper_bound,
               t = times[[j]],
               data = case_data_sub,
               params = params,
               init = init_cond[[j]],
               hessian = TRUE
  )
  
  # store MLE --------------------------------------------------------
  mles[[j]] <- c(beta = res$par[1]/10000, alpha = res$par[2])
  print(mles[[j]])
  
  # Run model --------------------------------------------------------
  params$beta <- res$par[1]
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_cond[[j]], times[[j]], age_struct_seir_ode_test,  
                  params, method = rk45) #, rtol = 1e-08, hmax = 0.02
  
  # store outputs ----------------------------------------------------
  out[[j]] <- as.data.frame(seir_out) 
  cases[[j]] <-  rowSums(params$sigma * out[[j]][c(paste0("E",1:9))] * params$p_report)
  
 
  # plot for quick check of fit --------------------------------------
  plot(case_data_sub$inc ~ times[[j]], pch = 16, col = "red", 
       ylim = c(0, max(case_data_sub$inc,cases[[j]])))
  lines(cases[[j]] ~ times[[j]]) 
  
  # update initial conditions for next time window
  s_vec   <- tail(out[[j]][,c(paste0("S",1:9))],1)
  e_vec   <- tail(out[[j]][,c(paste0("E",1:9))],1)
  i_vec   <- tail(out[[j]][,c(paste0("I",1:9))],1)
  h_vec   <- tail(out[[j]][,c(paste0("H",1:9))],1)
  ic_vec  <- tail(out[[j]][,c(paste0("IC",1:9))],1)
  hic_vec <- tail(out[[j]][,c(paste0("H_IC",1:9))],1)
  d_vec   <- tail(out[[j]][,c(paste0("D",1:9))],1)
  r_vec   <- tail(out[[j]][,c(paste0("R",1:9))],1)
  r_vec1  <- tail(out[[j]][,c(paste0("R_1w",1:9))],1)
  r_vec2  <- tail(out[[j]][,c(paste0("R_2w",1:9))],1)
  r_vec3  <- params$N - s_vec - e_vec - i_vec - h_vec - ic_vec - hic_vec - d_vec - r_vec - r_vec1 - r_vec2
  
  init_cond[[j+1]] <- c(t    = tail(times[[j]],1),
                        S    = s_vec,
                        E    = e_vec,
                        I    = i_vec,
                        H    = h_vec,
                        IC   = ic_vec,
                        H_IC = hic_vec,
                        D    = d_vec,
                        R    = r_vec,
                        R_1w = r_vec1,
                        R_2w = r_vec2,
                        R_3w = r_vec3 )
  
  # output error message if negative compartment values
  if(any(init_cond[[j+1]] < 0)){
    stop("Error: Negative compartment values")
  }
  
} # end of for loop over breakpoints


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plot output -------------------------------------------------------
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