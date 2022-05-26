# Test script for debugging ----------------------------------------

# Load required packages/functions ---------------------------------
library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
library(rARPACK)
library(readr)
library(lubridate)

devtools::load_all() # load vacamole
library(vacamole)
# -------------------------------------------------------------------

# Load data ---------------------------------------------------------
# if off the server, read in from inst/extdata/data
data_date <- "2022-03-12"
osiris1 <- readRDS(paste0("inst/extdata/data/case_data_upto_", data_date, ".rds"))

# read in transition rates -----------------------------------------
transition_rates <- readRDS("inst/extdata/inputs/transition_rates.rds")

# define population size (by age group)
n_vec <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463,
           0.14514332, 0.12092904, 0.08807406, 0.04622194) * 17407585 # Dutch population size

# contact matrices --------------------------------------------------
#path <- "/rivm/s/ainsliek/data/contact_matrices/converted/"
path <- "inst/extdata/inputs/contact_matrices/converted/"
april_2017     <- readRDS(paste0(path,"transmission_matrix_april_2017.rds"))
april_2020     <- readRDS(paste0(path,"transmission_matrix_april_2020.rds"))
june_2020      <- readRDS(paste0(path,"transmission_matrix_june_2020.rds"))
september_2020 <- readRDS(paste0(path,"transmission_matrix_september_2020.rds"))
february_2021  <- readRDS(paste0(path,"transmission_matrix_february_2021.rds"))
june_2021      <- readRDS(paste0(path,"transmission_matrix_june_2021.rds"))
november_2021  <- readRDS(paste0(path,"transmission_matrix_november_2021.rds"))

# put contact matrices into a list for input into fit_to_data_func()
cm_list <- list(
  april_2017 = april_2017,
  april_2020 = april_2020,
  june_2020 = june_2020,
  september_2020 = september_2020,
  february_2021 = february_2021,
  june_2021 = june_2021,
  november_2021 = november_2021
)

# Specify initial conditions --------------------------------------
empty_state <- c(rep(0, 9)) # vector of zeros

init <- c(
  t = 0,
  S = c(n_vec[1:4], n_vec[5]-1, n_vec[6:9]),
  E = empty_state,
  I = c(rep(0,4),1,rep(0,4)),
  H = empty_state,
  IC = empty_state,
  H_IC = empty_state,
  D = empty_state,
  R = empty_state #,
  # R_1w = empty_state,
  # R_2w = empty_state,
  # R_3w = empty_state
)


# specify initial model parameters ---------------------------------
# parameters must be in a named list
dt <- 1
params <- list(dt = dt,
               beta = 0.0004/dt,
               beta1 = 0.14/dt,
               gamma = 0.5/dt,
               sigma = 0.5/dt,
               epsilon = 0.005/dt,
               omega = 0.0038/dt,
               N = n_vec,
               h = transition_rates$h/dt,
               i1 = transition_rates$i1/dt,
               i2 = c(rep(0,9)), #transition_rates$i2/dt,
               d = transition_rates$d/dt,
               d_ic = c(rep(0,9)), #transition_rates$d_ic/dt,
               d_hic = transition_rates$d_hic/dt,
               r = transition_rates$r/dt,
               r_ic = transition_rates$r_ic/dt,
               p_report = 1/3,
               c_start = april_2017,
               #keep_cm_fixed = TRUE,
               #vac_inputs = vac_rates_wt,
               #use_cases = TRUE,  
               #no_vac = FALSE,
               calendar_start_date = as.Date("2020-01-01"), 
               beta_change = NULL,
               t_beta_change = NULL
)
# --------------------------------------------------------------------
# Run forward simulations --------------------------------------------
times <- seq(0,500, by = 1)

seir_out <- lsoda(init, times, age_struct_seir_ode_test, params)
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output2(seir_out)

# get number of people in each compartment
susceptibles <- rowSums(out$S)
exposed <- rowSums(out$E)
infected <- rowSums(out$I)
hospitalised <- rowSums(out$H)
ic <- rowSums(out$IC)
hosp_after_ic <- rowSums(out$H_IC)
deaths <- rowSums(out$D)
recovered <- rowSums(out$R) #+ rowSums(out$R_1w) + rowSums(out$R_2w) + rowSums(out$R_3w) 

cases <- params$sigma * rowSums(out$E) * params$p_report
plot(cases~times, type = "l")
# --------------------------------------------------------------------
# plot 
plot(susceptibles ~ times, type = "l", ylim = c(0, sum(params$N)))
abline(h = sum(params$N), lty = "dashed")
lines(recovered ~ times, type = "l", col = "blue") #, ylim = c(0,max(recovered))
lines(exposed ~ times, col = "green")
lines(infected ~ times, col = "red")
lines(hospitalised ~ times, col = "orange")
lines(ic ~ times, col = "pink")
lines(hosp_after_ic ~ times, col = "purple")
lines(deaths ~ times, col = "grey")
# --------------------------------------------------------------------
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# try model fits to first wave ----------------------------------------
df_breakpoints <- read_csv2("inst/extdata/inputs/breakpoints_for_model_fit_v3.csv") %>%
  mutate(date = as.Date(date, format = "%d-%m-%Y"),
         time = as.numeric(date - date[1])) %>%
  select(date, time, variant, contact_matrix)

# specify initial value and bounds for fitted parameters
fit_params <- list(
  init_value = c(2.3, 1),
  lower_bound = c(0.0001, 0.0001),
  upper_bound = c(Inf, Inf)
)

# run fit procedure (not using function wrapper) 
breakpoint_sub <- df_breakpoints[1:18,]

n_bp <- length(breakpoint_sub$date)-1
mles <- matrix(rep(NA, 2*n_bp), nrow = n_bp)
colnames(mles) <- c("beta", "alpha")

out_mle <- list()
parameter_draws <- list()
beta_draws <- list()
daily_cases <- list()
# susceptibles <- list()
# exposed <- list()
# infected <- list()
# recovered <- list()
# recovered1 <- list()
# recovered2 <- list()
# recovered3 <- list()

breakpoints <- breakpoint_sub
fit_pars <- fit_params
case_data <- osiris1
contact_matrices <- cm_list
# begin loop over breakpoints --------------------------------
for (j in 1:n_bp) {
  
  print(paste("Fitting from", breakpoints$date[j], "to", breakpoints$date[j+1]))
  
  # set contact matrix for time window
  if (breakpoints$contact_matrix[j+1] == "april_2017"){contact_matrix <- contact_matrices$april_2017
  } else if (breakpoints$contact_matrix[j+1] == "april_2020"){contact_matrix <- contact_matrices$april_2020
  } else if (breakpoints$contact_matrix[j+1] == "june_2020"){contact_matrix <- contact_matrices$june_2020
  } else if (breakpoints$contact_matrix[j+1] == "september_2020"){contact_matrix <- contact_matrices$september_2020
  } else if (breakpoints$contact_matrix[j+1] == "february_2021"){contact_matrix <- contact_matrices$february_2021
  } else if (breakpoints$contact_matrix[j+1] == "june_2021"){contact_matrix <- contact_matrices$june_2021
  } else {contact_matrix <- contact_matrices$november_2021} 
  
  params$c_start <- contact_matrix
  
  # set time sequence  
  times <- seq(breakpoints$time[j], breakpoints$time[j+1], by = 1)
  
  # update initial conditions based on last time window
  if (j == 1) {
    init_update <- init
    pars <- fit_pars$init_value
    S_diag <- diag(init_update[c(2:10)])
    rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  } else {
    init_update <- c(t = times[1], unlist(lapply(unname(out_mle[[j-1]]), tail,1)))
    beta_est <- (mles[j-1,1]/params$gamma)*rho
    pars <- c(beta_est, fit_pars$init_value[2]) # mles[j-1,-1]
    S_diag <- diag(init_update[c(2:10)])
    rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  }
  
  # subset data for time window
  case_data_sub <- case_data[times + 1, ]
  
  res <- optim(par = pars, 
               fn = likelihood_func_test,
               method = "L-BFGS-B",
               lower = fit_pars$lower_bound,
               upper = fit_pars$upper_bound,
               t = times,
               data = case_data_sub,
               params = params,
               init = init_update,
               model_func = age_struct_seir_ode_test,
               hessian = TRUE
  )
  
  # store MLE
  mles[j,1] <- (res$par[1] / rho) * params$gamma
  mles[j,2] <- res$par[2]
  
  print(mles[j,])
 # --------------------------------------------------
  # run for mle to get initial conditions for next timepoint
  params$beta <- mles[j,1]
  
  seir_out <- lsoda(init_update, times, age_struct_seir_ode_test, params)
  seir_out <- as.data.frame(seir_out)
  
  # store outputs
  out_mle[[j]] <- postprocess_age_struct_model_output2(seir_out)
  daily_cases[[j]] <- params$sigma * rowSums(out_mle[[j]]$E) * params$p_report #+ out_mle[[j]]$Ev_1d + out_mle[[j]]$Ev_2d + out_mle[[j]]$Ev_3d + out_mle[[j]]$Ev_4d + out_mle[[j]]$Ev_5d
  # susceptibles[[j]] <- rowSums(out_mle[[j]]$S)
  # exposed[[j]] <- rowSums(out_mle[[j]]$E)
  # infected[[j]] <- rowSums(out_mle[[j]]$I)
  # recovered[[j]] <- rowSums(out_mle[[j]]$R) 
  # recovered1[[j]] <- rowSums(out_mle[[j]]$R_1w)
  # recovered2[[j]] <- rowSums(out_mle[[j]]$R_2w)
  # recovered3[[j]] <- rowSums(out_mle[[j]]$R_3w)
  
  # plot for quick check of fit
  plot(daily_cases[[j]]~times, type = "l")
  points(times, case_data_sub$inc, pch = 16, col = "red")
  
} # end of for loop over breakpoints


# plot states ----------------------------------------
cases <- list()
susceptibles <- list()
exposed <- list()
infected <- list()
hospitalised <- list()
ic <- list()
hosp_after_ic <- list()
recovered <- list()
deaths <- list()
# recovered1 <- list()
# recovered2 <- list()
# recovered3 <- list()
for (j in 1:n_bp){
  susceptibles[[j]] <- rowSums(out_mle[[j]]$S)
  exposed[[j]] <- rowSums(out_mle[[j]]$E)
  infected[[j]] <- rowSums(out_mle[[j]]$I)
  hospitalised[[j]] <- rowSums(out_mle[[j]]$H)
  ic[[j]] <- rowSums(out_mle[[j]]$IC)
  hosp_after_ic[[j]] <- rowSums(out_mle[[j]]$H_IC)
  deaths[[j]] <- rowSums(out_mle[[j]]$D)
  recovered[[j]] <- rowSums(out_mle[[j]]$R)
  # recovered1[[j]] <- rowSums(out_mle[[j]]$R_1w)
  # recovered2[[j]] <- rowSums(out_mle[[j]]$R_2w)
  # recovered3[[j]] <- rowSums(out_mle[[j]]$R_3w)
}

# plot susceptibles
times_all <- seq(0,breakpoints$time[n_bp+1])
plot(unique(unlist(susceptibles)) ~ times_all, type = "l"#, ylim = c(0, sum(params$N))
     )
abline(h = sum(params$N), lty = "dashed")

# plot recovered
plot(unique(unlist(recovered)) ~ times_all, type = "l", col = "blue", ylim = c(0,max(unlist(recovered))))
lines(unique(unlist(exposed)) ~ times_all, col = "green")
lines(unique(unlist(infected)) ~ times_all, col = "red")
plot(unique(unlist(hospitalised)) ~ times_all, col = "orange")
plot(unique(unlist(ic)) ~ times_all, tpe = "l", col = "purple")
lines(unique(unlist(hosp_after_ic)) ~ times_all, col = "grey")
