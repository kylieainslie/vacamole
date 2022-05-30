# Test script for debugging ----------------------------------------

# Load required packages/functions ---------------------------------
library(deSolve)
#library(reshape2)
#library(ggplot2)
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

# probabilities -----------------------------------------------------
dons_probs <- read_xlsx("inst/extdata/inputs/ProbabilitiesDelays_20210107.xlsx")
p_infection2admission <- dons_probs$P_infection2admission
p_admission2death <- dons_probs$P_admission2death
p_admission2IC <- dons_probs$P_admission2IC
p_IC2hospital <- dons_probs$P_IC2hospital
p_hospital2death <- c(rep(0, 5), 0.01, 0.04, 0.12, 0.29) # (after ICU)
p_reported_by_age <- c(0.29, 0.363, 0.381, 0.545, 0.645, 0.564, 0.365, 0.33, 0.409) # from Jantien
# p_reported_all <- 0.428 # from Jantien
# p_inf_by_age <- c(0.018, 0.115, 0.156, 0.118, 0.142, 0.199, 0.114, 0.062, 0.054 + 0.023)
# p_recovered <- c(
#   0.01120993, 0.09663659, 0.24141186, 0.11004723, 0.10677859, 0.11977255,
#   0.11904044, 0.11714503, 0.11347191
# )

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

# specify initial model parameters ---------------------------------
# parameters must be in a named list
dt <- 1
params <- list(dt = dt,
               N = n_vec,
               beta = 0.0004/dt, #4.848224e-04
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
breakpoint_sub <- df_breakpoints[1:10,]

n_bp <- length(breakpoint_sub$date)-1
mles <- matrix(rep(NA, 2*n_bp), nrow = n_bp)
colnames(mles) <- c("beta", "alpha")

out_mle <- list()
parameter_draws <- list()
beta_draws <- list()
daily_cases <- list()
init_update <- list()

breakpoints <- breakpoint_sub
fit_pars <- fit_params
case_data <- osiris1
contact_matrices <- cm_list
# begin loop over breakpoints --------------------------------
for (j in 1:n_bp) {
  
  print(paste("Fitting from", breakpoints$date[j], "to", breakpoints$date[j+1]))
  
  # set contact matrix for time window
  if (breakpoints$contact_matrix[j+1] == "april_2017"){contact_matrix <- contact_matrices$april_2017; print("april_2017")
  } else if (breakpoints$contact_matrix[j+1] == "april_2020"){contact_matrix <- contact_matrices$april_2020; print("april_2020")
  } else if (breakpoints$contact_matrix[j+1] == "june_2020"){contact_matrix <- contact_matrices$june_2020; print("june_2020")
  } else if (breakpoints$contact_matrix[j+1] == "september_2020"){contact_matrix <- contact_matrices$september_2020; print("septemeber_2020")
  } else if (breakpoints$contact_matrix[j+1] == "february_2021"){contact_matrix <- contact_matrices$february_2021; print("february_2021")
  } else if (breakpoints$contact_matrix[j+1] == "june_2021"){contact_matrix <- contact_matrices$june_2021; print("june_2021")
  } else {contact_matrix <- contact_matrices$november_2021} 
  
  params$contact_mat <- contact_matrix
  
  # set time sequence  
  times <- seq(breakpoints$time[j], breakpoints$time[j+1], by = 1)
  
  # update initial conditions based on last time window
  g <- mean(params$gamma)
  if (j == 1) {
    init_update <- init
    pars <- fit_pars$init_value
    S_diag <- diag(init_update[c(2:10)])
    rho <- as.numeric(eigs(S_diag %*% params$contact_mat, 1)$values)
  } else {
    end_states <- unlist(lapply(unname(out_mle[[j-1]]), tail,1))
    s_vec   <- end_states[c(paste0("S",1:9))]
    e_vec   <- end_states[c(paste0("E",1:9))]
    i_vec   <- end_states[c(paste0("I",1:9))]
    h_vec   <- end_states[c(paste0("H",1:9))]
    ic_vec  <- end_states[c(paste0("IC",1:9))]
    hic_vec <- end_states[c(paste0("H_IC",1:9))]
    d_vec   <- end_states[c(paste0("D",1:9))]
    r_vec   <- end_states[c(paste0("R",1:9))]
    r_vec1  <- end_states[c(paste0("R_1w",1:9))]
    r_vec2  <- end_states[c(paste0("R_2w",1:9))]
    r_vec3  <- n_vec - s_vec - e_vec - i_vec - h_vec - ic_vec - hic_vec - d_vec - r_vec - r_vec1 - r_vec2
    
    init_update <- c(t    = times[1],
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
    
    init_update <- c(t = times[1], end_states)
    
    # output error message if sum of all compartments is not equal to the total population size
    if(!isTRUE(all.equal(sum(init_update[-1]),sum(params$N)))){
      stop("Error: sum of compartments is not equal to population size")
    }
    if(any(end_states < 0)){
      stop("Error: Negative compartment values")
    }
    
    beta_est <- (mles[j-1,1]/g)*rho
    pars <- c(beta_est, fit_pars$init_value[2]) # mles[j-1,-1]
    S_diag <- diag(init_update[c(2:10)])
    rho <- as.numeric(eigs(S_diag %*% params$contact_mat, 1)$values)
  }
  
  # print(init_update)
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
  mles[j,1] <- (res$par[1] / rho) * g
  mles[j,2] <- res$par[2]
  
  print(mles[j,])
 # --------------------------------------------------
  # run for mle to get initial conditions for next timepoint
  params$beta <- mles[j,1]
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_update, times, age_struct_seir_ode_test, params, method = rk45, rtol = 1e-08, hmax = 0.02) # , method = rk45, rtol = 1e-08, hmax = 0.02
  out <- as.data.frame(seir_out)
  
  # store outputs
  out_mle[[j]] <- postprocess_age_struct_model_output2(out)
  daily_cases[[j]] <-  rowSums(params$sigma * out_mle[[j]]$E * params$p_report) #+ out_mle[[j]]$Ev_1d + out_mle[[j]]$Ev_2d + out_mle[[j]]$Ev_3d + out_mle[[j]]$Ev_4d + out_mle[[j]]$Ev_5d
  
  # plot for quick check of fit
  plot(times, case_data_sub$inc, pch = 16, col = "red", ylim = c(0, max(case_data_sub$inc,daily_cases[[j]])))
  lines(daily_cases[[j]]~times) # , type = "l", ylim = c(0, max(daily_cases[[j]]))
  
  
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
recovered1 <- list()
recovered2 <- list()
recovered3 <- list()
for (j in 1:n_bp){
  susceptibles[[j]] <- rowSums(out_mle[[j]]$S)
  exposed[[j]] <- rowSums(out_mle[[j]]$E)
  infected[[j]] <- rowSums(out_mle[[j]]$I)
  hospitalised[[j]] <- rowSums(out_mle[[j]]$H)
  ic[[j]] <- rowSums(out_mle[[j]]$IC)
  hosp_after_ic[[j]] <- rowSums(out_mle[[j]]$H_IC)
  deaths[[j]] <- rowSums(out_mle[[j]]$D)
  recovered[[j]] <- rowSums(out_mle[[j]]$R)
  recovered1[[j]] <- rowSums(out_mle[[j]]$R_1w)
  recovered2[[j]] <- rowSums(out_mle[[j]]$R_2w)
  recovered3[[j]] <- rowSums(out_mle[[j]]$R_3w)
}

# plot susceptibles
times_all <- seq(0,breakpoints$time[n_bp+1])
plot(unique(unlist(susceptibles)) ~ times_all, type = "l"#, ylim = c(0, sum(params$N))
     )
abline(h = sum(params$N), lty = "dashed")

# plot recovered
plot(unique(unlist(recovered)) ~ times_all, type = "l", col = "blue",ylim = c(0,max(unlist(recovered))))
lines(unique(unlist(recovered1)) ~ times_all, col = "blue", lty = "dashed")
lines(unique(unlist(recovered2)) ~ times_all, col = "blue", lty = "dotted")
lines(unique(unlist(recovered3)) ~ times_all, col = "blue", lty = "twodash")
lines(unique(unlist(exposed)) ~ times_all, col = "green")
lines(unique(unlist(infected)) ~ times_all, col = "red")

plot(unique(unlist(hospitalised)) ~ times_all, col = "orange", type = "l", ylim = c(min(unlist(ic)), max(unlist(hospitalised))))
plot(unique(unlist(ic)) ~ times_all, tpe = "l", col = "purple")
lines(unique(unlist(hosp_after_ic)) ~ times_all, col = "grey")

# --------------------------------------------------------------------
# Run forward simulations --------------------------------------------
# simulate forward with initial conditions that result in negative IC values
params$beta <- 0.0003445653
#times <- seq(0, 119, by = 1)

#rk45 <- rkMethod("rk45dp7")
seir_out <- ode(init_update, times, age_struct_seir_ode_test, params)  # , method = rk45, rtol = 1e-08, hmax = 0.02
seir_out1 <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output2(seir_out1)

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
points(times, case_data_sub$inc, pch = 16, col = "red") # , ylim = c(0, max(case_data_sub$inc))
# --------------------------------------------------------------------
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
# --------------------------------------------------------------------