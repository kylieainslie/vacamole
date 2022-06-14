# Script for running scenarios for the Eupropean Scenario Hub
# URL: https://github.com/covid19-forecast-hub-europe/covid19-scenario-hub-europe#readme

# preamble ---------------------------------------------------------
# This script will load necessary packages, data sources, fit the model to data,
# and then run scenarios.
# All scenarios will be using Dutch data
# TODO: generalize to other European countries
# TODO: break up code into work chunks
# ------------------------------------------------------------------

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
data_date <- "2022-05-22"
osiris1 <- readRDS(paste0("inst/extdata/data/case_data_upto_", data_date, ".rds"))

# read in transition rates -----------------------------------------
transition_rates <- readRDS("inst/extdata/inputs/transition_rates.rds")

# define population size (by age group)
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 
              0.13083463,0.14514332, 0.12092904, 0.08807406, 
              0.04622194)
n <- 17407585 # Dutch population size
n_vec <- n * age_dist

# determine waning rate from Erlang distribution --------------------
# We want the rate that corresponds to a 60% reduction in immunity after 
#   - 3 months (92 days) or
#   - 8 months (244 days)

# we need to solve the following equation for lambda (waning rate)
# tau = time since recovery
# p = probability still immune
Fk <- function(lambda, tau, p){
  exp(-tau * lambda) * (6 + (6 * tau * lambda) + (3 * tau^2 * lambda^2) 
                        + (tau^3 * lambda^3)) - (p * 6)
}

wane_3months <- uniroot(Fk, c(0,1), tau = 92, p = 0.6)$root
wane_8months <- uniroot(Fk, c(0,1), tau = 244, p = 0.6)$root
# 50% reduction after 6 months (used for model fits)
wane_6months <- uniroot(Fk, c(0,1), tau = 182, p = 0.5)$root
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
# ve estimates ------------------------------------------------------
# list containing the following named lists:
# delays, ve_inf, ve_hosp, ve_trans
# each named list has the following named elements:
# pfizer, moderna, astrazeneca, jansen
ve_params <- readRDS("inst/extdata/inputs/ve_params.rds")

# vaccination schedule ----------------------------------------------
# read in vaccination schedule
vac_schedule <- read_csv("inst/extdata/inputs/vac_schedule_real_w_4th_and_5th_dose.csv") %>%
  select(-X1)

# convert vaccination schedule for input into model
vac_rates_wt <- convert_vac_schedule2(
  vac_schedule = vac_schedule,
  delay = ve_params$delays,
  ve = ve_params$ve_inf$wildtype,
  hosp_multiplier = ve_params$ve_hosp$wildtype,
  ve_trans = ve_params$ve_trans$wildtype,
  wane = FALSE)

vac_rates_alpha <- convert_vac_schedule2(
  vac_schedule = vac_schedule,
  delay = ve_params$delays,
  ve = ve_params$ve_inf$alpha,
  hosp_multiplier = ve_params$ve_hosp$alpha,
  ve_trans = ve_params$ve_trans$alpha,
  wane = FALSE)

vac_rates_delta <- convert_vac_schedule2(
  vac_schedule = vac_schedule,
  delay = ve_params$delays,
  ve = ve_params$ve_inf$delta,
  hosp_multiplier = ve_params$ve_hosp$delta,
  ve_trans = ve_params$ve_trans$delta,
  wane = FALSE)

vac_rates_omicron <- convert_vac_schedule2(
  vac_schedule = vac_schedule,
  delay = ve_params$delays,
  ve = ve_params$ve_inf$omicron,
  hosp_multiplier = ve_params$ve_hosp$omicron,
  ve_trans = ve_params$ve_trans$omicron,
  wane = FALSE)

# make into a list for input into fit_to_data_func()
vac_rates_list <- list(
  wildtype = vac_rates_wt,
  alpha = vac_rates_alpha,
  delta = vac_rates_delta,
  omicron = vac_rates_omicron
)
# specify initial model parameters ---------------------------------
# parameters must be in a named list
params <- list(beta = 0.0004,
               beta1 = 0.14,
               gamma = 0.5,
               sigma = 0.5,
               epsilon = 0.0,
               omega = 0.0038,
               N = n_vec,
               h = transition_rates$h,
               i1 = transition_rates$i1,
               i2 = transition_rates$i2,
               d = transition_rates$d, 
               d_ic = transition_rates$d_ic,
               d_hic = transition_rates$d_hic,
               r = transition_rates$r,
               r_ic = transition_rates$r_ic,
               p_report = 1/3,
               c_start = april_2017,
               keep_cm_fixed = TRUE,
               vac_inputs = vac_rates_wt,
               use_cases = TRUE,  
               no_vac = FALSE,
               t_calendar_start = yday(as.Date("2020-01-01")), 
               beta_change = 0.0001,
               t_beta_change = 165
              )

# Specify initial conditions --------------------------------------
empty_state <- c(rep(0, 9)) # vector of zeros

init <- c(
  t = 0,
  S = c(n_vec[1:4], n_vec[5]-1, n_vec[6:9]),
  Shold_1d = empty_state,
  Sv_1d = empty_state,
  Shold_2d = empty_state,
  Sv_2d = empty_state,
  Shold_3d = empty_state,
  Sv_3d = empty_state,
  Shold_4d = empty_state,
  Sv_4d = empty_state,
  Shold_5d = empty_state,
  Sv_5d = empty_state,
  E = empty_state,
  Ev_1d = empty_state,
  Ev_2d = empty_state,
  Ev_3d = empty_state,
  Ev_4d = empty_state,
  Ev_5d = empty_state,
  I = c(rep(0,4),1,rep(0,4)),
  Iv_1d = empty_state,
  Iv_2d = empty_state,
  Iv_3d = empty_state,
  Iv_4d = empty_state,
  Iv_5d = empty_state,
  H = empty_state,
  Hv_1d = empty_state,
  Hv_2d = empty_state,
  Hv_3d = empty_state,
  Hv_4d = empty_state,
  Hv_5d = empty_state,
  IC = empty_state,
  ICv_1d = empty_state,
  ICv_2d = empty_state,
  ICv_3d = empty_state,
  ICv_4d = empty_state,
  ICv_5d = empty_state,
  H_IC = empty_state,
  H_ICv_1d = empty_state,
  H_ICv_2d = empty_state,
  H_ICv_3d = empty_state,
  H_ICv_4d = empty_state,
  H_ICv_5d = empty_state,
  D = empty_state,
  R = empty_state,
  Rv_1d = empty_state,
  Rv_2d = empty_state,
  Rv_3d = empty_state,
  Rv_4d = empty_state,
  Rv_5d = empty_state,
  R_1w = empty_state, 
  Rv_1d_1w = empty_state, 
  Rv_2d_1w = empty_state, 
  Rv_3d_1w = empty_state, 
  Rv_4d_1w = empty_state, 
  Rv_5d_1w = empty_state,
  R_2w = empty_state, 
  Rv_1d_2w = empty_state, 
  Rv_2d_2w = empty_state, 
  Rv_3d_2w = empty_state, 
  Rv_4d_2w = empty_state, 
  Rv_5d_2w = empty_state,
  R_3w = empty_state, 
  Rv_1d_3w = empty_state, 
  Rv_2d_3w = empty_state, 
  Rv_3d_3w = empty_state, 
  Rv_4d_3w = empty_state, 
  Rv_5d_3w = empty_state
)

# Model fit ---------------------------------------------------------
# read in csv with breakpoints
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
# run fit procedure
breakpoint_sub <- df_breakpoints[1:7,]

fits <- fit_to_data_func(breakpoints = breakpoint_sub, params = params, 
                         init = init, fit_pars = fit_params,
                         case_data = osiris1, contact_matrices = cm_list,
                         vac_info = vac_rates_list,
                         save_output_to_file = FALSE, path_out = NULL)

# plot S and R as a check
s_comp <- unique(unlist(susceptibles))
r_comp <- unique(unlist(recovered))
r_comp1 <-unique(unlist(recovered1))
r_comp2 <- unique(unlist(recovered2))
t_all <- seq(0,breakpoint_sub$time[length(breakpoint_sub$time)], by = 1)

plot(s_comp~t_all, type = "l", ylim = c(0, s_comp[1]))
plot(r_comp ~ t_all, type = "l", col = "green")
lines(r_comp1[1:190]~t_all[1:190], col = "blue")
lines(r_comp2[1:190]~t_all[1:190], col = "red")
abline(h = sum(params$N), lty = "dashed")
abline(v = 152, lty = "dotted")
# Run forward simulations --------------------------------------------
times <- seq(0,250, by = 1)

seir_out <- lsoda(init, times, age_struct_seir_ode2, params)
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output2(seir_out)

susceptibles <- rowSums(out$S)
exposed <- rowSums(out$E)
infected <- rowSums(out$I)
hospitalised <- rowSums(out$H)
ic <- rowSums(out$IC)
hosp_after_ic <- rowSums(out$H_IC)
deaths <- rowSums(out$D)
recovered <- rowSums(out$R) + rowSums(out$R_1w) + rowSums(out$R_2w) + rowSums(out$R_3w) 
#cases <- params$sigma * rowSums(out$E + out$Ev_1d + out$Ev_2d + out$Ev_3d + out$Ev_4d + out$Ev_5d) * params$p_report

plot(susceptibles ~ times, type = "l") #, ylim = c(0, sum(params$N))
abline(h = sum(params$N), lty = "dashed")

plot(recovered ~ times, type = "l", col = "blue", ylim = c(0,max(recovered)))
lines(exposed ~ times, col = "green")
lines(infected ~ times, col = "red")
# lines(hospitalised ~ times, col = "orange")
# lines(ic ~ times, col = "pink")
# lines(hosp_after_ic ~ times, col = "purple")

plot((susceptibles + exposed + infected + recovered +
         hospitalised + ic + hosp_after_ic) ~ times, type = "l", lty = "dashed")

