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
# osiris case data --------------------------------------------------
# this must be run on the RIVM R servers
path <- "/rivm/r/COVID-19/Surveillance/Data/OSIRIS/Geschoond/"
file <- list.files(path, pattern = ".rds")
if (identical(file, character(0))) {
  path <- paste0(path,"Previous/")
  file <- list.files(path, pattern = ".rds") %>%
    max()
}

osiris <- readRDS(paste0(path,file)) # read in file from path

osiris_tally <- osiris %>%           # aggregate for number of cases per day
                                     # this removes any identifiable data
  select(OSIRISNR, INFECTIEZIEKTE, ZIE1eZiekteDt, Land) %>%
  filter(Land == "Nederland",
         INFECTIEZIEKTE %in% c("NCOV", "Weak Positive", 
                               "Antitgen Pos. + Symptoms", 
                               "PCR Positief", "Antigen Positief")) %>%
  select(-Land) %>%
  rename(date = ZIE1eZiekteDt) %>%
  group_by(date) %>%
  summarise(inc = n()) %>%
  filter(!is.na(date)) %>%
  complete(date = seq.Date(min(date), max(date), by="day"), fill = list(inc = 0))

cutoff_date <- as.Date("2022-03-12")

osiris1 <- osiris_tally %>%
  filter(date <= cutoff_date)

# if off the server, read in from inst/extdata/data
data_date <- "2022-03-12"
osiris1 <- readRDS(paste0("inst/extdata/data/case_data_upto_", data_date, ".rds"))

# read in transition rates -----------------------------------------
transition_rates <- readRDS("inst/extdata/inputs/transition_rates.rds")

# define population size (by age group)
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 
              0.13083463,0.14514332, 0.12092904, 0.08807406, 
              0.04622194)
n <- 17407585 # Dutch population size
n_vec <- n * age_dist

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
vac_schedule <- read_csv("inst/extdata/inputs/vac_schedule_real_w_booster.csv") %>%
  rename_with(~ gsub("B", "d3", .x, fixed = TRUE)) %>%
  select(-starts_with("...")) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

# create "empty" values from 1/1/2020 until start of vac sched (1/4/2021)
n_cols <- dim(vac_schedule)[2]-1 #exclude date column
n_rows <- vac_schedule$date[1] - osiris1$date[1]
empty_mat <- matrix(rep(0, n_cols * n_rows), nrow = n_rows)
dates <- seq.Date(osiris1$date[1], vac_schedule$date[1]-1, by = "day")
my_df <- data.frame(date = dates, empty_mat)
names(my_df) <- names(vac_schedule)
vac_sched <- bind_rows(my_df, vac_schedule)

vac_rates_wt <- convert_vac_schedule2(
  vac_schedule = vac_sched,
  delay = ve_params$delays,
  ve = ve_params$ve_inf$wildtype,
  hosp_multiplier = ve_params$ve_hosp$wildtype,
  ve_trans = ve_params$ve_trans$wildtype,
  wane = FALSE)

vac_rates_alpha <- convert_vac_schedule2(
  vac_schedule = vac_sched,
  delay = ve_params$delays,
  ve = ve_params$ve_inf$alpha,
  hosp_multiplier = ve_params$ve_hosp$alpha,
  ve_trans = ve_params$ve_trans$alpha,
  wane = FALSE)

vac_rates_delta <- convert_vac_schedule2(
  vac_schedule = vac_sched,
  delay = ve_params$delays,
  ve = ve_params$ve_inf$delta,
  hosp_multiplier = ve_params$ve_hosp$delta,
  ve_trans = ve_params$ve_trans$delta,
  wane = FALSE)

vac_rates_omicron <- convert_vac_schedule2(
  vac_schedule = vac_sched,
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
               epsilon = 0.01,
               omega = 0.0016,
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
               no_vac = TRUE,
               t_calendar_start = yday(as.Date("2020-01-01")), 
               beta_change = NULL 
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
  E = empty_state,
  Ev_1d = empty_state,
  Ev_2d = empty_state,
  Ev_3d = empty_state,
  I = c(rep(0,4),1,rep(0,4)),
  Iv_1d = empty_state,
  Iv_2d = empty_state,
  Iv_3d = empty_state,
  H = empty_state,
  Hv_1d = empty_state,
  Hv_2d = empty_state,
  Hv_3d = empty_state,
  IC = empty_state,
  ICv_1d = empty_state,
  ICv_2d = empty_state,
  ICv_3d = empty_state,
  H_IC = empty_state,
  H_ICv_1d = empty_state,
  H_ICv_2d = empty_state,
  H_ICv_3d = empty_state,
  D = empty_state,
  R = empty_state,
  Rv_1d = empty_state,
  Rv_2d = empty_state,
  Rv_3d = empty_state
)

#times <- seq(0,75, by = 1)
# Model fit ---------------------------------------------------------
# read in csv with breakpoints
breakpoints <- read_csv2("inst/extdata/inputs/breakpoints_for_model_fit_v3.csv") %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y"),
         time = as.numeric(date - date[1])) %>%
  select(date, time, variant, contact_matrix)

# run fit procedure
fits <- fit_to_data_func(breakpoints = breakpoints, params = params, init = init, 
                 case_data = osiris1, contact_matrices = cm_list,
                 vac_info = vac_rates_list, est_omega = FALSE,
                 save_output_to_file = FALSE, path_out = NULL)

# Run forward simulations --------------------------------------------
# seir_out <- lsoda(init, times, age_struct_seir_ode2, params)
# seir_out <- as.data.frame(seir_out)
# out <- postprocess_age_struct_model_output2(seir_out)
# cases <- params$sigma * rowSums(out$E + out$Ev_1d + out$Ev_2d + out$Ev_3d) * params$p_report
# plot(cases ~ times, type = "l")
