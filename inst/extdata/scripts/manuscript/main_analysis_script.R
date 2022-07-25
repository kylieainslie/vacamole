# ------------------------------------------------------------------
# Forward simulations to determine impacts of vaccinating
# 12-17 year olds for manuscript
# ------------------------------------------------------------------

# Options ----------------------------------------------------------
# suppress dplyr::summarise() warnings
options(dplyr.summarise.inform = FALSE)

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
library(foreach)
library(doParallel)

source("R/convert_vac_schedule2.R")
source("R/na_to_zero.R")
source("R/calc_waning.R")
source("R/age_struct_seir_ode2.R")
source("R/postprocess_age_struct_model_output2.R")
source("R/summarise_results.R")
source("R/get_foi.R")
# -------------------------------------------------------------------
# Define population size --------------------------------------------
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 
              0.13083463,0.14514332, 0.12092904, 0.08807406, 
              0.04622194)
n <- 17407585 # Dutch population size
n_vec <- n * age_dist

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

# 60% reduction after 8 months
wane_8months <- uniroot(Fk, c(0,1), tau = 244, p = 0.6)$root

# contact matrices --------------------------------------------------
path <- "/rivm/s/ainsliek/data/contact_matrices/converted/"
# path <- "inst/extdata/inputs/contact_matrices/converted/"
april_2017 <- readRDS(paste0(path,"transmission_matrix_april_2017.rds"))
june_2021  <- readRDS(paste0(path,"transmission_matrix_june_2021.rds"))
# VE estimates -----------------------------------------------------------------
# read in xlsx file with VEs (there is 1 sheet for each variant)
# we'll only use wildtype values for now
# wt_ve <- read_excel("inst/extdata/inputs/ve_estimates/ve_dat_round2_AB.xlsx", sheet = "wildtype") 
# alpha_ve <- read_excel("inst/extdata/inputs/ve_estimates/ve_dat_round2_AB.xlsx", sheet = "alpha")
delta_ve <- read_excel("inst/extdata/inputs/ve_estimates/ve_dat_round2_AB.xlsx", sheet = "delta")
# omicron_ve <- read_excel("inst/extdata/inputs/ve_estimates/ve_dat_round2_AB.xlsx", sheet = "omicron")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------
# read in vac schedules --------------------------------------------
# ------------------------------------------------------------------
basis_5plus <- read_csv("inst/extdata/inputs/vaccination_schedules/vac_schedule_5plus.csv") %>%
  select(-starts_with("X")) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

basis_12plus <- read_csv("inst/extdata/inputs/vaccination_schedules/vac_schedule_12plus.csv") %>%
  select(-starts_with("X")) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

basis_18plus <- read_csv("inst/extdata/inputs/vaccination_schedules/vac_schedule_18plus.csv") %>%
  select(-starts_with("X")) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

# add extra dates for forward sims ---------------------------------------------
extra_start_date <- tail(basis_12plus$date,1) + 1
extra_end_date <- as.Date("2022-03-31")
extra_dates <- seq.Date(from = as.Date(extra_start_date), 
                        to = as.Date(extra_end_date), by = 1)

vac_schedule_5plus <- data.frame(date = extra_dates) %>%
  full_join(basis_5plus, ., by = "date") %>%
  fill(-.data$date)

vac_schedule_12plus <- data.frame(date = extra_dates) %>%
  full_join(basis_12plus, ., by = "date") %>%
  fill(-.data$date)

vac_schedule_18plus <- data.frame(date = extra_dates) %>%
  full_join(basis_18plus, ., by = "date") %>%
  fill(-.data$date)

# add columns for 3rd, 4th, and 5th doses of pfizer and moderna vaccines -------
new_columns <- c(paste0("pf_d3_", 1:10),
                 paste0("pf_d4_", 1:10),
                 paste0("pf_d5_", 1:10),
                 paste0("mo_d3_", 1:10),
                 paste0("mo_d4_", 1:10),
                 paste0("mo_d5_", 1:10))

vac_schedule_5plus[,new_columns[c(1:9, 11:19, 21:29,
                                  31:39, 41:49, 51:59)]] <- 0
vac_schedule_12plus[,new_columns] <- 0
vac_schedule_18plus[,new_columns] <- 0

# rearrange columns to preserve correct order (also exclude Novovax doses)
vac_schedule_5plus1 <- vac_schedule_5plus %>%
  select(date, pf_d1_1:pf_d2_9, pf_d3_1:pf_d5_9,
         mo_d1_1:mo_d2_9, mo_d3_1:mo_d5_9,
         az_d1_1:az_d2_9, ja_d1_1:ja_d2_9)

vac_schedule_12plus1 <- vac_schedule_12plus %>%
  select(date, pf_d1_1:pf_d2_10, pf_d3_1:pf_d5_10,
         mo_d1_1:mo_d2_10, mo_d3_1:mo_d5_10,
         az_d1_1:az_d2_10, ja_d1_1:ja_d2_10)

vac_schedule_18plus1 <- vac_schedule_18plus %>%
  select(date, pf_d1_1:pf_d2_10, pf_d3_1:pf_d5_10,
         mo_d1_1:mo_d2_10, mo_d3_1:mo_d5_10,
         az_d1_1:az_d2_10, ja_d1_1:ja_d2_10)

# ------------------------------------------------------------------
# convert vac schedules --------------------------------------------
# ------------------------------------------------------------------
vac_rates_5plus <- convert_vac_schedule2(
  vac_schedule = vac_schedule_5plus1, ve_pars = delta_ve,
  wane = TRUE, k_inf = 0.006, k_sev = 0.012, t0 = 365)


vac_rates_12plus <- convert_vac_schedule2(
  vac_schedule = vac_schedule_12plus1, ve_pars = delta_ve,
  wane = TRUE, k_inf = 0.006, k_sev = 0.012, t0 = 365)


vac_rates_18plus <- convert_vac_schedule2(
  vac_schedule = vac_schedule_18plus1, ve_pars = delta_ve,
  wane = TRUE, k_inf = 0.006, k_sev = 0.012, t0 = 365)

# data wrangle for model input -------------------------------------------------
df_input_5plus <- pivot_wider(vac_rates_5plus %>%
  filter(param != "comp_ve") %>%
  mutate(param = ifelse(param == "comp_delay", "delay",
    param
  )),
names_from = c("param", "age_group"),
names_sep = "", values_from = "value")

df_input_12plus <- pivot_wider(vac_rates_12plus %>%
  filter(param != "comp_ve") %>%
  mutate(param = ifelse(param == "comp_delay", "delay",
    param
  )),
names_from = c("param", "age_group"),
names_sep = "", values_from = "value")

df_input_18plus <- pivot_wider(vac_rates_18plus %>%
  filter(param != "comp_ve") %>%
  mutate(param = ifelse(param == "comp_delay", "delay",
    param
  )),
names_from = c("param", "age_group"),
names_sep = "", values_from = "value")
# ------------------------------------------------------------------------------
# parameters lists--------------------------------------------------------------
# ------------------------------------------------------------------------------

# vac in 5+ --------------------------------------------------------------------
params_5plus <- list(
  N = n_vec,
  beta = 0.0004,
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
  omega = wane_8months,
  # --- daily vaccination rate ---
  alpha1 = df_input_5plus %>% 
    filter(dose == "d1", outcome == "infection") %>% 
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha2 = df_input_5plus %>% 
    filter(dose == "d2", outcome == "infection") %>% 
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha3 = df_input_5plus %>% 
    filter(dose == "d3", outcome == "infection") %>% 
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha4 = df_input_5plus %>%
    filter(dose == "d4", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha5 = df_input_5plus %>%
    filter(dose == "d5", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  # --- delay to protection ---
  delay1 = df_input_5plus %>% 
    filter(dose == "d1", outcome == "infection") %>% 
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay2 = df_input_5plus %>% 
    filter(dose == "d2", outcome == "infection") %>% 
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay3 = df_input_5plus %>% 
    filter(dose == "d3", outcome == "infection") %>% 
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay4 = df_input_5plus %>%
    filter(dose == "d4", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay5 = df_input_5plus %>%
    filter(dose == "d5", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  # --- protection against infection ---
  eta1 = df_input_5plus %>% 
    filter(dose == "d1", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta2 = df_input_5plus %>% 
    filter(dose == "d2", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta3 = df_input_5plus %>% 
    filter(dose == "d3", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta4 = df_input_5plus %>%
    filter(dose == "d4", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta5 = df_input_5plus %>%
    filter(dose == "d5", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  # protection from hospitalisation
  eta_hosp1 = df_input_5plus %>% 
    filter(dose == "d1", outcome == "hospitalisation") %>% 
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp2 = df_input_5plus %>%
    filter(dose == "d2", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp3 = df_input_5plus %>% 
    filter(dose == "d3", outcome == "hospitalisation") %>% 
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp4 = df_input_5plus %>%
    filter(dose == "d4", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp5 = df_input_5plus %>%
    filter(dose == "d5", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  # protection from transmission
  eta_trans1 = df_input_5plus %>% 
    filter(dose == "d1", outcome == "transmission") %>% 
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans2 = df_input_5plus %>% 
    filter(dose == "d2", outcome == "transmission") %>% 
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans3 = df_input_5plus %>% 
    filter(dose == "d3", outcome == "transmission") %>% 
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans4 = df_input_5plus %>%
    filter(dose == "d4", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans5 = df_input_5plus %>%
    filter(dose == "d5", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  p_report = p_reported_by_age,
  contact_mat = june_2021,
  calendar_start_date = as.Date("2020-01-01")
)

# vac in 12+ -------------------------------------------------------------------
params_12plus <- list(
  N = n_vec,
  beta = 0.0004,
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
  omega = wane_8months,
  # --- daily vaccination rate ---
  alpha1 = df_input_12plus %>%
    filter(dose == "d1", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha2 = df_input_12plus %>%
    filter(dose == "d2", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha3 = df_input_12plus %>%
    filter(dose == "d3", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha4 = df_input_12plus %>%
    filter(dose == "d4", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha5 = df_input_12plus %>%
    filter(dose == "d5", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  # --- delay to protection ---
  delay1 = df_input_12plus %>%
    filter(dose == "d1", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay2 = df_input_12plus %>%
    filter(dose == "d2", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay3 = df_input_12plus %>%
    filter(dose == "d3", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay4 = df_input_12plus %>%
    filter(dose == "d4", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay5 = df_input_12plus %>%
    filter(dose == "d5", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  # --- protection against infection ---
  eta1 = df_input_12plus %>%
    filter(dose == "d1", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta2 = df_input_12plus %>%
    filter(dose == "d2", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta3 = df_input_12plus %>%
    filter(dose == "d3", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta4 = df_input_12plus %>%
    filter(dose == "d4", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta5 = df_input_12plus %>%
    filter(dose == "d5", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  # protection from hospitalisation
  eta_hosp1 = df_input_12plus %>%
    filter(dose == "d1", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp2 = df_input_12plus %>%
    filter(dose == "d2", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp3 = df_input_12plus %>%
    filter(dose == "d3", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp4 = df_input_12plus %>%
    filter(dose == "d4", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp5 = df_input_12plus %>%
    filter(dose == "d5", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  # protection from transmission
  eta_trans1 = df_input_12plus %>%
    filter(dose == "d1", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans2 = df_input_12plus %>%
    filter(dose == "d2", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans3 = df_input_12plus %>%
    filter(dose == "d3", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans4 = df_input_12plus %>%
    filter(dose == "d4", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans5 = df_input_12plus %>%
    filter(dose == "d5", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  p_report = p_reported_by_age,
  contact_mat = june_2021,
  calendar_start_date = as.Date("2020-01-01")
)

# vac in 18+ -------------------------------------------------------------------
params_18plus <- list(
  N = n_vec,
  beta = 0.0004,
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
  omega = wane_8months,
  # daily vaccination rate
  alpha1 = df_input_18plus %>%
    filter(dose == "d1", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha2 = df_input_18plus %>%
    filter(dose == "d2", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha3 = df_input_18plus %>%
    filter(dose == "d3", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha4 = df_input_18plus %>%
    filter(dose == "d4", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  alpha5 = df_input_18plus %>%
    filter(dose == "d5", outcome == "infection") %>%
    select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
  # delay to protection
  delay1 = df_input_18plus %>%
    filter(dose == "d1", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay2 = df_input_18plus %>%
    filter(dose == "d2", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay3 = df_input_18plus %>%
    filter(dose == "d3", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay4 = df_input_18plus %>%
    filter(dose == "d4", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  delay5 = df_input_18plus %>%
    filter(dose == "d5", outcome == "infection") %>%
    select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
  # protection against infection
  eta1 = df_input_18plus %>%
    filter(dose == "d1", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta2 = df_input_18plus %>%
    filter(dose == "d2", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta3 = df_input_18plus %>%
    filter(dose == "d3", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta4 = df_input_18plus %>%
    filter(dose == "d4", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta5 = df_input_18plus %>%
    filter(dose == "d5", outcome == "infection") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  # protection from hospitalisation
  eta_hosp1 = df_input_18plus %>%
    filter(dose == "d1", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp2 = df_input_18plus %>%
    filter(dose == "d2", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp3 = df_input_18plus %>%
    filter(dose == "d3", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp4 = df_input_18plus %>%
    filter(dose == "d4", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_hosp5 = df_input_18plus %>%
    filter(dose == "d5", outcome == "hospitalisation") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  # protection from transmission
  eta_trans1 = df_input_18plus %>%
    filter(dose == "d1", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans2 = df_input_18plus %>%
    filter(dose == "d2", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans3 = df_input_18plus %>%
    filter(dose == "d3", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans4 = df_input_18plus %>%
    filter(dose == "d4", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  eta_trans5 = df_input_18plus %>%
    filter(dose == "d5", outcome == "transmission") %>%
    select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
  p_report = p_reported_by_age,
  contact_mat = june_2021,
  calendar_start_date = as.Date("2020-01-01")
)

# Specify initial conditions ---------------------------------------------------
init_cond_list <- readRDS("inst/extdata/results/model_fits/manuscript/initial_conditions_manuscript.rds")
init_cond <- unlist(init_cond_list[[length(init_cond_list)]])

# ------------------------------------------------------------------------------
# Run forward simulations ------------------------------------------------------
# ------------------------------------------------------------------------------
t_start <- init_cond[1]
t_end1 <- t_start + 132
t_end2 <- 820  # 31 March 2022
times1 <- as.integer(seq(t_start, t_end1, by = 1))
times2 <- as.integer(seq(t_end1, t_end2, by = 1))
betas <- readRDS("inst/extdata/results/model_fits/manuscript/beta_draws.rds")

# sample betas from last time window
n_sim <- 200
betas_sample <- sample(betas[[length(betas)]]$beta, n_sim)
betas_sample2<- rnorm(n_sim, mean = 0.00087, sd = 7.1145e-6)

# register parallel backend
n_cores <- ifelse(n_sim < 15, n_sim, 15)
registerDoParallel(cores = n_cores)

#-------------------------------------------------------------------------------
# Vaccination in 5+ ------------------------------------------------------------
scenario_5plus <- foreach(i = 1:n_sim) %dopar% {
  rk45 <- rkMethod("rk45dp7")
  
  # run model from 22 June 2021 to 1 November 2021
  params_5plus$beta <- betas_sample[i]
  params_5plus$contact_mat <- june_2021[[i]]
  seir_out1 <- ode(init_cond, times1, age_struct_seir_ode2, params_5plus, method = rk45)
  
  # run model from 1 November 2021 to 31 March 2022
  init_cond2 <-  unlist(tail(as.data.frame(seir_out1), 1)[-1])
  params_5plus$beta <- betas_sample2[i]
  params_5plus$contact_mat <- april_2017[[i]]
  seir_out2 <- ode(init_cond2, times2, age_struct_seir_ode2, params_5plus, method = rk45)
  
  bind_rows(as.data.frame(seir_out1),as.data.frame(seir_out2))
}
saveRDS(scenario_5plus, "/rivm/s/ainsliek/results/impact_vac/resubmission/results_5plus.rds")

# Vaccination in 12+ -----------------------------------------------------------
scenario_12plus <- foreach(i = 1:n_sim) %dopar% {
  rk45 <- rkMethod("rk45dp7")
  
  # run model from 22 June 2021 to 1 November 2021
  params_12plus$beta <- betas_sample[i]
  params_12plus$contact_mat <- june_2021[[i]]
  seir_out1 <- ode(init_cond, times1, age_struct_seir_ode2, params_12plus, method = rk45)
  
  # run model from 1 November 2021 to 31 March 2022
  init_cond2 <-  unlist(tail(as.data.frame(seir_out1), 1)[-1])
  params_12plus$beta <- betas_sample2[i]
  params_12plus$contact_mat <- april_2017[[i]]
  seir_out2 <- ode(init_cond2, times2, age_struct_seir_ode2, params_12plus, method = rk45)
  
  bind_rows(as.data.frame(seir_out1),as.data.frame(seir_out2))
}
saveRDS(scenario_12plus, "/rivm/s/ainsliek/results/impact_vac/resubmission/results_12plus.rds")

# Vaccination in 18+ -----------------------------------------------------------
scenario_18plus <- foreach(i = 1:n_sim) %dopar% {
  rk45 <- rkMethod("rk45dp7")
  
  # run model from 22 June 2021 to 1 November 2021
  params_18plus$beta <- betas_sample[i]
  params_18plus$contact_mat <- june_2021[[i]]
  seir_out1 <- ode(init_cond, times1, age_struct_seir_ode2, params_18plus, method = rk45)
  
  # run model from 1 November 2021 to 31 March 2022
  init_cond2 <-  unlist(tail(as.data.frame(seir_out1), 1)[-1])
  params_18plus$beta <- betas_sample2[i]
  params_18plus$contact_mat <- april_2017[[i]]
  seir_out2 <- ode(init_cond2, times2, age_struct_seir_ode2, params_18plus, method = rk45)
  
  bind_rows(as.data.frame(seir_out1),as.data.frame(seir_out2))
}
saveRDS(scenario_18plus, "/rivm/s/ainsliek/results/impact_vac/resubmission/results_18plus.rds")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Post-process scenario runs ---------------------------------------------------
#-------------------------------------------------------------------------------
# Results must be in a csv file that contains only the following columns (in any
# order). No additional columns are allowed.
# - origin_date (date):	Date as YYYY-MM-DD, last day (Monday) of submission window
# - scenario_id	(string):	A specified "scenario ID"
# - target_variable	(string):	"inc case", "inc death", "inc hosp", "inc icu", "inc infection"
# - horizon	(string):	The time horizon for the projection relative to the origin_date (e.g., X wk)
# - target_end_date	(date):	Date as YYYY-MM-DD, the last day (Saturday) of the target week
# - location (string): An ISO-2 country code or "H0" ("NL" for the Netherlands)
# - sample	(numeric):	A integer corresponding to the sample index (used to plot trajectories)
# - value	(numeric):	The projected count, a non-negative integer number of new cases or deaths in the epidemiological week

# wrangle output ---------------------------------------------------------------
p_report_vec <- c(rep(as.numeric(params_5plus$p_report),6))

# Vaccination in 5+ ------------------------------------------------------------
# read in saved output from model runs
scenario_5plus <- readRDS("/rivm/s/ainsliek/results/impact_vac/resubmission/results_5plus.rds")
sim <- length(scenario_5plus)
# loop over samples and summarise results
out_5plus <- list()
for(s in 1:sim){
  out1 <- scenario_5plus[[s]][which(scenario_5plus[[s]]$time %in% times1),] %>%
    distinct()
  out2 <- scenario_5plus[[s]][which(scenario_5plus[[s]]$time %in% times2),] %>%
    distinct
  
  seir_output1 <- postprocess_age_struct_model_output2(out1)
  params_5plus$beta <- betas_sample[s]
  params_5plus$contact_mat <- june_2021[[s]]
  seir_outcomes1 <- summarise_results(seir_output = seir_output1, params = params_5plus, t_vec = times1) %>%
    mutate(sample = s)
  
  seir_output2 <- postprocess_age_struct_model_output2(out2)
  params_5plus$beta <- betas_sample2[s]
  params_5plus$contact_mat <- april_2017[[s]]
  seir_outcomes2 <- summarise_results(seir_output2, params = params_5plus, t_vec = times2) %>%
    mutate(sample = s)
  
  out_5plus[[s]] <- bind_rows(seir_outcomes1, seir_outcomes2)
}
df_5plus <- bind_rows(out_5plus) %>%
  mutate(scenario_id = "Vaccination in 5+") #%>%
  #filter(horizon != "53 wk")


# Vaccination in 12+ -----------------------------------------------------------
# read in saved output from model runs
scenario_12plus <- readRDS("/rivm/s/ainsliek/results/impact_vac/resubmission/results_12plus.rds")
sim <- length(scenario_12plus)
# loop over samples and summarise results
out_12plus <- list()
for(s in 1:sim){
  out1 <- scenario_12plus[[s]][which(scenario_12plus[[s]]$time %in% times1),] %>%
    distinct()
  out2 <- scenario_12plus[[s]][which(scenario_12plus[[s]]$time %in% times2),] %>%
    distinct
  
  seir_output1 <- postprocess_age_struct_model_output2(out1)
  params_12plus$beta <- betas_sample[s]
  params_12plus$contact_mat <- june_2021[[s]]
  seir_outcomes1 <- summarise_results(seir_output = seir_output1, params = params_12plus, t_vec = times1) %>%
    mutate(sample = s)
  
  seir_output2 <- postprocess_age_struct_model_output2(out2)
  params_12plus$beta <- betas_sample2[s]
  params_12plus$contact_mat <- april_2017[[s]]
  seir_outcomes2 <- summarise_results(seir_output2, params = params_12plus, t_vec = times2) %>%
    mutate(sample = s)
  
  out_12plus[[s]] <- bind_rows(seir_outcomes1, seir_outcomes2)
}
df_12plus <- bind_rows(out_12plus) %>%
  mutate(scenario_id = "Vaccination in 12+")# %>%
  #filter(horizon != "53 wk")

# Vaccination ni 18+ -----------------------------------------------------------
# read in saved output from model runs
scenario_18plus <- readRDS("/rivm/s/ainsliek/results/impact_vac/resubmission/results_18plus.rds")
sim <- length(scenario_18plus)
# loop over samples and summarise results
out_18plus <- list()
for(s in 1:sim){
  out1 <- scenario_18plus[[s]][which(scenario_18plus[[s]]$time %in% times1),] %>%
    distinct()
  out2 <- scenario_18plus[[s]][which(scenario_18plus[[s]]$time %in% times2),] %>%
    distinct
  
  seir_output1 <- postprocess_age_struct_model_output2(out1)
  params_18plus$beta <- betas_sample[s]
  params_18plus$contact_mat <- june_2021[[s]]
  seir_outcomes1 <- summarise_results(seir_output = seir_output1, params = params_18plus, t_vec = times1) %>%
    mutate(sample = s)
  
  seir_output2 <- postprocess_age_struct_model_output2(out2)
  params_18plus$beta <- betas_sample2[s]
  params_18plus$contact_mat <- april_2017[[s]]
  seir_outcomes2 <- summarise_results(seir_output2, params = params_18plus, t_vec = times2) %>%
    mutate(sample = s)
  
  out_18plus[[s]] <- bind_rows(seir_outcomes1, seir_outcomes2)
}
df_18plus <- bind_rows(out_18plus) %>%
  mutate(scenario_id = "Vaccination in 18+") #%>%
  #filter(horizon != "53 wk")

# ------------------------------------------------------------------------------
# join all scenarios in a single data frame
df_all <- bind_rows(df_5plus, df_12plus, df_18plus) 

# output for plotting
saveRDS(df_all, "/rivm/s/ainsliek/results/impact_vac/resubmission/results_all.rds")
