# Script for running scenarios for the Eupropean Scenario Hub
# URL: https://github.com/covid19-forecast-hub-europe/covid19-scenario-hub-europe#readme

# preamble ---------------------------------------------------------
# This script will load necessary packages, data sources, fit the model to data,
# and then run scenarios.
# All scenarios will be using Dutch data
# TODO: generalize to other European countries
# TODO: break up code into work chunks
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
library(here)

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
april_2017     <- readRDS(paste0(path,"transmission_matrix_april_2017.rds"))

# vaccination schedule ----------------------------------------------
# read in vaccination schedule
vac_schedule_A <- readRDS("inst/extdata/inputs/vaccination_schedules/vac_schedule_60plus_0.5_60minus_0.5.rds") 
vac_schedule_B <- readRDS("inst/extdata/inputs/vaccination_schedules/vac_schedule_60plus_0.5_60minus_0.2.rds") 
vac_schedule_C <- readRDS("inst/extdata/inputs/vaccination_schedules/vac_schedule_60plus_0.8_60minus_0.5.rds")
vac_schedule_D <- readRDS("inst/extdata/inputs/vaccination_schedules/vac_schedule_60plus_0.8_60minus_0.2.rds")
# read in xlsx file with VEs (there is 1 sheet for each variant)
ve_dat <- read_excel("inst/extdata/inputs/ve_estimates/ve_dat.xlsx", sheet = "omicron")

# specify initial model parameters ---------------------------------
# convert vaccination schedule to vaccination rates
vac_ratesA <- convert_vac_schedule2(
  vac_schedule = vac_schedule_A, ve_pars = ve_dat,
  wane = TRUE, k_inf = 0.006, k_sev = 0.012, t0 = 365)

vac_ratesB <- convert_vac_schedule2(
  vac_schedule = vac_schedule_B, ve_pars = ve_dat,
  wane = TRUE, k_inf = 0.006, k_sev = 0.012, t0 = 365)

vac_ratesC <- convert_vac_schedule2(
  vac_schedule = vac_schedule_C, ve_pars = ve_dat,
  wane = TRUE, k_inf = 0.006, k_sev = 0.012, t0 = 365)

vac_ratesD <- convert_vac_schedule2(
  vac_schedule = vac_schedule_D, ve_pars = ve_dat,
  wane = TRUE, k_inf = 0.006, k_sev = 0.012, t0 = 365)

# data wrangle for model input
df_inputA <- pivot_wider(vac_ratesA %>% 
                         filter(param != "comp_ve") %>%
                         mutate(param = ifelse(param == "comp_delay", "delay", param)), 
                        names_from = c("param", "age_group"), 
                        names_sep = "", values_from = "value")

df_inputB <- pivot_wider(vac_ratesB %>% 
                            filter(param != "comp_ve") %>%
                            mutate(param = ifelse(param == "comp_delay", "delay", param)), 
                          names_from = c("param", "age_group"), 
                          names_sep = "", values_from = "value")

df_inputC <- pivot_wider(vac_ratesC %>% 
                           filter(param != "comp_ve") %>%
                           mutate(param = ifelse(param == "comp_delay", "delay", param)), 
                         names_from = c("param", "age_group"), 
                         names_sep = "", values_from = "value")

df_inputD <- pivot_wider(vac_ratesD %>% 
                           filter(param != "comp_ve") %>%
                           mutate(param = ifelse(param == "comp_delay", "delay", param)), 
                         names_from = c("param", "age_group"), 
                         names_sep = "", values_from = "value")


# parameters must be in a named list
# scenarios A 
paramsA <- list(N = n_vec,
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
               alpha1 = df_inputA %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha2 = df_inputA %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha3 = df_inputA %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha4 = df_inputA %>%
                 filter(dose == "d4", outcome == "infection") %>%
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha5 = df_inputA %>%
                 filter(dose == "d5", outcome == "infection") %>%
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               # delay to protection
               delay1 = df_inputA %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay2 = df_inputA %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay3 = df_inputA %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay4 = df_inputA %>%
                 filter(dose == "d4", outcome == "infection") %>%
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay5 = df_inputA %>%
                 filter(dose == "d5", outcome == "infection") %>%
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               # protection against infection
               eta1 = df_inputA %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta2 = df_inputA %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta3 = df_inputA %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta4 = df_inputA %>%
                 filter(dose == "d4", outcome == "infection") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta5 = df_inputA %>%
                 filter(dose == "d5", outcome == "infection") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               # protection from hospitalisation
               eta_hosp1 = df_inputA %>% 
                 filter(dose == "d1", outcome == "hospitalisation") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp2 = df_inputA %>% 
                 filter(dose == "d2", outcome == "hospitalisation") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp3 = df_inputA %>% 
                 filter(dose == "d3", outcome == "hospitalisation") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp4 = df_inputA %>%
                 filter(dose == "d4", outcome == "hospitalisation") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp5 = df_inputA %>%
                 filter(dose == "d5", outcome == "hospitalisation") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               # protection from transmission
               eta_trans1 = df_inputA %>% 
                 filter(dose == "d1", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans2 = df_inputA %>% 
                 filter(dose == "d2", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans3 = df_inputA %>% 
                 filter(dose == "d3", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans4 = df_inputA %>%
                 filter(dose == "d4", outcome == "transmission") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans5 = df_inputA %>%
                 filter(dose == "d5", outcome == "transmission") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               p_report = p_reported_by_age,
               contact_mat = april_2017,
               calendar_start_date = as.Date("2020-01-01")
)

# scenarios B
paramsB <- list(N = n_vec,
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
                 alpha1 = df_inputB %>% 
                   filter(dose == "d1", outcome == "infection") %>% 
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha2 = df_inputB %>% 
                   filter(dose == "d2", outcome == "infection") %>% 
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha3 = df_inputB %>% 
                   filter(dose == "d3", outcome == "infection") %>% 
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha4 = df_inputB %>%
                   filter(dose == "d4", outcome == "infection") %>%
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha5 = df_inputB %>%
                   filter(dose == "d5", outcome == "infection") %>%
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 # delay to protection
                 delay1 = df_inputB %>% 
                   filter(dose == "d1", outcome == "infection") %>% 
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay2 = df_inputB %>% 
                   filter(dose == "d2", outcome == "infection") %>% 
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay3 = df_inputB %>% 
                   filter(dose == "d3", outcome == "infection") %>% 
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay4 = df_inputB %>%
                   filter(dose == "d4", outcome == "infection") %>%
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay5 = df_inputB %>%
                   filter(dose == "d5", outcome == "infection") %>%
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 # protection against infection
                 eta1 = df_inputB %>% 
                   filter(dose == "d1", outcome == "infection") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta2 = df_inputB %>% 
                   filter(dose == "d2", outcome == "infection") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta3 = df_inputB %>% 
                   filter(dose == "d3", outcome == "infection") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta4 = df_inputB %>%
                   filter(dose == "d4", outcome == "infection") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta5 = df_inputB %>%
                   filter(dose == "d5", outcome == "infection") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 # protection from hospitalisation
                 eta_hosp1 = df_inputB %>% 
                   filter(dose == "d1", outcome == "hospitalisation") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp2 = df_inputB %>% 
                   filter(dose == "d2", outcome == "hospitalisation") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp3 = df_inputB %>% 
                   filter(dose == "d3", outcome == "hospitalisation") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp4 = df_inputB %>%
                   filter(dose == "d4", outcome == "hospitalisation") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp5 = df_inputB %>%
                   filter(dose == "d5", outcome == "hospitalisation") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 # protection from transmission
                 eta_trans1 = df_inputB %>% 
                   filter(dose == "d1", outcome == "transmission") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans2 = df_inputB %>% 
                   filter(dose == "d2", outcome == "transmission") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans3 = df_inputB %>% 
                   filter(dose == "d3", outcome == "transmission") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans4 = df_inputB %>%
                   filter(dose == "d4", outcome == "transmission") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans5 = df_inputB %>%
                   filter(dose == "d5", outcome == "transmission") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 p_report = p_reported_by_age,
                 contact_mat = april_2017,
                 calendar_start_date = as.Date("2020-01-01")
)

# scenarios C
paramsC <- list(N = n_vec,
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
                alpha1 = df_inputC %>% 
                  filter(dose == "d1", outcome == "infection") %>% 
                  select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                alpha2 = df_inputC %>% 
                  filter(dose == "d2", outcome == "infection") %>% 
                  select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                alpha3 = df_inputC %>% 
                  filter(dose == "d3", outcome == "infection") %>% 
                  select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                alpha4 = df_inputC %>%
                  filter(dose == "d4", outcome == "infection") %>%
                  select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                alpha5 = df_inputC %>%
                  filter(dose == "d5", outcome == "infection") %>%
                  select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                # delay to protection
                delay1 = df_inputC %>% 
                  filter(dose == "d1", outcome == "infection") %>% 
                  select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                delay2 = df_inputC %>% 
                  filter(dose == "d2", outcome == "infection") %>% 
                  select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                delay3 = df_inputC %>% 
                  filter(dose == "d3", outcome == "infection") %>% 
                  select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                delay4 = df_inputC %>%
                  filter(dose == "d4", outcome == "infection") %>%
                  select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                delay5 = df_inputC %>%
                  filter(dose == "d5", outcome == "infection") %>%
                  select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                # protection against infection
                eta1 = df_inputC %>% 
                  filter(dose == "d1", outcome == "infection") %>% 
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta2 = df_inputC %>% 
                  filter(dose == "d2", outcome == "infection") %>% 
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta3 = df_inputC %>% 
                  filter(dose == "d3", outcome == "infection") %>% 
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta4 = df_inputC %>%
                  filter(dose == "d4", outcome == "infection") %>%
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta5 = df_inputC %>%
                  filter(dose == "d5", outcome == "infection") %>%
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                # protection from hospitalisation
                eta_hosp1 = df_inputC %>% 
                  filter(dose == "d1", outcome == "hospitalisation") %>% 
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta_hosp2 = df_inputC %>% 
                  filter(dose == "d2", outcome == "hospitalisation") %>% 
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta_hosp3 = df_inputC %>% 
                  filter(dose == "d3", outcome == "hospitalisation") %>% 
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta_hosp4 = df_inputC %>%
                  filter(dose == "d4", outcome == "hospitalisation") %>%
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta_hosp5 = df_inputC %>%
                  filter(dose == "d5", outcome == "hospitalisation") %>%
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                # protection from transmission
                eta_trans1 = df_inputC %>% 
                  filter(dose == "d1", outcome == "transmission") %>% 
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta_trans2 = df_inputC %>% 
                  filter(dose == "d2", outcome == "transmission") %>% 
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta_trans3 = df_inputC %>% 
                  filter(dose == "d3", outcome == "transmission") %>% 
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta_trans4 = df_inputC %>%
                  filter(dose == "d4", outcome == "transmission") %>%
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                eta_trans5 = df_inputC %>%
                  filter(dose == "d5", outcome == "transmission") %>%
                  select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                p_report = p_reported_by_age,
                contact_mat = april_2017,
                calendar_start_date = as.Date("2020-01-01")
)

# scenarios D
paramsD <- list(N = n_vec,
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
                 alpha1 = df_inputD %>% 
                   filter(dose == "d1", outcome == "infection") %>% 
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha2 = df_inputD %>% 
                   filter(dose == "d2", outcome == "infection") %>% 
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha3 = df_inputD %>% 
                   filter(dose == "d3", outcome == "infection") %>% 
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha4 = df_inputD %>%
                   filter(dose == "d4", outcome == "infection") %>%
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha5 = df_inputD %>%
                   filter(dose == "d5", outcome == "infection") %>%
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 # delay to protection
                 delay1 = df_inputD %>% 
                   filter(dose == "d1", outcome == "infection") %>% 
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay2 = df_inputD %>% 
                   filter(dose == "d2", outcome == "infection") %>% 
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay3 = df_inputD %>% 
                   filter(dose == "d3", outcome == "infection") %>% 
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay4 = df_inputD %>%
                   filter(dose == "d4", outcome == "infection") %>%
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay5 = df_inputD %>%
                   filter(dose == "d5", outcome == "infection") %>%
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 # protection against infection
                 eta1 = df_inputD %>% 
                   filter(dose == "d1", outcome == "infection") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta2 = df_inputD %>% 
                   filter(dose == "d2", outcome == "infection") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta3 = df_inputD %>% 
                   filter(dose == "d3", outcome == "infection") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta4 = df_inputD %>%
                   filter(dose == "d4", outcome == "infection") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta5 = df_inputD %>%
                   filter(dose == "d5", outcome == "infection") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 # protection from hospitalisation
                 eta_hosp1 = df_inputD %>% 
                   filter(dose == "d1", outcome == "hospitalisation") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp2 = df_inputD %>% 
                   filter(dose == "d2", outcome == "hospitalisation") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp3 = df_inputD %>% 
                   filter(dose == "d3", outcome == "hospitalisation") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp4 = df_inputD %>%
                   filter(dose == "d4", outcome == "hospitalisation") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp5 = df_inputD %>%
                   filter(dose == "d5", outcome == "hospitalisation") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 # protection from transmission
                 eta_trans1 = df_inputD %>% 
                   filter(dose == "d1", outcome == "transmission") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans2 = df_inputD %>% 
                   filter(dose == "d2", outcome == "transmission") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans3 = df_inputD %>% 
                   filter(dose == "d3", outcome == "transmission") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans4 = df_inputD %>%
                   filter(dose == "d4", outcome == "transmission") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans5 = df_inputD %>%
                   filter(dose == "d5", outcome == "transmission") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 p_report = p_reported_by_age,
                 contact_mat = april_2017,
                 calendar_start_date = as.Date("2020-01-01")
)
# Specify initial conditions ---------------------------------------------------
init_cond_list <- readRDS("inst/extdata/results/model_fits/initial_conditions_2022-08-19.rds")
init_cond <- unlist(init_cond_list[[length(init_cond_list)]])
#init_cond[1] <- 872

# ------------------------------------------------------------------------------
# Run forward simulations ------------------------------------------------------
t_start <- init_cond[1]
t_end <- 1306#t_start + 365
times <- as.integer(seq(t_start, t_end, by = 1))
betas <- readRDS("inst/extdata/results/model_fits/beta_draws_2022-08-19.rds")
# sample 100 betas from last time window
betas100 <- sample(betas[[length(betas)]]$beta, 100)

# register parallel backend
n_cores = 15
registerDoParallel(cores=n_cores)
n_sim <- 100
# Scenario A
# Fall booster campaign in 60+, optimistic VE
scenarioA <- foreach(i = 1:n_sim) %dopar% {
  paramsA$beta <- betas100[i]
  paramsA$contact_mat <- april_2017[[i]]
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_cond, times, age_struct_seir_ode2, paramsA, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioA, "/rivm/s/ainsliek/results/scenarioA.rds")
doParallel::stopImplicitCluster()
# Scenario B
# Fall booster campaign in 18+, optimistic VE
# register parallel backend
registerDoParallel(cores=n_cores)
scenarioB <- foreach(i = 1:n_sim) %dopar% {
  paramsB$beta <- betas100[i]
  paramsB$contact_mat <- april_2017[[i]]
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_cond, times, age_struct_seir_ode2, paramsB, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioB, "/rivm/s/ainsliek/results/scenarioB.rds")
doParallel::stopImplicitCluster()
# Scenario C
# Fall booster campaign in 60+, pessimistic VE
# register parallel backend
registerDoParallel(cores=n_cores)
scenarioC <- foreach(i = 1:n_sim) %dopar% {
  paramsC$beta <- betas100[i]
  paramsC$contact_mat <- april_2017[[i]]
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_cond, times, age_struct_seir_ode2, paramsC, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioC, "/rivm/s/ainsliek/results/scenarioC.rds")
doParallel::stopImplicitCluster()
# Scenario D
# Fall booster campaign in 18+, pessimistic VE
# register parallel backend
registerDoParallel(cores=n_cores)
scenarioD <- foreach(i = 1:n_sim) %dopar% {
  paramsD$beta <- betas100[i]
  paramsD$contact_mat <- april_2017[[i]]
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_cond, times, age_struct_seir_ode2, paramsD, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioD, "/rivm/s/ainsliek/results/scenarioD.rds")
doParallel::stopImplicitCluster()
#-------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Post-process scenario runs ---------------------------------------------------
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

# wrangle Scenario A output ----------------------------------------------------
p_report_vec <- c(rep(as.numeric(paramsA$p_report),6))

# read in saved output from model runs
scenarioA <- readRDS("/rivm/s/ainsliek/results/scenarioA.rds")
#scenarioA <- readRDS("C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/ECDC Scenario Modelling Hub/round 1/scenarioA.rds")
sim <- length(scenarioA)
# loop over samples and summarise results
outA <- list()
for(s in 1:sim){
  seir_output <- postprocess_age_struct_model_output2(scenarioA[[s]])
  paramsA$beta <- betas100[s]
  paramsA$contact_mat <- april_2017[[s]]
  seir_outcomes <- summarise_results(seir_output, params = paramsA, t_vec = times) %>%
    mutate(sample = s)
  outA[[s]] <- seir_outcomes
}
dfA <- bind_rows(outA) %>%
  mutate(scenario_id = "A-2022-07-24") 

# wrangle Scenario B output ----------------------------------------------------
# read in saved output from model runs
scenarioB <- readRDS("/rivm/s/ainsliek/results/scenarioB.rds")
#scenarioB <- readRDS("C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/ECDC Scenario Modelling Hub/round 1/scenarioB.rds")
sim <- length(scenarioB)
# loop over samples and summarise results
outB <- list()
for(s in 1:sim){
  seir_output <- postprocess_age_struct_model_output2(scenarioB[[s]])
  paramsB$beta <- betas100[s]
  paramsB$contact_mat <- april_2017[[s]]
  seir_outcomes <- summarise_results(seir_output, params = paramsB, t_vec = times) %>%
    mutate(sample = s)
  outB[[s]] <- seir_outcomes
}
dfB <- bind_rows(outB) %>%
  mutate(scenario_id = "B-2022-07-24") 

# wrangle Scenario C output ----------------------------------------------------
# read in saved output from model runs
scenarioC <- readRDS("/rivm/s/ainsliek/results/scenarioC.rds")
#scenarioC <- readRDS("C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/ECDC Scenario Modelling Hub/round 1/scenarioC.rds")
sim <- length(scenarioC)
# loop over samples and summarise results
outC <- list()
for(s in 1:sim){
  seir_output <- postprocess_age_struct_model_output2(scenarioC[[s]])
  paramsC$beta <- betas100[s]
  paramsC$contact_mat <- april_2017[[s]]
  seir_outcomes <- summarise_results(seir_output, params = paramsC, t_vec = times) %>%
    mutate(sample = s)
  outC[[s]] <- seir_outcomes
}
dfC <- bind_rows(outC) %>%
  mutate(scenario_id = "C-2022-07-24") 

# wrangle Scenario D output ----------------------------------------------------
# read in saved output from model runs
scenarioD <- readRDS("/rivm/s/ainsliek/results/scenarioD.rds")
#scenarioD <- readRDS("C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/ECDC Scenario Modelling Hub/round 1/scenarioD.rds")
sim <- length(scenarioD)
# loop over samples and summarise results
outD <- list()
for(s in 1:sim){
  seir_output <- postprocess_age_struct_model_output2(scenarioD[[s]])
  paramsD$beta <- betas100[s]
  paramsD$contact_mat <- april_2017[[s]]
  seir_outcomes <- summarise_results(seir_output, params = paramsD, t_vec = times) %>%
    mutate(sample = s)
  outD[[s]] <- seir_outcomes
}
dfD <- bind_rows(outD) %>%
  mutate(scenario_id = "D-2022-07-24") 

# ------------------------------------------------------------------------------
# join all scenarios in a single data frame
df_round2 <- bind_rows(dfA, dfB, dfC, dfD) 

# output for plotting
saveRDS(df_round2, "inst/extdata/results/2022-08-19-rivm-vacamole.rds")

# put all scenarios together into single data frame and sum over epiweek & 
# age groups
df_round2_sh <- df_round2 %>%
  group_by(scenario_id, sample, epiweek, horizon, target_variable) %>%
  summarise_at(.vars = "value", .funs = "sum") %>%
  ungroup() %>%
  mutate(value = round(value),
         origin_date = as.Date("2022-07-24"),
         horizon_ = as.numeric(substr(horizon, start = 1, stop = 2)),
         target_end_date = origin_date + weeks(horizon_),
         horizon_ = NULL,
         location = "NL") %>%
  select(-epiweek) %>%
  arrange(scenario_id, sample, horizon)

# output for submission to scenario hub
#write_csv(df_round2_sh, "C:/Users/ainsliek/Documents/covid19-scenario-hub-europe/data-processed/RIVM-vacamole/2022-07-24-RIVM-vacamole.csv")
write_csv(df_round2_sh, "/rivm/s/ainsliek/results/scenario_hub/round2/2022-07-24-RIVM-vacamole.csv")
