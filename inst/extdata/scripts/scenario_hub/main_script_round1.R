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

wane_3months <- uniroot(Fk, c(0,1), tau = 92, p = 0.6)$root
wane_8months <- uniroot(Fk, c(0,1), tau = 244, p = 0.6)$root
# 50% reduction after 6 months (used for model fits)
wane_6months <- uniroot(Fk, c(0,1), tau = 182, p = 0.5)$root
# contact matrices --------------------------------------------------
# path <- "/rivm/s/ainsliek/data/contact_matrices/converted/"
path <- "inst/extdata/inputs/contact_matrices/converted/"
april_2017     <- readRDS(paste0(path,"transmission_matrix_april_2017.rds"))
# april_2020     <- readRDS(paste0(path,"transmission_matrix_april_2020.rds"))
# june_2020      <- readRDS(paste0(path,"transmission_matrix_june_2020.rds"))
# september_2020 <- readRDS(paste0(path,"transmission_matrix_september_2020.rds"))
# february_2021  <- readRDS(paste0(path,"transmission_matrix_february_2021.rds"))
# june_2021      <- readRDS(paste0(path,"transmission_matrix_june_2021.rds"))
# november_2021  <- readRDS(paste0(path,"transmission_matrix_november_2021.rds"))

# put contact matrices into a list for input into fit_to_data_func()
# cm_list <- list(
#   april_2017 = april_2017,
#   april_2020 = april_2020,
#   june_2020 = june_2020,
#   september_2020 = september_2020,
#   february_2021 = february_2021,
#   june_2021 = june_2021,
#   november_2021 = november_2021
# )

# vaccination schedule ----------------------------------------------
# read in vaccination schedule
vac_scheduleAC <- readRDS("inst/extdata/inputs/vac_schedule_scenario_hub_round1_AC.rds") 
vac_scheduleBD <- readRDS("inst/extdata/inputs/vac_schedule_scenario_hub_round1_BD.rds") 

# read in xlsx file with VEs (there is 1 sheet for each variant)
omicron_ve <- read_excel("inst/extdata/inputs/ve_dat.xlsx", sheet = "omicron")

# specify initial model parameters ---------------------------------
# parameters must be in a named list
vac_ratesAC <- convert_vac_schedule2(
  vac_schedule = vac_scheduleAC,
  ve_pars = omicron_ve,
  wane = TRUE)

vac_ratesBD <- convert_vac_schedule2(
  vac_schedule = vac_scheduleBD,
  ve_pars = omicron_ve,
  wane = TRUE)

# data wrangle for model input
df_inputAC <- pivot_wider(vac_ratesAC %>% 
                          filter(param != "comp_ve") %>%
                          mutate(param = ifelse(param == "comp_delay", "delay", param)), 
                        names_from = c("param", "age_group"), 
                        names_sep = "", values_from = "value")

df_inputBD <- pivot_wider(vac_ratesBD %>% 
                            filter(param != "comp_ve") %>%
                            mutate(param = ifelse(param == "comp_delay", "delay", param)), 
                          names_from = c("param", "age_group"), 
                          names_sep = "", values_from = "value")


# parameters must be in a named list
# scenarios A & C
paramsAC <- list(N = n_vec,
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
               omega = 0.02017514,
               # daily vaccination rate
               alpha1 = df_inputAC %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha2 = df_inputAC %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha3 = df_inputAC %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha4 = df_inputAC %>%
                 filter(dose == "d4", outcome == "infection") %>%
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               alpha5 = df_inputAC %>%
                 filter(dose == "d5", outcome == "infection") %>%
                 select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
               # delay to protection
               delay1 = df_inputAC %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay2 = df_inputAC %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay3 = df_inputAC %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay4 = df_inputAC %>%
                 filter(dose == "d4", outcome == "infection") %>%
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               delay5 = df_inputAC %>%
                 filter(dose == "d5", outcome == "infection") %>%
                 select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
               # protection against infection
               eta1 = df_inputAC %>% 
                 filter(dose == "d1", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta2 = df_inputAC %>% 
                 filter(dose == "d2", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta3 = df_inputAC %>% 
                 filter(dose == "d3", outcome == "infection") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta4 = df_inputAC %>%
                 filter(dose == "d4", outcome == "infection") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta5 = df_inputAC %>%
                 filter(dose == "d5", outcome == "infection") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               # protection from hospitalisation
               eta_hosp1 = df_inputAC %>% 
                 filter(dose == "d1", outcome == "hospitalisation") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp2 = df_inputAC %>% 
                 filter(dose == "d2", outcome == "hospitalisation") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp3 = df_inputAC %>% 
                 filter(dose == "d3", outcome == "hospitalisation") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp4 = df_inputAC %>%
                 filter(dose == "d4", outcome == "hospitalisation") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_hosp5 = df_inputAC %>%
                 filter(dose == "d5", outcome == "hospitalisation") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               # protection from transmission
               eta_trans1 = df_inputAC %>% 
                 filter(dose == "d1", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans2 = df_inputAC %>% 
                 filter(dose == "d2", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans3 = df_inputAC %>% 
                 filter(dose == "d3", outcome == "transmission") %>% 
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans4 = df_inputAC %>%
                 filter(dose == "d4", outcome == "transmission") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               eta_trans5 = df_inputAC %>%
                 filter(dose == "d5", outcome == "transmission") %>%
                 select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
               p_report = p_reported_by_age,
               contact_mat = april_2017,
               calendar_start_date = as.Date("2020-01-01")
)

# scenarios B & D
paramsBD <- list(N = n_vec,
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
                 omega = 0.02017514,
                 # daily vaccination rate
                 alpha1 = df_inputBD %>% 
                   filter(dose == "d1", outcome == "infection") %>% 
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha2 = df_inputBD %>% 
                   filter(dose == "d2", outcome == "infection") %>% 
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha3 = df_inputBD %>% 
                   filter(dose == "d3", outcome == "infection") %>% 
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha4 = df_inputBD %>%
                   filter(dose == "d4", outcome == "infection") %>%
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 alpha5 = df_inputBD %>%
                   filter(dose == "d5", outcome == "infection") %>%
                   select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
                 # delay to protection
                 delay1 = df_inputBD %>% 
                   filter(dose == "d1", outcome == "infection") %>% 
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay2 = df_inputBD %>% 
                   filter(dose == "d2", outcome == "infection") %>% 
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay3 = df_inputBD %>% 
                   filter(dose == "d3", outcome == "infection") %>% 
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay4 = df_inputBD %>%
                   filter(dose == "d4", outcome == "infection") %>%
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 delay5 = df_inputBD %>%
                   filter(dose == "d5", outcome == "infection") %>%
                   select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
                 # protection against infection
                 eta1 = df_inputBD %>% 
                   filter(dose == "d1", outcome == "infection") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta2 = df_inputBD %>% 
                   filter(dose == "d2", outcome == "infection") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta3 = df_inputBD %>% 
                   filter(dose == "d3", outcome == "infection") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta4 = df_inputBD %>%
                   filter(dose == "d4", outcome == "infection") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta5 = df_inputBD %>%
                   filter(dose == "d5", outcome == "infection") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 # protection from hospitalisation
                 eta_hosp1 = df_inputBD %>% 
                   filter(dose == "d1", outcome == "hospitalisation") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp2 = df_inputBD %>% 
                   filter(dose == "d2", outcome == "hospitalisation") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp3 = df_inputBD %>% 
                   filter(dose == "d3", outcome == "hospitalisation") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp4 = df_inputBD %>%
                   filter(dose == "d4", outcome == "hospitalisation") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_hosp5 = df_inputBD %>%
                   filter(dose == "d5", outcome == "hospitalisation") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 # protection from transmission
                 eta_trans1 = df_inputBD %>% 
                   filter(dose == "d1", outcome == "transmission") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans2 = df_inputBD %>% 
                   filter(dose == "d2", outcome == "transmission") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans3 = df_inputBD %>% 
                   filter(dose == "d3", outcome == "transmission") %>% 
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans4 = df_inputBD %>%
                   filter(dose == "d4", outcome == "transmission") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 eta_trans5 = df_inputBD %>%
                   filter(dose == "d5", outcome == "transmission") %>%
                   select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
                 p_report = p_reported_by_age,
                 contact_mat = april_2017,
                 calendar_start_date = as.Date("2020-01-01")
)

# Specify initial conditions ---------------------------------------------------
init_cond_list <- readRDS("inst/extdata/results/model_fits/initial_conditions.rds")
init_cond <- unlist(init_cond_list[[length(init_cond_list)]])
init_cond[1] <- 872

# ------------------------------------------------------------------------------
# Run forward simulations ------------------------------------------------------
t_start <- init_cond[1]
t_end <- t_start + 365
times <- seq(t_start, t_end, by = 1)
betas <- readRDS("inst/extdata/results/model_fits/beta_draws.rds")
# sample 100 betas from last time window
betas100 <- sample(betas[[length(betas)]]$beta, 100)

# register parallel backend
registerDoParallel(cores=15)

# Scenario A
# Slow waning, Summer booster campaign (increase coverage 4th dose)
scenarioA <- foreach(i = 1:100) %dopar% {
  paramsAC$beta <- betas100[i]
  paramsAC$contact_mat <- april_2017[[i]]
  paramsAC$omega <- wane_8months
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_cond, times, age_struct_seir_ode2, paramsAC, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioA, "/rivm/s/ainsliek/results/scenario_hub/round1/scenarioA.rds")
# Scenario B
# Slow waning, autumn booster campaign (5th dose)
scenarioB <- foreach(i = 1:100) %dopar% {
  paramsBD$beta <- betas100[i]
  paramsBD$contact_mat <- april_2017[[i]]
  paramsBD$omega <- wane_8months
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_cond, times, age_struct_seir_ode2, paramsBD, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioB, "/rivm/s/ainsliek/results/scenario_hub/round1/scenarioB.rds")
# Scenario C
# Fast waning, summer booster campaign (increase coverage of 4th dose)
scenarioC <- foreach(i = 1:100) %dopar% {
  paramsAC$beta <- betas100[i]
  paramsAC$contact_mat <- april_2017[[i]]
  paramsAC$omega <- wane_3months
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_cond, times, age_struct_seir_ode2, paramsAC, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioC, "/rivm/s/ainsliek/results/scenario_hub/round1/scenarioC.rds")
# Scenario D
# Fast waning, autumn booster campaign (5th dose)
scenarioD <- foreach(i = 1:100) %dopar% {
  paramsBD$beta <- betas100[i]
  paramsBD$contact_mat <- april_2017[[i]]
  paramsBD$omega <- wane_3months
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_cond, times, age_struct_seir_ode2, paramsBD, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioD, "/rivm/s/ainsliek/results/scenario_hub/round1/scenarioD.rds")
#-------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Post-process scenario runs ---------------------------------------------------
# Results must be in a csv file that containa only the following columns (in any
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
p_report_vec <- c(rep(as.numeric(paramsAC$p_report),6))

# read in saved output from model runs
# df_scenA <- readRDS("/rivm/s/ainsliek/results/scenario_hub/round1/scenarioA.rds")
df_scenA <- readRDS("C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/ECDC Scenario Modelling Hub/round 1/scenarioA.rds")
sim <- length(df_scenA)
# loop over samples and summarise results
outA <- list()
for(s in 1:sim){
  seir_output <- postprocess_age_struct_model_output2(df_scenA[[s]])
  seir_outcomes <- summarise_results(seir_output, params = paramsAC, t_vec = times) %>%
    mutate(sample = s)
  outA[[s]] <- seir_outcomes
}
dfA <- bind_rows(outA) %>%
  mutate(scenario_id = "A-2022-05-22") %>%
  filter(date <= as.Date("2023-05-20"))

# wrangle Scenario B output ----------------------------------------------------
# read in saved output from model runs
# df_scenB <- readRDS("/rivm/s/ainsliek/results/scenario_hub/round1/scenarioB.rds")
df_scenB <- readRDS("C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/ECDC Scenario Modelling Hub/round 1/scenarioB.rds")
sim <- length(df_scenB)
# loop over samples and summarise results
outB <- list()
for(s in 1:sim){
  seir_output <- postprocess_age_struct_model_output2(df_scenB[[s]])
  seir_outcomes <- summarise_results(seir_output, params = paramsBD, t_vec = times) %>%
    mutate(sample = s)
  outB[[s]] <- seir_outcomes
}
dfB <- bind_rows(outB) %>%
  mutate(scenario_id = "B-2022-05-22") %>%
  filter(date <= as.Date("2023-05-20"))

# wrangle Scenario C output ----------------------------------------------------
# read in saved output from model runs
# df_scenC <- readRDS("/rivm/s/ainsliek/results/scenario_hub/round1/scenarioC.rds")
df_scenC <- readRDS("C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/ECDC Scenario Modelling Hub/round 1/scenarioC.rds")
sim <- length(df_scenC)
# loop over samples and summarise results
outC <- list()
for(s in 1:sim){
  seir_output <- postprocess_age_struct_model_output2(df_scenC[[s]])
  seir_outcomes <- summarise_results(seir_output, params = paramsAC, t_vec = times) %>%
    mutate(sample = s)
  outC[[s]] <- seir_outcomes
}
dfC <- bind_rows(outC) %>%
  mutate(scenario_id = "C-2022-05-22") %>%
  filter(date <= as.Date("2023-05-20"))

# wrangle Scenario D output ----------------------------------------------------
# read in saved output from model runs
# df_scenD <- readRDS("/rivm/s/ainsliek/results/scenario_hub/round1/scenarioD.rds")
df_scenD <- readRDS("C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/ECDC Scenario Modelling Hub/round 1/scenarioD.rds")
sim <- length(df_scenD)
# loop over samples and summarise results
outD <- list()
for(s in 1:sim){
  seir_output <- postprocess_age_struct_model_output2(df_scenD[[s]])
  seir_outcomes <- summarise_results(seir_output, params = paramsBD, t_vec = times) %>%
    mutate(sample = s)
  outD[[s]] <- seir_outcomes
}
dfD <- bind_rows(outD) %>%
  mutate(scenario_id = "D-2022-05-22") %>%
  filter(date <= as.Date("2023-05-20"))

# join all scenarios in a single data frame
df_round1 <- bind_rows(dfA, dfB, dfC, dfD) 

# output for plotting
saveRDS(df_round1, "inst/extdata/results/scenario_hub/2022-05-22-rivm-vacamole.rds")

# put all scenarios together into single data frame and sum over epiweek & 
# age groups
df_round1_sh <- df_round1 %>%
  group_by(scenario_id, sample, epiweek, horizon, target_variable) %>%
  summarise_at(.vars = "value", .funs = "sum") %>%
  ungroup() %>%
  mutate(value = round(value),
         origin_date = as.Date("2022-05-22"),
         target_end_date = as.Date("2023-05-20"),
         location = "NL") %>%
  select(-epiweek)

# output for submission to scenario hub
write_csv(df_round1_sh, "C:/Users/ainsliek/Documents/covid19-scenario-hub-europe/data-processed/RIVM-vacamole/2022-05-22-RIVM-vacamole.csv")
