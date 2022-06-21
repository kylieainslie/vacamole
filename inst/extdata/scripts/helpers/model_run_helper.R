# Load packages/functions and define parameter values to run
# age-structured SEIR model
# Load packages ----------------------------------------------------
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

# Source functions -------------------------------------------------
# source("R/age_struct_seir_ode.R")
# source("R/stochastic_age_struct_seir_ode.R")
# source("R/postprocess_age_struct_model_output.R")
# source("R/choose_contact_matrix.R")
# source("R/get_foi.R")
# source("R/summarise_results.R")
# source("R/convert_vac_schedule.R")
# source("R/calc_ve_w_waning.R")
# source("R/my_rmultinom.R")

# load data ---------------------------------------------------------
# probabilities -----------------------------------------------------
dons_probs <- read_xlsx("inst/extdata/inputs/ProbabilitiesDelays_20210107.xlsx")
p_infection2admission <- dons_probs$P_infection2admission
p_admission2death <- dons_probs$P_admission2death
p_admission2IC <- dons_probs$P_admission2IC
p_IC2hospital <- dons_probs$P_IC2hospital
p_IC2death <- 1 - p_IC2hospital
p_hospital2death <- c(rep(0, 5), 0.01, 0.04, 0.12, 0.29) # (after ICU)
p_reported_by_age <- c(0.29, 0.363, 0.381, 0.545, 0.645, 0.564, 0.365, 0.33, 0.409) # from Jantien
p_reported_all <- 0.428 # from Jantien
p_inf_by_age <- c(0.018, 0.115, 0.156, 0.118, 0.142, 0.199, 0.114, 0.062, 0.054 + 0.023)
p_recovered <- c(
  0.01120993, 0.09663659, 0.24141186, 0.11004723, 0.10677859, 0.11977255,
  0.11904044, 0.11714503, 0.11347191
)

# delays ------------------------------------------------------------
time_symptom2admission <- c(2.29, 5.51, 5.05, 5.66, 6.55, 5.88, 5.69, 5.09, 4.33) # assume same as infectious2admission
time_admission2discharge <- 7.9
time_admission2IC <- 2.28
time_IC2hospital <- 15.6
time_hospital2discharge <- 10.1 # (after ICU)
time_admission2death <- 7
time_IC2death <- 19
time_hospital2death <- 10 # (after ICU)

# age distribution and pop size -------------------------------------
age_dist <- c(
  0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463,
  0.14514332, 0.12092904, 0.08807406, 0.04622194
)
n <- 17407585 # Dutch population size
n_vec <- n * age_dist

# contact matrices --------------------------------------------------
#path <- "inst/extdata/inputs/contact_matrices/converted/"
path <- "/rivm/s/ainsliek/data/contact_matrices/converted/"

april_2017     <- readRDS(paste0(path,"transmission_matrix_april_2017.rds"))
april_2020     <- readRDS(paste0(path,"transmission_matrix_april_2020.rds"))
june_2020      <- readRDS(paste0(path,"transmission_matrix_june_2020.rds"))
september_2020 <- readRDS(paste0(path,"transmission_matrix_september_2020.rds"))
february_2021  <- readRDS(paste0(path,"transmission_matrix_february_2021.rds"))
june_2021      <- readRDS(paste0(path,"transmission_matrix_june_2021.rds"))

# parameter inputs -------------------------------------------------
s <- 0.5
g <- 0.5
r0 <- 5.75

# determine transmission rate (beta) for r0 ------------------------
S <- diag(n_vec - 1)
rho <- as.numeric(eigs(S %*% april_2017, 1)$values)
beta <- (r0 / rho) * g
# check
K <- (1 / g) * beta * S %*% april_2017
as.numeric(eigs(K, 1)$values) # this should be r0

# define state transition rates ------------------------------------
gamma <- (1-p_infection2admission) / 2                 # I -> R
h <- p_infection2admission / time_symptom2admission    # I -> H

i1 <- p_admission2IC / time_admission2IC               # H -> IC
d <- p_admission2death / time_admission2death          # H -> D
r <- (1 - (p_admission2IC + p_admission2death)) / time_admission2discharge# H -> R

i2 <- p_IC2hospital / time_IC2hospital                 # IC -> H_IC
d_ic <- p_IC2death / time_IC2death                     # IC -> D

d_hic <- p_hospital2death / time_hospital2death        # H_IC -> D
r_ic <- (1 - p_hospital2death) / time_hospital2discharge # H_IC -> R

transition_rates <- list(gamma = gamma,
                         h = h,
                         i1 = i1,
                         i2 = i2,
                         d = d,
                         d_ic = d_ic,
                         d_hic = d_hic,
                         r = r,
                         r_ic = r_ic)
saveRDS(transition_rates, "inst/extdata/inputs/transition_rates.rds")
# vaccinations params ----------------------------------------------
delays <- list(
  pfizer = c(14, 7, 7),
  moderna = c(14, 7, 7), 
  astrazeneca = c(14, 7),
  jansen = c(14, 7)
)

ve_inf_list <- list(
  wildtype = list(
    pfizer = c(0.926, 0.948), # from clinical trial
    moderna = c(0.896, 0.941), # from clinical trial
    astrazeneca = c(0.583, 0.621), # from clinical trial
    jansen = c(0.661, 0.661) # from clinical trial
  ),
  alpha = list(
    pfizer = c(0.66, 0.8), # from Pritchard et al. 2021 Nature
    moderna = c(0.66, 0.8), # assumed to be the same as pfizer
    astrazeneca = c(0.61, 0.79), # from Pritchard et al. 2021 Nature 
    jansen = c(0.767, 0.767) # from Corchado-Garcia et al. 2021 medRxiv (need to check if this is against alpha!)
  ),
  delta = list( # from Dutch data sources
    pfizer = c(0.57, 0.69, 0.93), 
    moderna = c(0.66, 0.82, 0.93), 
    astrazeneca = c(0.41, 0.54), 
    jansen = c(0.5, 0.5)
  ),
  omicron = list( # from https://www.medrxiv.org/content/10.1101/2022.02.06.22270457v2
    pfizer = c(0, 0.33, 0.68, 0.65, 0.65),  # Grewal et al. 2022 (4th dose) - assuming same for 5th dose
    moderna = c(0, 0.33, 0.68, 0.65, 0.65), # Grewal et al. 2022 (4th dose) - assuming same for 5th dose
    astrazeneca = c(0, 0.33), 
    jansen = c(0, 0.33)
  )
)


ve_trans_list <- list(
  wildtype = list( # same as alpha !!!
    pfizer = c(0.26, 0.70),      
    moderna = c(0.51, 0.88),     
    astrazeneca = c(0.15, 0.58),
    jansen = c(0.77)
  ),
  alpha = list(  # de Gier et al.
    pfizer = c(0.26, 0.70),      
    moderna = c(0.51, 0.88),     
    astrazeneca = c(0.15, 0.58),
    jansen = c(0.77)
  ),
  delta = list( # de Gier et al. (updated)
    pfizer = c(0.46, 0.52),      
    moderna = c(0.66, 0.24),     
    astrazeneca = c(0, 0.25),
    jansen = c(0.42) 
  ),
  omicron = list( # MADE-UP VALUES
    pfizer = c(0.25, 0.33, 0.4, 0.4, 0.4),      
    moderna = c(0.25, 0.33, 0.4, 0.4, 0.4),     
    astrazeneca = c(0, 0.25),
    jansen = c(0.25, 0.33) 
  )
)


ve_hosp_list <- list(
  wildtype = list( # same as alpha variant!!!
    pfizer = c(0.81, 0.95),      # Dutch data
    moderna = c(0.81, 0.95),     # assumed same as pfizer because Dutch estimates were weird
    astrazeneca = c(0.83, 0.95), # Dutch data
    jansen = c(0.85)             # from RIVM website: https://www.rivm.nl/en/covid-19-vaccination/vaccines/efficacy-and-protection
  ),
  alpha = list(
    pfizer = c(0.81, 0.95),      # Dutch data
    moderna = c(0.81, 0.95),     # assumed same as pfizer because Dutch estimates were weird
    astrazeneca = c(0.83, 0.95), # Dutch data
    jansen = c(0.85)             # from RIVM website: https://www.rivm.nl/en/covid-19-vaccination/vaccines/efficacy-and-protection
  ),
  delta = list(
   pfizer = c(0.89, 0.96, 0.98),      # from Brechje (pre-print)
   moderna = c(0.95, 0.85, 0.98),     # from Brechje (pre-print)
   astrazeneca = c(0.88, 0.94),       # from Brechje (pre-print)
   jansen = c(0.92)                   # from Brechje (pre-print)
  ),
  omicron = list( # from https://www.rivm.nl/documenten/effectiviteit-van-covid-19-vaccinatie-tegen-ziekenhuis-en-intensive-care-opname-in-8
    pfizer = c(0, 0.56, 0.88, 0.92, 0.92),   # Grewal et al. 2022 (4th dose) - assuming same for 5th dose
    moderna = c(0, 0.49, 0.82, 0.92, 0.92),  # Grewal et al. 2022 (4th dose) - assuming same for 5th dose
    astrazeneca = c(0, 0.39),
    jansen = c(0.73, 0.73)
  )
)
# hospitalisations multiplier
# calculated as (1-ve_hosp)/(1-ve_inf)
calc_mult_fun <- function(ve_hosp,ve_inf){
  h_multiplier <- list()
  
  for (i in 1:length(ve_hosp)){
    
    h_multiplier[[i]] <- list(
      pfizer = (1-ve_hosp[[i]]$pfizer)/(1-ve_inf[[i]]$pfizer),
      moderna = (1-ve_hosp[[i]]$moderna)/(1-ve_inf[[i]]$moderna),
      astrazeneca = (1-ve_hosp[[i]]$astrazeneca)/(1-ve_inf[[i]]$astrazeneca),
      jansen = (1-ve_hosp[[i]]$jansen)/(1-ve_inf[[i]]$jansen)
    )
  }
  names(h_multiplier) <- names(ve_hosp)
  return(h_multiplier)
}

h_mult_list <- calc_mult_fun(ve_hosp = ve_hosp_list, ve_inf = ve_inf_list)

ve_list <- list(delays = delays,
                ve_inf = ve_inf_list, 
                ve_hosp = h_mult_list,
                ve_trans = ve_trans_list)
saveRDS(ve_list, "inst/extdata/inputs/ve_params.rds")

# read in vac schedules --------------------------------------------
vac_path <- "C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/manuscripts/impact_vac/data/vaccination_scenarios/"
# basis_12plus <- read_csv(paste0(vac_path,"vac_schedule_12plus.csv")) %>%
#   select(-starts_with("X"))

# basis_18plus <- read_csv("inst/extdata/data/vaccination_scenarios/Cum_upt20210701 BASIS 75% in 18+ KA.csv") %>%
#   select(-starts_with("X"))

vac_sched <- read_csv("inst/extdata/inputs/vac_schedule_real_w_booster.csv") %>%
  rename_with(~ gsub("B1", "d3", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("B2", "d4", .x, fixed = TRUE)) %>%
  select(-starts_with("X"))

basis1 <- convert_vac_schedule(
  vac_schedule = vac_sched,
  ve = ve,
  hosp_multiplier = h_multiplier,
  delay = delays,
  ve_trans = ve_trans,
  wane = FALSE,
  add_child_vac = FALSE,
  add_extra_dates = TRUE,
  extra_end_date = "2022-03-31"
)

# baseline parameter inputs
params <- list(beta = 0.0003934816 * 2 ,  # transmission rate
               beta1 = 0.14,              # amplitude of seasonal forcing
               gamma = g,                 # 1/gamma = infectious period
               sigma = s,                 # 1/sigma = latent period
               epsilon = 0.01,            # import case
               N = n_vec,                 # Population (no need to change)
               h = h,                     # Rate from infection to hospital admission/ time from infection to hosp admission
               i1 = i1,
               i2 = i2,
               d = d, 
               d_ic = d_ic,
               d_hic = d_hic,
               r = r,
               r_ic = r_ic,
               p_report = 1/3, #p_reported_by_age,
               c_start = june_2021,
               c_lockdown = february_2021,
               c_relaxed = june_2020,
               c_very_relaxed = june_2021,
               c_normal = april_2017,
               keep_cm_fixed = FALSE,
               vac_inputs = basis1,
               use_cases = TRUE,                           # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),         # somewhat arbitrary cut-off ***need to check if realistic
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 100000/100000 * sum(n_vec),        #35.7  # 20 for IC admissions
               no_vac = FALSE,
               t_calendar_start = yday(as.Date("2020-01-01")),   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
               beta_change = NULL 
)

# initial values
empty_state <- c(rep(0, 9))
init <- c(
  t = 0,
  S = c(n_vec[1:4], n_vec[5]-1, n_vec[6:9]),
  Shold_1d = empty_state,
  Sv_1d = empty_state,
  Shold_2d = empty_state,
  Sv_2d = empty_state,
  E = empty_state,
  Ev_1d = empty_state,
  Ev_2d = empty_state,
  I = c(rep(0,4),1,rep(0,4)),
  Iv_1d = empty_state,
  Iv_2d = empty_state,
  H = empty_state,
  Hv_1d = empty_state,
  Hv_2d = empty_state,
  H_IC = empty_state,
  H_ICv_1d = empty_state,
  H_ICv_2d = empty_state,
  IC = empty_state,
  ICv_1d = empty_state,
  ICv_2d = empty_state,
  D = empty_state,
  R = empty_state,
  Rv_1d = empty_state,
  Rv_2d = empty_state
)

# times vector
times <- seq(0, nrow(vac_sched)-1, by = 1)
