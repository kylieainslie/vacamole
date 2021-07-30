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
library(lazymcmc)
library(lubridate)

# Source functions -------------------------------------------------
source("R/age_struct_seir_ode.R")
source("R/stochastic_age_struct_seir_ode.R")
source("R/postprocess_age_struct_model_output.R")
source("R/choose_contact_matrix.R")
source("R/get_foi.R")
source("R/summarise_results.R")
source("R/convert_vac_schedule.R")
source("R/calc_ve_w_waning.R")
source("R/my_rmultinom.R")

# load data ---------------------------------------------------------
# probabilities -----------------------------------------------------
dons_probs <- read_xlsx("inst/extdata/data/ProbabilitiesDelays_20210107.xlsx")
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
baseline_2017 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_baseline_2017.rds")
april_2020 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_april_2020.rds")
june_2020 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_june_2020.rds")
september_2020 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_september_2020.rds")
february_2021 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_february_2021.rds")
june_2021 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_june_2021.rds")

# parameter inputs -------------------------------------------------
s <- 0.5
g <- 0.5
r0 <- 2.3

# determine transmission rate (beta) for r0 ------------------------
S <- diag(n_vec - 1)
rho <- as.numeric(eigs(S %*% baseline_2017$mean, 1)$values)
beta <- (r0 / rho) * g
# check
K <- (1 / g) * beta * S %*% baseline_2017$mean
as.numeric(eigs(K, 1)$values) # this should be r0

days <- seq(yday(as.Date("2021-01-31")), yday(as.Date("2021-12-31")), by = 1)
breakpoints <- yday(c(
  as.Date("2021-02-28"),
  as.Date("2021-03-21"),
  as.Date("2021-04-01"),
  as.Date("2021-04-08"),
  as.Date("2021-04-15"),
  as.Date("2021-04-22"),
  as.Date("2021-04-30"),
  as.Date("2021-05-25")
))
beta_mle_vals <- c(0.00045, 0.00038, 0.00041, 0.00052, 0.00057, 0.00061, 0.00062, 0.00061)
betas <- c(
  rep(beta_mle_vals[1], length(days[1]:breakpoints[1])),
  rep(beta_mle_vals[2], length(breakpoints[1]:breakpoints[2]) - 1),
  rep(beta_mle_vals[3], length(breakpoints[2]:breakpoints[3]) - 1),
  rep(beta_mle_vals[4], length(breakpoints[3]:breakpoints[4]) - 1),
  rep(beta_mle_vals[5], length(breakpoints[4]:breakpoints[5]) - 1),
  rep(beta_mle_vals[6], length(breakpoints[5]:breakpoints[6]) - 1),
  rep(beta_mle_vals[7], length(breakpoints[6]:breakpoints[7]) - 1),
  rep(beta_mle_vals[8], length(breakpoints[7]:days[length(days)]) - 1)
)
beta_t_ <- betas * (1 + 0.15 * cos(2 * pi * days / 365.24))
R0 <- (beta_t_ / g) * rho
R0_dat <- date <- data.frame(
  date = seq.Date(as.Date("2021-01-31"), as.Date("2021-12-31"), by = 1),
  R0 = R0,
  roll_mean_R0 = zoo::rollmean(R0, k = 7, fill = NA)
)
r0_plot <- ggplot(data = R0_dat, aes(x = date, y = roll_mean_R0)) +
  geom_line() +
  labs(y = "Basic reproduction number (R0)", x = "Date") +
  theme( # legend.position = "bottom",
    panel.background = element_blank(),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )
r0_plot

# ggsave("inst/extdata/results/r0_plot_16june.jpg",
#   plot = r0_plot,
#   # height = 8,
#   # width = 12,
#   dpi = 300
# )

# define state transition rates ------------------------------------
h <- p_infection2admission / time_symptom2admission
i1 <- p_admission2IC / time_admission2IC
i2 <- p_IC2hospital / time_IC2hospital
d <- p_admission2death / time_admission2death
d_ic <- p_IC2death / time_IC2death
d_hic <- p_hospital2death / time_hospital2death
r <- (1 - p_admission2death) / time_admission2discharge
r_ic <- (1 - p_IC2death) / time_hospital2discharge

# vaccinations params ----------------------------------------------
ve <- list(
  pfizer = c(0.7, 0.85), # from SIREN study, estimates from clinical trial: 0.926, 0.948
  moderna = c(0.896, 0.941), # from clinical trial
  astrazeneca = c(0.583, 0.621), # from clinical trial
  jansen = c(0.661)
)

delays <- list(
  pfizer = c(21, 7),
  moderna = c(14, 14),
  astrazeneca = c(21, 14),
  jansen = c(14)
)

ve_sa20 <- list(
  pfizer = c(0.741, 0.758),
  moderna = c(0.717, 0.753),
  astrazeneca = c(0.466, 0.497),
  jansen = c(0.529)
)

ve_sa50 <- list(
  pfizer = c(0.463, 0.474),
  moderna = c(0.448, 0.471),
  astrazeneca = c(0.269, 0.311),
  jansen = c(0.331)
)

ve_vasileiou <- list(
  pfizer = c(0.645, 0.85),
  moderna = c(0.896, 0.941),
  astrazeneca = c(0.583, 0.621),
  jansen = c(0.661)
)
delays_vasileiou <- list(
  pfizer = c(7, 7),
  moderna = c(14, 14),
  astrazeneca = c(21, 14),
  jansen = c(14)
)

ve_trans <- list(
  pfizer = c(0.3, 0.54), # from Shah et al. 2021
  moderna = c(0.3, 0.54), # assumed to be the same as pfizer
  astrazeneca = c(0, 0), # no available data
  jansen = c(0)
) # no available data

# hospitalisations multiplier
# calculated as (1-ve_hosp)/(1-ve)
h_multiplier <- list(
  pfizer = c(0, 0),
  moderna = c(0, 0),
  astrazeneca = c(0.384, 0.384),
  jansen = c(0)
)
# initial states ---------------------------------------------------
# Jacco's suggested way to determine initial conditions
init_states_dat <- data.frame(
  age_group = c(
    "0-9", "10-19", "20-29", "30-39", "40-49",
    "50-59", "60-69", "70-79", "80+"
  ),
  n = n_vec,
  # from Scott
  n_recovered = c(
    30720, 397100, 642600, 419000, 412200, 505900,
    349100, 206800, 115900 + 33200
  ),
  # from sitrep for 26 januari tot 2 februari:
  # https://www.rivm.nl/coronavirus-covid-19/actueel/wekelijkse-update-epidemiologische-situatie-covid-19-in-nederland)
  n_cases = c(
    835, 2851, 4591, 3854, 3925, 5191, 3216, 1819,
    1376 + 485
  ),
  # from NICE data (n_hosp/n_ic refers to occupancy on 1 Feb 2021)
  n_hosp = c(2, 1, 8, 19, 29, 76, 142, 159, 186),
  n_ic = c(0, 2, 6, 9, 25, 83, 181, 150, 12),
  # from NICE data: people with length of stay >= 9 days
  n_hosp_after_ic = c(2, 2, 9, 15, 45, 158, 321, 392, 266)
) %>%
  mutate(
    n_infections = n_cases * 2.2, # calibrated to match osiris data
    init_E = n_infections * (2 / 7),
    init_I = n_infections * (2 / 7),
    init_S = n - n_recovered - init_E - init_I - n_hosp - n_ic - n_hosp_after_ic
  )

# determine transmission rate for reff ------------------------------
reff <- 1.04 # from RIVM open data for 1 Feb 2021 (midpoint between
#              0.94 (wt) and 1.13 (UK variant))
S2 <- diag(init_states_dat$init_S)
rho2 <- as.numeric(eigs(S2 %*% t5, 1)$values)
beta2 <- reff / rho2 * g
# check
B <- t5
K2 <- beta2 * (1 / g) * S2 %*% B
as.numeric(eigs(K2, 1)$values) # this should be r0

# callibrate distribution of cases across age groups ----------------
w <- eigen(K2)$vectors[, 1] / sum(eigen(K2)$vectors[, 1]) # should match dist_cases (I think!)
x <- init_states_dat$n_cases / sum(init_states_dat$n_cases)
A <- diag(x / w)
B_prime <- A %*% B
rho2_prime <- as.numeric(eigs(S2 %*% B_prime, 1)$values)
beta2_prime <- reff / rho2_prime * g
K2_prime <- beta2_prime * (1 / g) * S2 %*% B_prime
dom_eig_vec <- eigen(K2_prime)$vectors[, 1]
w_prime <- dom_eig_vec / sum(dom_eig_vec)
as.numeric(eigs(K2_prime, 1)$values)

# Specify initial values -------------------------------------------
empty_state <- c(rep(0, 9))
init <- c(
  t = 0,
  S = init_states_dat$init_S,
  Shold_1d = empty_state,
  Sv_1d = empty_state,
  Shold_2d = empty_state,
  Sv_2d = empty_state,
  E = init_states_dat$init_E,
  Ev_1d = empty_state,
  Ev_2d = empty_state,
  I = init_states_dat$init_I,
  Iv_1d = empty_state,
  Iv_2d = empty_state,
  H = init_states_dat$n_hosp,
  Hv_1d = empty_state,
  Hv_2d = empty_state,
  H_IC = init_states_dat$n_hosp_after_ic,
  H_ICv_1d = empty_state,
  H_ICv_2d = empty_state,
  IC = init_states_dat$n_ic,
  ICv_1d = empty_state,
  ICv_2d = empty_state,
  D = empty_state,
  R = init_states_dat$n_recovered,
  Rv_1d = empty_state,
  Rv_2d = empty_state
)

# read in vac schedules --------------------------------------------
basis <- read_csv("inst/extdata/data/Cum_upt20210603 BASIS.csv") %>%
  select(-starts_with("X"))

# no childhood vaccination, waning
basis1 <- convert_vac_schedule(
  vac_schedule = basis,
  ve = ve,
  hosp_multiplier = h_multiplier,
  delay = delays,
  ve_trans = ve_trans,
  wane = TRUE,
  add_child_vac = FALSE,
  add_extra_dates = TRUE,
  extra_end_date = "2022-03-31"
)

# childhood vaccination, waning
basis1_child_vac <- convert_vac_schedule(
  vac_schedule = basis,
  ve = ve,
  hosp_multiplier = h_multiplier,
  delay = delays,
  ve_trans = ve_trans,
  wane = TRUE,
  add_child_vac = TRUE,
  child_vac_coverage = 0.7,
  child_doses_per_day = 50000,
  child_vac_start_date = "2021-07-15",
  add_extra_dates = TRUE,
  extra_end_date = "2022-03-31"
)

# no childhood vaccination, no waning
basis1_no_wane <- convert_vac_schedule(
  vac_schedule = basis,
  ve = ve,
  hosp_multiplier = h_multiplier,
  delay = delays,
  ve_trans = ve_trans,
  wane = FALSE,
  add_child_vac = FALSE,
  add_extra_dates = TRUE,
  extra_end_date = "2022-03-31"
)

# childhood vaccination, no waning
basis1_no_wane_child_vac <- convert_vac_schedule(
  vac_schedule = basis,
  ve = ve,
  hosp_multiplier = h_multiplier,
  delay = delays,
  ve_trans = ve_trans,
  wane = FALSE,
  add_child_vac = TRUE,
  child_vac_coverage = 0.7,
  child_doses_per_day = 50000,
  child_vac_start_date = "2021-07-15",
  add_extra_dates = TRUE,
  extra_end_date = "2022-03-31"
)


# load model fits --------------------------------------------------
daily_cases_from_fit <- readRDS("inst/extdata/results/model_fit_df_2021-05-25.rds")
# mle_betas <- read_csv("inst/extdata/results/mle_betas_2021-04-26.csv")
initial_conditions <- readRDS("inst/extdata/results/init_conditions_2021-05-25.rds")
# osiris_dat <- readRDS("inst/extdata/data/Osiris_Data_20210505_1034.rds")
