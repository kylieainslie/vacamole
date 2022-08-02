# ------------------------------------------------------------------
# Forward simulations to determine impacts of vaccinating
# 12-17 year olds for manuscript
# ------------------------------------------------------------------

# Data and model parameters are loaded/defined in the script 
# inst/extdata/scripts/model_run_helper.R
source("inst/extdata/scripts/helpers/model_run_helper.R")
source("R/forward_sim_func_wrap.R")

# ------------------------------------------------------------------
# read in vac schedules --------------------------------------------
# ------------------------------------------------------------------
basis_18plus <- read_csv("inst/extdata/data/vaccination_scenarios/vac_schedule_18plus.csv") %>%
  select(-starts_with("X"))

basis_12plus_jun <- read_csv("inst/extdata/data/vaccination_scenarios/vac_schedule_12plus.csv") %>%
  select(-starts_with("X"))

# model wrapper function inputs ------------------------------------
date_of_fit <- "2021-10-01"

output_from_model_fit <- readRDS(paste0("inst/extdata/results/model_fits/output_from_fits_", date_of_fit, ".rds"))
init_cond_22june2021 <- unlist(lapply(unname(output_from_model_fit$`end_date_2021-06-22`), tail,1))
beta_mles <- data.frame(beta = readRDS(paste0("inst/extdata/results/model_fits/mles_from_fits_", date_of_fit,".rds"))) %>%
  mutate(end_date = names(output_from_model_fit))
beta_draws <- readRDS(paste0("inst/extdata/results/model_fits/beta_draws_from_fits_", date_of_fit, ".rds"))

index <- which(beta_mles$end_date == "end_date_2021-06-22")
cm <- list(baseline_2017 = baseline_2017,
           april_2020 = april_2020,
           june_2020 = june_2020,
           september_2020 = september_2020,
           february_2021 = february_2021,
           june_2021 = june_2021)
# for sd, use sqrt(var) of beta from fits to last time window before
# forward simulations (variance is for Rt, so need to convert to beta)
my_sd <- sqrt(0.0000945) 
S_diag <- diag(init[c(2:10)])
rho <- as.numeric(eigs(S_diag %*% params$c_normal, 1)$values)

# R draws
rt_for_delta_period <- rnorm(n = 200, mean = 4.6, sd = my_sd)  # beta = 0.0003934816 * 2 = 0.0007869632
rt_for_alpha_period <- rnorm(n = 200, mean = 3.45, sd = my_sd) # beta = 0.0005902224

# convert to beta
betas_for_delta_period <- c(0.0007869632, (rt_for_delta_period / rho) * params$gamma)
betas_for_alpha_period <- c(0.0005902224, (rt_for_alpha_period / rho) * params$gamma)

# -------------------------------------------------------------------
# run models --------------------------------------------------------
# -------------------------------------------------------------------
todays_date <- Sys.Date()


# no waning ---------------------------------------------------------
# with vaccination of 5-11 year olds
basis_5plus1_delta <- convert_vac_schedule(vac_schedule = basis_12plus_jun,
  ve = ve_delta, hosp_multiplier = h_multiplier_delta, ve_trans = ve_trans_delta, delay = delays,
  wane = FALSE, add_child_vac = TRUE, child_vac_start_date = "2021-10-01",
  add_extra_dates = TRUE, extra_end_date = "2022-03-31")

basis_5plus1_delta_wane <- convert_vac_schedule(vac_schedule = basis_12plus_jun,
  ve = ve_delta, hosp_multiplier = h_multiplier_delta, ve_trans = ve_trans_delta, delay = delays,
  wane = TRUE, add_child_vac = TRUE, child_vac_start_date = "2021-10-01",
  add_extra_dates = TRUE, extra_end_date = "2022-03-31")

# delta ---------------------------------------------------------------
# 5+, no waning
forward_sim_func_wrap(start_date = "2021-06-22",
                      end_date = "2021-03-31",
                      init_cond = init_cond_22june2021,
                      beta_m = beta_mles[index,1],
                      vac_inputs = basis_5plus1_delta,
                      beta_c = betas_for_delta_period,
                      beta_d = beta_draws[[index]][,1],
                      t_normal = yday(as.Date("2021-11-01")) + 365,
                      contact_matrices = cm,
                      tag = paste0("results_5plus_delta_",todays_date)
)
# 5+, waning
forward_sim_func_wrap(start_date = "2021-06-22",
                      end_date = "2021-03-31",
                      init_cond = init_cond_22june2021,
                      beta_m = beta_mles[index,1],
                      vac_inputs = basis_5plus1_delta_wane,
                      beta_c = betas_for_delta_period,
                      beta_d = beta_draws[[index]][,1],
                      t_normal = yday(as.Date("2021-11-01")) + 365,
                      contact_matrices = cm,
                      tag = paste0("results_5plus_delta_wane",todays_date)
)
