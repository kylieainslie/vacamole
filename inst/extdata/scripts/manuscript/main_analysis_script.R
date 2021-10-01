# ------------------------------------------------------------------
# Forward simulations to determine impacts of vaccinating
# 12-17 year olds for manuscript
# ------------------------------------------------------------------

# Data and model parameters are loaded/defined in the script 
# inst/extdata/scripts/model_run_helper.R
source("inst/extdata/scripts/helpers/model_run_helper.R")
source("R/forward_sim_func_wrap.R")
#library(foreach)
# read in vac schedules --------------------------------------------
basis_12plus <- read_csv("inst/extdata/data/vaccination_scenarios/Cum_upt20210701 BASIS 75% in 12+ KA.csv") %>%
  select(-starts_with("X"))

basis_18plus <- read_csv("inst/extdata/data/vaccination_scenarios/Cum_upt20210701 BASIS 75% in 18+ KA.csv") %>%
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
           june_2021 = june_2021
           )
# for sd, use sqrt(var) of beta from fits to last time window before
# forward simulations
my_sd <- 0.0000945 
sigma <- matrix(c(0.0000945, -0.0000261,
                  -0.0000261, 0.00176), nrow = 2)
#test <- mvrnorm(n = 200, mu = c(0.0003934816 * 2, 0.01), Sigma = sigma)
betas_for_delta_period <- truncnorm::rtruncnorm(n = 200, 
                                                a = 0.000001, 
                                                b = Inf, 
                                                mean = 0.0003934816 * 2, 
                                                sd = my_sd) # R0 = 4.6
betas_for_alpha_period <- truncnorm::rtruncnorm(n = 200, 
                                                a = 0.000001,
                                                b = Inf,
                                                mean = 0.0005902224, 
                                                sd = my_sd)   # R0 = 3.45
#betas_upper <- rnorm(n = 200, mean = 0.0009837041, sd = my_sd)   # R0 = 5.75
  #0.001539711    # R0 = 9

# -------------------------------------------------------------------
# run models --------------------------------------------------------
# -------------------------------------------------------------------
todays_date <- Sys.Date()

# register clusters for parallel computing
# cl <- parallel::makeCluster(6)
# doParallel::registerDoParallel(cl)
# parallel::stopCluster(cl)
# no waning ---------------------------------------------------------
basis_12plus1 <- convert_vac_schedule(
  vac_schedule = basis_12plus,
  ve = ve,
  hosp_multiplier = h_multiplier,
  delay = delays,
  ve_trans = ve_trans,
  wane = FALSE,
  before_feb = FALSE,
  add_child_vac = FALSE,
  add_extra_dates = TRUE,
  extra_end_date = "2022-03-31"
)

basis_18plus1 <- convert_vac_schedule(
  vac_schedule = basis_18plus,
  ve = ve,
  hosp_multiplier = h_multiplier,
  delay = delays,
  ve_trans = ve_trans,
  wane = FALSE,
  before_feb = FALSE,
  add_child_vac = FALSE,
  add_extra_dates = TRUE,
  extra_end_date = "2022-03-31"
)

# 12+ alpha
forward_sim_func_wrap(start_date = "2021-06-22",
                      end_date = "2021-03-31",
                      init_cond = init_cond_22june2021,
                      beta_m = beta_mles[index,1],
                      vac_inputs = basis_12plus1,
                      beta_c = betas_for_alpha_period,
                      beta_d = beta_draws[[index]][,1],
                      t_normal = NULL,
                      contact_matrices = cm,
                      tag = paste0("results_12plus_mle_beta_",todays_date)
                      )

# 12+ delta
forward_sim_func_wrap(start_date = "2021-06-22",
                      end_date = "2021-03-31",
                      init_cond = init_cond_22june2021,
                      beta_m = beta_mles[index,1],
                      vac_inputs = basis_12plus1,
                      beta_c = betas_for_delta_period,
                      beta_d = beta_draws[[index]][,1],
                      t_normal = NULL,
                      contact_matrices = cm,
                      tag = paste0("results_12plus_lower_beta_",todays_date)
)

# # 12+ upper
# forward_sim_func_wrap(start_date = "2021-06-22",
#                       end_date = "2021-03-31",
#                       init_cond = init_cond_22june2021,
#                       beta_m = beta_mles[index,1],
#                       vac_inputs = basis_12plus1,
#                       beta_c = beta_delta_upper,
#                       beta_d = beta_draws[[index]][,1],
#                       t_normal = NULL,
#                       contact_matrices = cm,
#                       tag = paste0("results_12plus_upper_beta_",todays_date)
# )

# 18+ alpha
forward_sim_func_wrap(start_date = "2021-06-22",
                      end_date = "2021-03-31",
                      init_cond = init_cond_22june2021,
                      beta_m = beta_mles[index,1],
                      vac_inputs = basis_18plus1,
                      beta_c = betas_for_alpha_period,
                      beta_d = beta_draws[[index]][,1],
                      t_normal = NULL,
                      contact_matrices = cm,
                      tag = paste0("results_18plus_mle_beta_",todays_date)
)
# 18+ delta
forward_sim_func_wrap(start_date = "2021-06-22",
                      end_date = "2021-03-31",
                      init_cond = init_cond_22june2021,
                      beta_m = beta_mles[index,1],
                      vac_inputs = basis_18plus1,
                      beta_c = betas_for_delta_period,
                      beta_d = beta_draws[[index]][,1],
                      t_normal = NULL,
                      contact_matrices = cm,
                      tag = paste0("results_18plus_lower_beta_",todays_date)
)

# # 18+ upper
# forward_sim_func_wrap(start_date = "2021-06-22",
#                       end_date = "2021-03-31",
#                       init_cond = init_cond_22june2021,
#                       beta_m = beta_mles[index,1],
#                       vac_inputs = basis_18plus1,
#                       beta_c = beta_delta_upper,
#                       beta_d = beta_draws[[index]][,1],
#                       t_normal = NULL,
#                       contact_matrices = cm,
#                       tag = paste0("results_18plus_upper_beta_",todays_date)
# )

# waning ------------------------------------------------------------
basis_12plus1_wane <- convert_vac_schedule(
  vac_schedule = basis_12plus,
  ve = ve,
  hosp_multiplier = h_multiplier,
  delay = delays,
  ve_trans = ve_trans,
  wane = TRUE,
  before_feb = FALSE,
  add_child_vac = FALSE,
  add_extra_dates = TRUE,
  extra_end_date = "2022-03-31"
)

basis_18plus1_wane <- convert_vac_schedule(
  vac_schedule = basis_18plus,
  ve = ve,
  hosp_multiplier = h_multiplier,
  delay = delays,
  ve_trans = ve_trans,
  wane = TRUE,
  before_feb = FALSE,
  add_child_vac = FALSE,
  add_extra_dates = TRUE,
  extra_end_date = "2022-03-31"
)

# 12+ alpha
forward_sim_func_wrap(start_date = "2021-06-22",
                      end_date = "2021-03-31",
                      init_cond = init_cond_22june2021,
                      beta_m = beta_mles[index,1],
                      vac_inputs = basis_12plus1_wane,
                      beta_c = betas_for_alpha_period,
                      beta_d = beta_draws[[index]][,1],
                      t_normal = NULL,
                      contact_matrices = cm,
                      tag = paste0("results_12plus_wane_mle_beta_",todays_date)
)
# 12+ delta
forward_sim_func_wrap(start_date = "2021-06-22",
                      end_date = "2021-03-31",
                      init_cond = init_cond_22june2021,
                      beta_m = beta_mles[index,1],
                      vac_inputs = basis_12plus1_wane,
                      beta_c = betas_for_delta_period,
                      beta_d = beta_draws[[index]][,1],
                      t_normal = NULL,
                      contact_matrices = cm,
                      tag = paste0("results_12plus_wane_lower_beta_",todays_date)
)

# # 12+ upper
# forward_sim_func_wrap(start_date = "2021-06-22",
#                       end_date = "2021-03-31",
#                       init_cond = init_cond_22june2021,
#                       beta_m = beta_mles[index,1],
#                       vac_inputs = basis_12plus1_wane,
#                       beta_c = beta_delta_upper,
#                       beta_d = beta_draws[[index]][,1],
#                       t_normal = NULL,
#                       contact_matrices = cm,
#                       tag = paste0("results_12plus_wane_upper_beta_",todays_date)
# )

# 18+ alpha
forward_sim_func_wrap(start_date = "2021-06-22",
                      end_date = "2021-03-31",
                      init_cond = init_cond_22june2021,
                      beta_m = beta_mles[index,1],
                      vac_inputs = basis_18plus1_wane,
                      beta_c = betas_for_alpha_period,
                      beta_d = beta_draws[[index]][,1],
                      t_normal = NULL,
                      contact_matrices = cm,
                      tag = paste0("results_18plus_wane_mle_beta_",todays_date)
)
# 18+ delta
forward_sim_func_wrap(start_date = "2021-06-22",
                      end_date = "2021-03-31",
                      init_cond = init_cond_22june2021,
                      beta_m = beta_mles[index,1],
                      vac_inputs = basis_18plus1_wane,
                      beta_c = betas_for_delta_period,
                      beta_d = beta_draws[[index]][,1],
                      t_normal = NULL,
                      contact_matrices = cm,
                      tag = paste0("results_18plus_wane_lower_beta_",todays_date)
)

# # 18+ upper
# forward_sim_func_wrap(start_date = "2021-06-22",
#                       end_date = "2021-03-31",
#                       init_cond = init_cond_22june2021,
#                       beta_m = beta_mles[index,1],
#                       vac_inputs = basis_18plus1_wane,
#                       beta_c = beta_delta_upper,
#                       beta_d = beta_draws[[index]][,1],
#                       t_normal = NULL,
#                       contact_matrices = cm,
#                       tag = paste0("results_18plus_wane_upper_beta_",todays_date)
# )
