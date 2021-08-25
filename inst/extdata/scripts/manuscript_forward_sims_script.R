# Forward simulations to determine impacts of vaccinating
# 12-17 year olds for manuscript
# ------------------------------------------------------------------

# Data and model parameters are loaded/defined in the script 
# inst/extdata/scripts/model_run_helper.R
source("inst/extdata/scripts/model_run_helper.R")

# read in vac schedules --------------------------------------------
basis_12plus <- read_csv("inst/extdata/data/vaccination_scenarios/Cum_upt20210701 BASIS 75% in 12+ KA.csv") %>%
  select(-starts_with("X"))

basis_18plus <- read_csv("inst/extdata/data/vaccination_scenarios/Cum_upt20210701 BASIS 75% in 18+ KA.csv") %>%
  select(-starts_with("X"))

# no childhood vaccination, waning
vac_sched <- basis_18plus

basis1 <- convert_vac_schedule(
  vac_schedule = vac_sched,
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


# load initial conditions from model fits --------------------------
last_date_of_fit <- "2021-07-27"

output_from_model_fit <- readRDS(paste0("inst/extdata/results/model_fits/output_from_fits_", last_date_in_osiris, ".rds"))
init_cond_22june2021 <- unlist(lapply(unname(output_from_model_fit$`end_date_2021-06-22`), tail,1))
beta_mles <- data.frame(beta = readRDS(paste0("inst/extdata/results/model_fits/mles_from_fits_",last_date_of_fit,".rds"))) %>%
  mutate(end_date = names(output_from_model_fit))
parameter_draws <-  readRDS(paste0(path, "parameter_draws_from_fits_", last_date_in_osiris, ".rds"))

index <- which(beta_mles$end_date == "end_date_2021-06-22")

tag <- "df_basis_18plus_mle_beta_23aug"
start_date <- "2021-06-22"
end_date <- "2022-03-30"
init_cond <- init_cond_22june2021 # without first time point
beta_m <- beta_mles[index,1]
vac_inputs <- basis1
beta_c <- 0.0003934816 * 2
  # mle = 0.0003934816 * 2 (R0 = 4.6); 
  # lower = 0.0005902224 (R0 = 3.45); 
  # upper = 0.001539711 (R0 = 9)
contact_matrices <- list(baseline_2017 = baseline_2017,
                         april_2020 = april_2020,
                         june_2020 = june_2020,
                         september_2020 = september_2020,
                         february_2021 = february_2021,
                         june_2021 = june_2021
                         )
# ------------------------------------------------------------------
# call forward simulation function wrapper -------------------------
#-------------------------------------------------------------------

