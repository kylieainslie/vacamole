# ------------------------------------------------------------------
# compare running time of forward simulation wrapper functions
# ------------------------------------------------------------------

# load packages ----------------------------------------------------
library(tictoc)
library(foreach)
# source functions -------------------------------------------------
source("R/forward_sim_func_wrap.R")
source("R/forward_sim_func_wrap_parallel.R")
source("inst/extdata/scripts/helpers/model_run_helper.R")

# model inputs -----------------------------------------------------
date_of_fit <- "2021-08-25"

output_from_model_fit <- readRDS(paste0("inst/extdata/results/model_fits/output_from_fits_", date_of_fit, ".rds"))
init_cond_22june2021 <- unlist(lapply(unname(output_from_model_fit$`end_date_2021-06-22`), tail,1))
beta_mles <- data.frame(beta = readRDS(paste0("inst/extdata/results/model_fits/mles_from_fits_", date_of_fit,".rds"))) %>%
  mutate(end_date = names(output_from_model_fit))
beta_draws <-  readRDS(paste0("inst/extdata/results/model_fits/beta_draws_from_fits_", date_of_fit, ".rds"))

index <- which(beta_mles$end_date == "end_date_2021-06-22")
cm <- list(baseline_2017 = baseline_2017,
           april_2020 = april_2020,
           june_2020 = june_2020,
           september_2020 = september_2020,
           february_2021 = february_2021,
           june_2021 = june_2021
)

beta_delta_mle <- 0.0003934816 * 2 # R0 = 4.6

# -------------------------------------------------------------------
# run model ---------------------------------------------------------
# -------------------------------------------------------------------
todays_date <- Sys.Date()

# for loop version
tic("for loop")
forward_sim_func_wrap(start_date = "2021-06-22",
                      end_date = "2021-03-31",
                      init_cond = init_cond_22june2021,
                      beta_m = beta_mles[index,1],
                      vac_inputs = basis1,
                      beta_c = beta_delta_mle,
                      beta_draws = beta_draws[[index]][1:3,],
                      contact_matrices = cm,
                      tag = paste0("benchmark_test_for_loop",todays_date))

# parallel version
tic("foreach")
forward_sim_func_wrap_parallel(start_date = "2021-06-22",
                                end_date = "2021-03-31",
                                init_cond = init_cond_22june2021,
                                beta_m = beta_mles[index,1],
                                vac_inputs = basis1,
                                beta_c = beta_delta_mle,
                                beta_draws = beta_draws[[index]][1:3,],
                                contact_matrices = cm,
                                tag = paste0("benchmark_test_foreach",todays_date))
toc()
