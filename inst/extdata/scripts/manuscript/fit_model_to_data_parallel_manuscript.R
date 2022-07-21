# ------------------------------------------------------------------
# Fit model to data script
# ------------------------------------------------------------------
# Getting negative compartment values, so re-writing the entire fit 
# script. The negative compartment values are not reproduced in 
# minimal working example script with the same initial conditions/
# parameter values. Therefore, problem is likely somewhere in the 
# original script.
# ------------------------------------------------------------------

# Load required packages -------------------------------------------
library(deSolve)
library(lubridate)
library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(lubridate)
library(ggplot2)
library(stringr)
library(optimParallel)

source("R/convert_vac_schedule2.R")
source("R/na_to_zero.R")
source("R/calc_waning.R")
source("R/age_struct_seir_ode2.R")
# suppress dplyr::summarise() warnings
options(dplyr.summarise.inform = FALSE)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# Specify initial conditions ----------------------------------------
# define population size (by age group)
n_vec <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463,
           0.14514332, 0.12092904, 0.08807406, 0.04622194) * 17407585 # Dutch population size
empty_state <- c(rep(0, 9)) # vector of zeros
seed_age_group <- sample(1:9,1)
inf_seed_vec <- empty_state
inf_seed_vec[seed_age_group] <- 1

s_vec   <- n_vec - inf_seed_vec; sv1_vec <- empty_state; sv2_vec <- empty_state; sv3_vec <- empty_state; sv4_vec <- empty_state; sv5_vec <- empty_state
shold1_vec <- empty_state; shold2_vec <- empty_state; shold3_vec <- empty_state; shold4_vec <- empty_state; shold5_vec <- empty_state
e_vec <- empty_state; ev1_vec <- empty_state; ev2_vec <- empty_state; ev3_vec <- empty_state; ev4_vec <- empty_state; ev5_vec <- empty_state
i_vec <- inf_seed_vec; iv1_vec <- empty_state; iv2_vec <- empty_state; iv3_vec <- empty_state; iv4_vec <- empty_state; iv5_vec <- empty_state
h_vec <- empty_state; hv1_vec <- empty_state; hv2_vec <- empty_state; hv3_vec <- empty_state; hv4_vec <- empty_state; hv5_vec <- empty_state
ic_vec <- empty_state; icv1_vec <- empty_state; icv2_vec <- empty_state; icv3_vec <- empty_state; icv4_vec <- empty_state; icv5_vec <- empty_state
hic_vec <- empty_state; hicv1_vec <- empty_state; hicv2_vec <- empty_state; hicv3_vec <- empty_state; hicv4_vec <- empty_state; hicv5_vec <- empty_state
d_vec <- empty_state
r_vec <- empty_state; rv1_vec <- empty_state; rv2_vec <- empty_state; rv3_vec <- empty_state; rv4_vec <- empty_state; rv5_vec <- empty_state
r_vec1   <- empty_state; rv1_vec1 <- empty_state; rv2_vec1 <- empty_state; rv3_vec1 <- empty_state; rv4_vec1 <- empty_state; rv5_vec1 <- empty_state
r_vec2   <- empty_state; rv1_vec2 <- empty_state; rv2_vec2 <- empty_state; rv3_vec2 <- empty_state; rv4_vec2 <- empty_state; rv5_vec2 <- empty_state
r_vec3   <- empty_state; rv1_vec3 <- empty_state; rv2_vec3 <- empty_state; rv3_vec3 <- empty_state; rv4_vec3 <- empty_state; rv5_vec3 <- empty_state
# r_vec3 <- n_vec - s_vec - e_vec - i_vec - h_vec - hic_vec - ic_vec - d_vec - r_vec - r_vec1 - r_vec2

init_t0 <- c(t = 0,
               S = s_vec, Sv_1d = sv1_vec, Sv_2d = sv2_vec, Sv_3d = sv3_vec, Sv_4d = sv4_vec, Sv_5d = sv5_vec,
               Shold_1d = shold1_vec, Shold_2d = shold2_vec, Shold_3d = shold3_vec, Shold_4d = shold4_vec, Shold_5d = shold5_vec, 
               E = e_vec, Ev_1d = ev1_vec, Ev_2d = ev2_vec, Ev_3d = ev3_vec, Ev_4d = ev4_vec, Ev_5d = ev5_vec,
               I = i_vec, Iv_1d = iv1_vec, Iv_2d = iv2_vec, Iv_3d = iv3_vec, Iv_4d = iv4_vec, Iv_5d = iv5_vec,
               H = h_vec, Hv_1d = hv1_vec, Hv_2d = hv2_vec, Hv_3d = hv3_vec, Hv_4d = hv4_vec, Hv_5d = hv5_vec,
               IC = ic_vec, ICv_1d = icv1_vec, ICv_2d = icv2_vec, ICv_3d = icv3_vec, ICv_4d = icv4_vec, ICv_5d = icv5_vec,
               H_IC = hic_vec, H_ICv_1d = hicv1_vec, H_ICv_2d = hicv2_vec, H_ICv_3d = hicv3_vec, H_ICv_4d = hicv4_vec, H_ICv_5d = hicv5_vec,
               D = d_vec,
               R = r_vec, Rv_1d = rv1_vec, Rv_2d = rv2_vec, Rv_3d = rv3_vec, Rv_4d = rv4_vec, Rv_5d = rv5_vec,
               R_1w = r_vec1, Rv_1d_1w = rv1_vec1, Rv_2d_1w = rv2_vec1, Rv_3d_1w = rv3_vec1, Rv_4d_1w = rv4_vec1, Rv_5d_1w = rv5_vec1,
               R_2w = r_vec2, Rv_1d_2w = rv1_vec2, Rv_2d_2w = rv2_vec2, Rv_3d_2w = rv3_vec2, Rv_4d_2w = rv4_vec2, Rv_5d_2w = rv5_vec2,
               R_3w = r_vec3, Rv_1d_3w = rv1_vec3, Rv_2d_3w = rv2_vec3, Rv_3d_3w = rv3_vec3, Rv_4d_3w = rv4_vec3, Rv_5d_3w = rv5_vec3
)

# Specify model parameters ------------------------------------------
# define contact/transmission matrix --------------------------------
path <- "inst/extdata/inputs/contact_matrices/converted/"
# path <- "/rivm/s/ainsliek/data/contact_matrices/converted/" 
april_2017     <- readRDS(paste0(path,"transmission_matrix_april_2017.rds"))
april_2020     <- readRDS(paste0(path,"transmission_matrix_april_2020.rds"))
june_2020      <- readRDS(paste0(path,"transmission_matrix_june_2020.rds"))
september_2020 <- readRDS(paste0(path,"transmission_matrix_september_2020.rds"))
february_2021  <- readRDS(paste0(path,"transmission_matrix_february_2021.rds"))
june_2021      <- readRDS(paste0(path,"transmission_matrix_june_2021.rds"))
november_2021  <- readRDS(paste0(path,"transmission_matrix_november_2021.rds"))

# put contact matrices into a list for input into fit_to_data_func()
contact_matrices <- list(
  april_2017 = april_2017,
  april_2020 = april_2020,
  june_2020 = june_2020,
  september_2020 = september_2020,
  february_2021 = february_2021,
  june_2021 = june_2021,
  november_2021 = november_2021
)

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

# vaccination schedule ----------------------------------------------
# read in vaccination schedule
raw_vac_schedule <- read_csv("inst/extdata/inputs/vaccination_schedules/vac_schedule_real_20220709.csv") #%>%
  # select(-X1)
vac_schedule <- raw_vac_schedule[,-1]

# read in xlsx file with VEs (there is 1 sheet for each variant)
# we'll only use wildtype values for now
wt_ve <- read_excel("inst/extdata/inputs/ve_estimates/ve_dat_round2_AB.xlsx", sheet = "wildtype") 
alpha_ve <- read_excel("inst/extdata/inputs/ve_estimates/ve_dat_round2_AB.xlsx", sheet = "alpha")
delta_ve <- read_excel("inst/extdata/inputs/ve_estimates/ve_dat_round2_AB.xlsx", sheet = "delta")
omicron_ve <- read_excel("inst/extdata/inputs/ve_estimates/ve_dat_round2_AB.xlsx", sheet = "omicron")

# -------------------------------------------------------------------
# Define likelihood function ----------------------------------------
# we're fitting the transmission probability (beta) and an
# over-dispersion parameter (alpha)
likelihood_func_test <- function(x, t, data, params, init) {
  library(tidyr)
  # parameters to be estimated
  beta <- x[1]
  alpha <- x[2]
  print(x)
  # observed daily cases
  inc_obs <- data$inc
  
  #print(init)
  # run model with current parameter values
  params$beta <- x[1]/10000
  rk45 <- deSolve::rkMethod("rk45dp7")
  seir_out <- deSolve::ode(init, t, age_struct_seir_ode2, params, method = rk45)  # , rtol = 1e-08, hmax = 0.02
  out <- as.data.frame(seir_out)
  
  print(paste("Negative values?:", any(tail(seir_out, 1) < 0)))
  # modeled cases
  e_comps <- out %>% 
    dplyr::select(starts_with("E"))
  daily_cases <- rowSums(params$sigma * e_comps * params$p_report)
  daily_cases <- ifelse(daily_cases == 0, 0.0001, daily_cases) # prevent likelihood function function from being Inf

  # log-likelihood function
  # lik <- sum(dpois(x = inc_obs,lambda = daily_cases,log=TRUE))
  lik <- -sum(stats::dnbinom(x = inc_obs, mu = daily_cases, size = alpha, log = TRUE))
  
  #print(lik)
  lik
}
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# define fit windows ------------------------------------------------
df_breakpoints <- read_csv2("inst/extdata/inputs/breakpoints_for_model_fit_v3.csv") %>%
  mutate(date = as.Date(date, format = "%d-%m-%Y"),
         time = as.numeric(date - date[1])) %>%
  select(date, time, variant, contact_matrix)

bp_for_fit <- df_breakpoints[1:28,]
n_bp <- dim(bp_for_fit)[1] - 1

# specify initial values and bounds for fitted parameters -----------
fit_params <- list(
  init_value = c(4, 1),
  lower_bound = c(0.5, 0.1),
  upper_bound = c(Inf, Inf)
)

# create empty containers to store outputs from fit -----------------
init_cond <- list()     # store initial conditions for each window
out <- list()           # model output for each time point
times <- list()         # time points
cases <- list()         # daily cases to plot against real data
mles <- list()          # MLEs for each time window
beta_draws <- list()    # store 200 parameter draws
# ci_out <- list()        # store model outputs for confidence bounds
# ci_cases <- list()

# load case data ----------------------------------------------------
data_date <- "2022-07-10"
case_data <- readRDS(paste0("inst/extdata/data/case_data_upto_", data_date, ".rds"))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Fit model to data -------------------------------------------------
init_cond[[1]] <- init_t0

# read in initial conditions if not starting at iteration 1
# init_cond <- readRDS("inst/extdata/results/model_fits/initial_conditions.rds")

# loop over time windows --------------------------------------------
for (j in 1:n_bp) {
  print(paste(paste0(j, ")"), "Fitting from", bp_for_fit$date[j], "to", bp_for_fit$date[j + 1]))

  # set contact matrix for time window ------------------------------
  if (bp_for_fit$contact_matrix[j + 1] == "april_2017") {
    contact_matrix <- contact_matrices$april_2017 # ; print("april_2017")
  } else if (bp_for_fit$contact_matrix[j + 1] == "april_2020") {
    contact_matrix <- contact_matrices$april_2020 # ; print("april_2020")
  } else if (bp_for_fit$contact_matrix[j + 1] == "june_2020") {
    contact_matrix <- contact_matrices$june_2020 # ; print("june_2020")
  } else if (bp_for_fit$contact_matrix[j + 1] == "september_2020") {
    contact_matrix <- contact_matrices$september_2020 # ; print("septemeber_2020")
  } else if (bp_for_fit$contact_matrix[j + 1] == "february_2021") {
    contact_matrix <- contact_matrices$february_2021 # ; print("february_2021")
  } else if (bp_for_fit$contact_matrix[j + 1] == "june_2021") {
    contact_matrix <- contact_matrices$june_2021 # ; print("june_2021")
  } else {
    contact_matrix <- contact_matrices$november_2021
  }

  # has vaccination started? ----------------------------------------
  # nv <- ifelse(bp_for_fit$date[j+1] >= as.Date("2021-01-04"), TRUE, FALSE)

  # convert vaccination schedule for input into model ---------------
  if (bp_for_fit$variant[j + 1] == "wildtype") {
    ve_params <- wt_ve
  } else if (bp_for_fit$variant[j + 1] == "alpha") {
    ve_params <- alpha_ve
  } else if (bp_for_fit$variant[j + 1] == "delta") {
    ve_params <- delta_ve
  } else if (bp_for_fit$variant[j + 1] == "omicron") {
    ve_params <- omicron_ve
  }

  vac_rates <- convert_vac_schedule2(
    vac_schedule = vac_schedule,
    ve_pars = ve_params,
    wane = TRUE,
    k_inf = 0.006,
    k_sev = 0.012
  )

  # data wrangle for model input
  df_input <- pivot_wider(vac_rates %>%
    filter(param != "comp_ve") %>%
    mutate(param = ifelse(param == "comp_delay", "delay", param)),
  names_from = c("param", "age_group"),
  names_sep = "", values_from = "value"
  )

  # parameters must be in a named list
  params <- list(
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
    omega = 0.02017514,
    # daily vaccination rate
    alpha1 = df_input %>%
      filter(dose == "d1", outcome == "infection") %>%
      select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
    alpha2 = df_input %>%
      filter(dose == "d2", outcome == "infection") %>%
      select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
    alpha3 = df_input %>%
      filter(dose == "d3", outcome == "infection") %>%
      select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
    alpha4 = df_input %>%
      filter(dose == "d4", outcome == "infection") %>%
      select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
    alpha5 = df_input %>%
      filter(dose == "d5", outcome == "infection") %>%
      select(date, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9),
    # delay to protection
    delay1 = df_input %>%
      filter(dose == "d1", outcome == "infection") %>%
      select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
    delay2 = df_input %>%
      filter(dose == "d2", outcome == "infection") %>%
      select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
    delay3 = df_input %>%
      filter(dose == "d3", outcome == "infection") %>%
      select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
    delay4 = df_input %>%
      filter(dose == "d4", outcome == "infection") %>%
      select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
    delay5 = df_input %>%
      filter(dose == "d5", outcome == "infection") %>%
      select(date, delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8, delay9),
    # protection against infection
    eta1 = df_input %>%
      filter(dose == "d1", outcome == "infection") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta2 = df_input %>%
      filter(dose == "d2", outcome == "infection") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta3 = df_input %>%
      filter(dose == "d3", outcome == "infection") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta4 = df_input %>%
      filter(dose == "d4", outcome == "infection") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta5 = df_input %>%
      filter(dose == "d5", outcome == "infection") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    # protection from hospitalisation
    eta_hosp1 = df_input %>%
      filter(dose == "d1", outcome == "hospitalisation") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta_hosp2 = df_input %>%
      filter(dose == "d2", outcome == "hospitalisation") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta_hosp3 = df_input %>%
      filter(dose == "d3", outcome == "hospitalisation") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta_hosp4 = df_input %>%
      filter(dose == "d4", outcome == "hospitalisation") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta_hosp5 = df_input %>%
      filter(dose == "d5", outcome == "hospitalisation") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    # protection from transmission
    eta_trans1 = df_input %>%
      filter(dose == "d1", outcome == "transmission") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta_trans2 = df_input %>%
      filter(dose == "d2", outcome == "transmission") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta_trans3 = df_input %>%
      filter(dose == "d3", outcome == "transmission") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta_trans4 = df_input %>%
      filter(dose == "d4", outcome == "transmission") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    eta_trans5 = df_input %>%
      filter(dose == "d5", outcome == "transmission") %>%
      select(date, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9),
    p_report = p_reported_by_age,
    contact_mat = contact_matrix, # $mean,  # change contact matrix
    calendar_start_date = as.Date("2020-01-01") # ,
    # no_vac = nv
  )

  # set time sequence -----------------------------------------------
  times[[j]] <- seq(bp_for_fit$time[j], bp_for_fit$time[j + 1], by = 1)

  # subset data for time window -------------------------------------
  case_data_sub <- case_data[times[[j]] + 1, ]

  # Start cluster
  ## - use all avilable processor cores
  ## - return cat() output to R prompt
  ## (may have issues on Windows)
  if (tolower(.Platform$OS.type) != "windows") {
    cl <- makeCluster(spec = detectCores(), type = "FORK", outfile = "")
  } else {
    cl <- makeCluster(spec = detectCores(), outfile = "")
  }
  setDefaultCluster(cl = cl)
  clusterExport(cl = cl, c("age_struct_seir_ode2"))
  ## return log information
  options(optimParallel.loginfo = TRUE)
  ## stop if change of f(x) is smaller than 0.01
  control <- list(factr = .001 / .Machine$double.eps)

  # run optimization procedure --------------------------------------
  res <- optimParallel(
    par = fit_params$init_value,
    fn = likelihood_func_test,
    method = "L-BFGS-B",
    lower = fit_params$lower_bound,
    upper = fit_params$upper_bound,
    t = times[[j]],
    data = case_data_sub,
    params = params,
    init = unlist(init_cond[[j]]),
    hessian = TRUE,
    control = control
  )

  setDefaultCluster(cl = NULL)
  stopCluster(cl)
  # store MLE --------------------------------------------------------
  mles[[j]] <- c(beta = res$par[1] / 10000, alpha = res$par[2])
  print(mles[[j]])
  saveRDS(mles, "inst/extdata/results/model_fits/manuscript/mle_list.rds")
  # Run model --------------------------------------------------------
  params$beta <- res$par[1] / 10000
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(unlist(init_cond[[j]]), times[[j]], age_struct_seir_ode2,
    params,
    method = rk45
  ) # , rtol = 1e-08, hmax = 0.02

  # checks -----------------------------------------------------------
  # output error message if negative compartment values
  if (any(tail(seir_out, 1) < 0)) {
    stop("Negative compartment values")
  }

  # check population size
  if (!all.equal(sum(unlist(tail(seir_out, 1)[-c(1:2)])), sum(params$N))) {
    stop("Number of individuals in compartments does not sum to population size")
  }
  # store outputs ----------------------------------------------------
  out[[j]] <- as.data.frame(seir_out)
  e_comps <- out[[j]] %>% dplyr::select(starts_with("E"))
  cases[[j]] <- rowSums(params$sigma * e_comps * params$p_report)
  # saveRDS(cases, "inst/extdata/results/model_fits/modelled_daily_cases.rds")
  # plot for quick check of fit --------------------------------------
  plot(case_data_sub$inc ~ times[[j]],
    pch = 16, col = "red",
    ylim = c(0, max(case_data_sub$inc, cases[[j]]))
  )
  lines(cases[[j]] ~ times[[j]])

  #------------------------------------------------------------------
  # get confidence bounds -------------------------------------------
  # draw 200 parameter values
  parameter_draws <- mvtnorm::rmvnorm(200, res$par, solve(res$hessian))
  beta_draws[[j]] <- data.frame(beta = (parameter_draws[, 1] / 10000)) %>%
    mutate(index = 1:200)
  saveRDS(beta_draws, "inst/extdata/results/model_fits/manuscript/beta_draws.rds")
  # run model for each beta draw (with different contact matrix) ----
  # ci_out[[j]] <- list()
  # ci_cases[[j]] <- list()
  # for(i in 1:200){
  #   params$beta <- beta_draws[[j]][i,1]
  #   params$contact_mat <- contact_matrix[[i]]
  #   seir_out_ci <- ode(init_cond[[j]], times[[j]], age_struct_seir_ode2,
  #                      params, method = rk45, rtol = 1e-08, hmax = 0.02)
  #   seir_out_ci1 <- as.data.frame(seir_out_ci)
  #   ci_cases[[j]][[i]] <-  rowSums(params$sigma * seir_out_ci1[c(paste0("E",1:9))] * params$p_report)
  # }
  # ci_out[[j]] <- do.call("rbind", ci_cases[[j]])
  # -----------------------------------------------------------------

  # update initial conditions for next time window
  init_cond[[j + 1]] <- tail(out[[j]], 1)[-1]
  saveRDS(init_cond, "inst/extdata/results/model_fits/manuscript/initial_conditions_manuscript.rds")
  # ------------------------------------------------------------------
} # end of for loop over breakpoints


# -------------------------------------------------------------------
# -------------------------------------------------------------------

# Plot output -------------------------------------------------------
# plot all cases with confidence bounds
# ci_out_wide <- do.call("cbind", ci_out)
# bounds <- apply(ci_out_wide, 2, quantile, probs = c(0.025, 0.975)) # get quantiles

df_model_fit <- data.frame(time = x_axis, 
                           date = params$calendar_start_date + x_axis,
                           obs = case_data$inc[x_axis + 1], 
                           mle = unlist(cases)#, 
                           # lower = bounds[1,], 
                           # upper = bounds[2,]
                           )

path_out <- "inst/extdata/results/model_fits/manuscript/"
saveRDS(df_model_fit,
        file = paste0(path_out, "model_fit_df_from_", df_model_fit$date[1],"_to_",
                      tail(df_model_fit$date,1), ".rds"))

p <- ggplot(data = df_model_fit, aes(x = date, y = mle, linetype="solid")) +
  geom_point(data = df_model_fit, aes(x = date, y = obs, color = "Osiris notifications")) +
  geom_line() +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Confidence bounds"), alpha = 0.3) +
  scale_color_manual(values = c("red"),
                     labels = c("Osiris notifications")) +
  scale_fill_manual(values = c("grey70")) +
  scale_linetype_manual(values=c(1), labels = c("Model Fit")) +
  #scale_shape_manual(values=c(NA,20)) +
  labs(y = "Daily Cases", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title=element_text(size=14))
p
# --------------------------------------------------------------------

# start_time <- Sys.time()
# seir_out <- ode(unlist(init_cond[[j]]), times[[j]], age_struct_seir_ode2,  
#                 params, method = rk45) # , rtol = 1e-08, hmax = 0.02
# end_time <- Sys.time()
# end_time - start_time
