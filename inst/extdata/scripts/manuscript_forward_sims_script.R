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

# specify time points ----------------------------------------------
start_date <- lubridate::yday(as.Date("2021-06-22") + 365) + 365
end_date <- lubridate::yday(as.Date("2022-03-30")) + (365*2) 
times <- seq(start_date, end_date, by = 1)

initial_conditions <- c(t=times[1], init_cond_22june2021)

# Create list of parameter values for input into model solver ------
params <- list(beta = beta_mles[index,1], # transmission rate
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
               c_start = june_2021$mean,
               c_lockdown = february_2021$mean,
               c_relaxed = june_2020$mean,
               c_very_relaxed = june_2021$mean,
               c_normal = baseline_2017$mean,
               keep_cm_fixed = FALSE,
               vac_inputs = basis1,
               use_cases = TRUE,                           # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),         # somewhat arbitrary cut-off ***need to check if realistic
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 100000/100000 * sum(n_vec),        #35.7  # 20 for IC admissions
               no_vac = FALSE,
               t_calendar_start = yday(as.Date("2020-01-01")),   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
               beta_change =  0.0003934816 * 2  
                                                            # mle = 0.0003934816 * 2 (R0 = 4.6); 
                                                            # lower = 0.0005902224 (R0 = 3.45); 
                                                            # upper = 0.001539711 (R0 = 9)
)


# if time doesn't start at 0 we need to initialise the contact 
# matrices flags ---------------------------------------------------
flag_relaxed <- 0 
flag_very_relaxed <- 0
flag_normal <- 0

# Solve model ------------------------------------------------------
# mle
seir_out <- lsoda(initial_conditions,times,age_struct_seir_ode,params) 
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)

# get outcomes ------------------------------------------------------
cases_mle <- params$sigma * rowSums(out$E + out$Ev_1d + out$Ev_2d) * params$p_report
infectious_mle <- (out$I + out$Iv_1d + out$Iv_2d)
hosp_admissions_mle <- rowSums(sweep(infectious_mle, 2, h, "*"))
hosp_occ_mle <- (out$H + out$Hv_1d + out$Hv_2d)
ic_mle <- rowSums(sweep(hosp_occ_mle, 2, i1, "*"))
ic_occ_mle <- (out$IC + out$ICv_1d + out$ICv_2d)
hosp_after_ic_mle <- sweep(ic_occ_mle, 2, i2, "*")
hosp_after_ic_occ_mle <- (out$H_IC + out$H_ICv_1d + out$H_ICv_2d)
deaths_mle <- rowSums(sweep(ic_occ_mle, 2, d_ic, "*") + sweep(hosp_occ_mle, 2, d, "*") + sweep(hosp_after_ic_occ_mle, 2, d_hic, "*"))

# calculate beta ----------------------------------------------------
# init_cond_05june2021 <- c(t = 521, unlist(lapply(unname(output_from_model_fit$`end_date_2021-06-05`), tail,1)))
# S_diag <- diag(init_cond_05june2021[c(2:10)])
# rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
# beta_draws <- (parameter_draws[[index]][,1] / rho) * params$gamma

# run for beta draws ------------------------------------------------
rtn_cases <- matrix(,nrow = length(times), ncol = length(beta_draws))
rtn_hosp <- rtn_cases
rtn_ic <- rtn_cases
rtn_deaths <- rtn_cases

for(i in 1:length(beta_draws)){
  print(i)
  
  flag_relaxed <- 0
  flag_very_relaxed <- 0
  flag_normal <- 0
  
  params$beta <- beta_draws[i]
  params$c_start <- june_2021[[i]]
  params$normal <- baseline_2017[[i]]
  
  seir_out <- lsoda(initial_conditions, times, age_struct_seir_ode, params) #
  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output(seir_out)
  # get outcomes
  daily_cases <- params$sigma * rowSums(out$E + out$Ev_1d + out$Ev_2d) * params$p_report
  infectious <- (out$I + out$Iv_1d + out$Iv_2d)
  hosp_admissions <- rowSums(sweep(infectious, 2, h, "*"))
  hosp_occ <- (out$H + out$Hv_1d + out$Hv_2d)
  ic <- rowSums(sweep(hosp_occ, 2, i1, "*"))
  ic_occ <- (out$IC + out$ICv_1d + out$ICv_2d)
  hosp_after_ic <- sweep(ic_occ, 2, i2, "*")
  hosp_after_ic_occ <- (out$H_IC + out$H_ICv_1d + out$H_ICv_2d)
  daily_deaths <- rowSums(sweep(ic_occ, 2, d_ic, "*") + sweep(hosp_occ, 2, d, "*") + sweep(hosp_after_ic_occ, 2, d_hic, "*"))
  
  rtn_cases[,i] <- daily_cases
  rtn_hosp[,i]   <- hosp_admissions
  rtn_ic[,i]     <- ic
  rtn_deaths[,i] <- daily_deaths
}

# run model for each combination of parameters - DOES NOT WORK
# out_mult <- apply(beta_draws[1:5,], 1, function_wrapper_2,
#              contact_matrix = june_2021,
#              init = initial_conditions, t = times) # rows are time points, columns are different simulations

# get confidence bounds of model runs
bounds_cases  <- apply(rtn_cases,1,function(x) quantile(x, c(0.025,0.975)))
bounds_hosp   <- apply(rtn_hosp,1,function(x) quantile(x, c(0.025,0.975)))
bounds_ic     <- apply(rtn_ic,1,function(x) quantile(x, c(0.025,0.975)))
bounds_deaths <- apply(rtn_deaths,1,function(x) quantile(x, c(0.025,0.975)))

# quick check ------------------------------------------------------
# cases_from_fit <- (params$sigma * (init_from_fit1$E + init_from_fit1$Ev_1d + init_from_fit1$Ev_2d)) / 3
df_check <- data.frame(time = times,
                       cases_mle = cases_mle,
                       cases_lower = bounds_cases[1,],
                       cases_upper = bounds_cases[2,],
                       hosp_mle = hosp_admissions_mle,
                       hosp_lower = bounds_hosp[1,],
                       hosp_upper = bounds_hosp[2,],
                       ic_mle = ic_mle,
                       ic_lower = bounds_ic[1,],
                       ic_upper = bounds_ic[2,],
                       deaths_mle = deaths_mle,
                       deaths_lower = bounds_deaths[1,],
                       deaths_upper = bounds_deaths[2,]
) 

df_plot <- df_check %>%
  pivot_longer(!time, names_to = c("outcome", "estimate"), names_sep = "_", values_to = "value") %>%
  pivot_wider( names_from = "estimate", values_from = "value") %>%
  mutate(outcome = factor(outcome, levels = c("cases", "hosp", "ic", "deaths")))

p <- ggplot(data = df_plot, aes(x = time, y = mle)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.3) +
  theme_bw() +
  facet_wrap(~outcome, scales = "free")
p

# Summarise results ------------------------------------------------
tag <- "df_basis_18plus_mle_beta_23aug"
# results <- summarise_results(out, params, start_date = "2021-01-31", 
#                              times = times, vac_inputs = params$vac_inputs)
saveRDS(df_plot, paste0("inst/extdata/results/",tag,".rds"))
