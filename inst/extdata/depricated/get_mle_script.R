# --------------------------------------------------
# find mle of beta from SEIR model fit to case data
# --------------------------------------------------

# load packages
devtools::load_all()
library(vacamole)

source("inst/extdata/scripts/helpers/model_run_helper.R")
source("inst/extdata/scripts/helpers/read_case_data_from_server.R")
# source("R/model_run_wrapper.R")
# source("R/likelihood_func.R")

# Read in OSIRIS data ------------------------------
last_date_in_osiris <- "2021-11-28"
case_data <- readRDS(paste0("inst/extdata/data/case_data_upto_", last_date_in_osiris,".rds"))
 
# plot data
p <- ggplot(case_data, aes(x = date, y = inc)) +
  geom_line() +
  geom_line(aes(x = date, y = roll_avg, color = "red")) +
  theme(panel.background = element_blank())
p

# --------------------------------------------------
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

# read in vac schedule
vac_schedule <- read_csv("/rivm/EPI/MOD/Projects/NovelCoronaWuhan/vaccineallocation/Cum_upt20211201.csv") %>%
  select(-starts_with("X"))

# convert vac schedule
basis_alpha <- convert_vac_schedule(
  vac_schedule = vac_schedule,
  ve = ve,
  hosp_multiplier = h_multiplier,
  delay = delays,
  ve_trans = ve_trans,
  wane = FALSE,
  add_child_vac = FALSE,
  add_extra_dates = TRUE,
  extra_start_date = "2020-01-01",
  extra_end_date = "2021-01-03"
)

basis_delta <- convert_vac_schedule(
  vac_schedule = vac_schedule,
  ve = ve_delta,
  hosp_multiplier = h_multiplier_delta,
  delay = delays,
  ve_trans = ve_trans_delta,
  wane = FALSE,
  add_child_vac = FALSE,
  add_extra_dates = TRUE,
  extra_start_date = "2020-01-01",
  extra_end_date = "2021-01-03"
)


# specify model parameters
params <- list(beta = 0.0003934816,       # transmission rate
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
               p_report = 1/3,
               c_start = april_2017$mean,
               c_lockdown = april_2020$mean,
               c_relaxed = september_2020$mean,
               c_very_relaxed = june_2020$mean,
               c_normal = april_2017$mean,
               keep_cm_fixed = TRUE,
               vac_inputs = basis_alpha,
               use_cases = TRUE,                           # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 35.7/100000 * sum(n_vec),        # 35.7  # 20 for IC admissions
               no_vac = FALSE,
               t_calendar_start = yday(as.Date("2020-01-01")),   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
               beta_change = NULL
               )


# --------------------------------------------------
# fit model to data for each time window
# --------------------------------------------------

breakpoints <- list( 
  date = c( 
    as.Date("2020-03-16"),  # 1) closure of hospitality, schools, daycares
    as.Date("2020-03-24"),  # 2) casinos subject to same measures as food/beverage outlets
    as.Date("2020-04-29"),  # 3) Children and young people can participate in outdoor activities
    as.Date("2020-05-11"),  # 4) reopening of primary schools and daycare @ 50% capacity
    as.Date("2020-06-01"),  # 5) cafes/restaurants reopen with max 30 visitor,
                            # primary schools and daycares open at 100% 
    as.Date("2020-07-01"),  # 6) secondary schools reopen @ 100% capacity
    as.Date("2020-08-06"),  # 7) Small groups allowed at universities
    as.Date("2020-08-18"),  # 8) maximum 6 people in household
    as.Date("2020-09-20"),  # 9) groups <= 50 people, some hospitality closes at midnight or 1 AM 
    as.Date("2020-09-29"),  # 10) max 3 visitors at home, max group size 30, face masks in public areas
    as.Date("2020-10-14"),  # 11) cafes/restaurants close, no team sports, 3 visitors at home per day
    as.Date("2020-10-23"),  # 12) hotels can't sell alcohol after 20:00
    as.Date("2020-11-04"),  # 13) 2 visitors per day, groups <= 20
    as.Date("2020-11-19"),  # 14) 3 visitors per day, groups <= 30
    as.Date("2020-12-01"),  # 15) masks mandatory in all public and indoor areas
    as.Date("2020-12-15"),  # 16) non-essential shops close
    as.Date("2021-01-01"),  # 17) end of year
    as.Date("2021-01-20"),  # 18) 1 visitor per day/curfew (23/01/2021)
    as.Date("2021-02-08"),  # 19) Primary schools, child care, special ed reopen
    as.Date("2021-03-01"),  # 20) secondary schools partially reopen, contact professions reopen (3/3/2021)
    as.Date("2021-03-16"),  # 21) some retail reopens
    as.Date("2021-03-31"),  # 22) curfew to start at 22:00 instead of 21:00
    as.Date("2021-04-19"),  # 23) out of school care fully reopens
    as.Date("2021-04-28"),  # 24) curfew canceled
    as.Date("2021-05-19"),  # 25) <27 can play outdoor sports, groups <= 30, non-essential travel allowed within NL 
    as.Date("2021-06-05"),  # 26) 4 visitors per day, museums reopen, group <= 50, restaurants reopen
    as.Date("2021-06-22"),  # 27) vaccination of 12-17 year olds starts, Delta becomes dominant strain
    as.Date("2021-06-26"),  # 28) all restrictions relaxed, except masks on public transport, nightclubs reopen
    as.Date("2021-07-10"),  # 29) catering industry reopens, test for entry with large events, nightclubs close
    as.Date("2021-07-19"),  # 30) work from home advisory re-instated
    as.Date("2021-08-30"),  # 31) Physical classes resume at MBO schools, colleges and universities. 
                            #     The 1.5 meters no longer applies. Maximum group size = 75, face masks 
                            #     must be worn outside the classrooms or lecture halls.
    as.Date("2021-09-25"),  # 32) 1.5 m distance no longer required, coronapass needed in horeca and events
    as.Date("2021-11-03"),  # 33) coronapass use expanded, work from home half the time, masks in public spaces
    as.Date("2021-11-13"),  # 34) non-essential shops close at 5pm, essential shops close at 8pm
    as.Date("2021-11-28"),  # 35) stay at home advise, 1.5 m distance, work from home, self-tests advised before visiting other, masks in schools
    as.Date("2021-12-19"),  # 36) 2 people in house per day, non-essential shops closed, non-medical contact professions closed,
                            #     schools closed, all events banned
    as.Date("2022-01-15"),  # 37) schools re-open, non-essential shops re-open until 5PM, face masks advised everywhere
    as.Date("2022-01-26"),  # 38) shops can stay open until 10PM, max number visitors increased to 4 per day, catering,/theaters/museums re-open, quesrantine not required for people with booster
    as.Date("2022-02-15"),  # 39) people can work from the office half the time, unlimited visitors at home
    as.Date("2022-02-25"),  # 40) most measures lifted, except corona pass, masks on public transport, testing for large events
    as.Date("2022-03-12")   # last date of fit.
  ),  
  contact_matrix = list( april_2017, 
                         april_2017, 
                         april_2020,
                         april_2020, 
                         april_2020, 
                         june_2020, 
                         june_2020,
                         june_2020,
                         september_2020,
                         september_2020,
                         september_2020, 
                         september_2020,
                         september_2020,
                         september_2020,
                         september_2020,
                         september_2020,
                         september_2020,
                         september_2020,
                         february_2021,
                         february_2021,
                         february_2021,
                         february_2021,
                         february_2021,
                         february_2021,
                         february_2021,
                         june_2021,
                         june_2021,
                         june_2021,
                         june_2021,
                         june_2021,
                         june_2021,
                         june_2021,
                         june_2021,
                         june_2021,
                         june_2021
                         )
)

# TODO wrap this into function!!!
n_bp <- length(breakpoints$date)
mles <- matrix(rep(NA, 2*n_bp), nrow = n_bp)
colnames(mles) <- c("beta", "alpha")
out_mle <- list()
parameter_draws <- list()
beta_draws <- list()
daily_cases_mle <- list()

for (j in 1:n_bp) {
  print(j)
  # set contact matrix for time window
  if (j == n_bp){
    params$c_start <- breakpoints$contact_matrix[[j-1]]$mean
  } else {
    params$c_start <- breakpoints$contact_matrix[[j]]$mean
  }
  
  # set VE for time window
  if (breakpoints$date[j] > as.Date("2021-06-21")){
    params$vac_inputs <- basis_delta
  } else {params$vac_inputs <- basis_alpha}
  
  if (j == 1) {
    # if first time window, start time at 0
    end_day <- yday(breakpoints$date[j]) - 1
    
    times <- seq(0, end_day, by = 1)
    # set initial conditions to those specified earlier in script
    init_update <- init
    pars <- c(2.3, 0.01)
    S_diag <- diag(init_update[c(2:10)])
    rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  } else {
    if (breakpoints$date[j] > as.Date("2020-12-31")){
      if(breakpoints$date[j-1] > as.Date("2020-12-31")){ # wait for two consecutive dates in 2021
        start_day <- yday(breakpoints$date[j-1]) - 1 + 366 # shift days by 1 because we start time at 0 (not 1)
      } else {
        start_day <- yday(breakpoints$date[j-1]) - 1
      }
      end_day <- yday(breakpoints$date[j]) - 1 + 366
    } else {
      start_day <- yday(breakpoints$date[j-1]) - 1 
      end_day <- yday(breakpoints$date[j]) - 1
    }
    times <- seq(start_day, end_day, by = 1)
    # update initial conditions based on last time window
    init_update <- c(t = times[1], unlist(lapply(unname(out_mle[[j-1]]), tail,1)))
    
    pars <- c((mles[j-1,1]/params$gamma)*rho, mles[j-1,2])
    S_diag <- diag(init_update[c(2:10)])
    rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  }
  
  # subset data for time window
  case_data_sub <- case_data[times + 1, ]
  
  # optimize
  res <- optim(par = pars, 
               fn = likelihood_func,
               method = "L-BFGS-B",
               lower = c(0,0.005),
               upper = c(10,1),
               t = times,
               data = case_data_sub,
               params = params,
               init = init_update,
               stochastic = FALSE,
               hessian = TRUE
  )
  
  # store MLE
  mles[j,1] <- (res$par[1] / rho) * params$gamma
  mles[j,2] <- res$par[2]
  
  # draw 200 parameter values
  parameter_draws[[j]] <- mvtnorm::rmvnorm(200, res$par, solve(res$hessian))
  print(solve(res$hessian))
  beta_draws[[j]] <- data.frame(beta = (parameter_draws[[j]][,1] / rho) * params$gamma) %>%
     mutate(index = 1:200)
  # --------------------------------------------------
  # run for mle to get initial conditions for next timepoint
  params$beta <- mles[j,1]
  seir_out <- lsoda(init_update, times, age_struct_seir_ode, params)
  seir_out <- as.data.frame(seir_out)
  out_mle[[j]] <- postprocess_age_struct_model_output(seir_out)
  cases <- params$sigma * rowSums(out_mle[[j]]$E + out_mle[[j]]$Ev_1d + out_mle[[j]]$Ev_2d) * params$p_report
  
  # plot for quick check of fit
  plot(cases~times, type = "l")
  points(times, case_data_sub$inc, pch = 16, col = "red")

} # end of for loop over breakpoints

todays_date <- Sys.Date()
path_out <- "/rivm/s/ainsliek/code/vacamole/inst/extdata/results/model_fits/"
# save outputs
saveRDS(mles, file = paste0(path_out, "mles_from_fits_", todays_date, ".rds"))
saveRDS(beta_draws, file = paste0(path_out, "beta_draws_from_fits_", todays_date, ".rds"))
names(out_mle) <- paste0("end_date_", breakpoints$date) # name list elements for easier indexing
saveRDS(out_mle, file = paste0(path_out, "output_from_fits_", todays_date, ".rds"))

# ----------------------------------------------------
# run simulations for mle, lower, and upper bounds 
# of beta
# ----------------------------------------------------
fit_date <- "2021-10-01"
beta_mles <- readRDS(paste0(path_out,"mles_from_fits_",fit_date,".rds"))
beta_mles_list <- split(beta_mles, seq(nrow(beta_mles)))
beta_draws <- readRDS(paste0(path_out,"beta_draws_from_fits_",fit_date,".rds"))

# run for 200 contact matrices
mle_run <- model_run_wrapper(breakpoints = breakpoints, beta_values = beta_mles_list, init_conditions = init, params = params)
ci_run  <- model_run_wrapper(breakpoints = breakpoints, beta_values = beta_draws, init_conditions = init, params = params, mle = FALSE)
ci_out <- list()
for (i in 1:n_bp){
  ci_out[[i]] <- do.call("rbind", ci_run[[i]])
}
ci_out_wide <- do.call("cbind", ci_out)
matplot(t(ci_out_wide), type = "l")

# try getting quantiles
bounds <- apply(ci_out_wide, 2, quantile, probs = c(0.025, 0.975))
matplot(t(bounds), type = "l")

# save outputs -------------------------------------
# --------------------------------------------------
#  combine all piecewise results to plot together
cases_mle <- unique(unlist(mle_run))
cases_lower <- unique(bounds[1,])
cases_upper <- unique(bounds[2,])
times_all <- 1:length(cases_mle)

model_fit <- data.frame(time = times_all, date = case_data$date, real = case_data$inc, mle = cases_mle, lower = cases_lower, upper = cases_upper)
saveRDS(model_fit, file = paste0(path_out, "model_fit_df_", todays_date, ".rds"))
# --------------------------------------------------
