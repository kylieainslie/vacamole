# find mle of beta from SEIR model fit to real data

# --------------------------------------------------
# load packages
#library(segmented)
library(ggplot2)
library(tidyverse)
library(lubridate)

# Read in OSIRIS data ------------------------------
source("inst/extdata/scripts/model_run_helper.R")
# read in OSIRIS data
osiris <- readRDS("inst/extdata/data/real_data/Osiris_Data_20210715_1054.rds")

osiris1 <- osiris %>%
  filter(!is.na(date)) %>%
  complete(date = seq.Date(min(date), max(date), by="day"), fill = list(inc = 0)) %>%
  mutate(roll_avg = zoo::rollmean(inc, k = 7, fill = 0))

osiris2 <- osiris1 %>%
  filter(date <= as.Date("2020-03-16")) # date of first lockdown

time_vec <- seq(0, nrow(osiris2)-1, by = 1)

# plot data
p <- ggplot(osiris1, aes(x = date, y = inc)) +
  geom_line() +
  geom_line(aes(x = date, y = roll_avg, color = "red")) +
  theme(panel.background = element_blank())
p

# Read in contact matrices -------------------------
baseline_2017 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_baseline_2017.rds")
april_2020 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_april_2020.rds")
june_2020 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_june_2020.rds")
september_2020 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_september_2020.rds")
february_2021 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_february_2021.rds")
june_2021 <- readRDS("inst/extdata/data/contact_matrices/contact_matrices_june_2021.rds")

# --------------------------------------------------
# specify model parameters
params <- list(dt = 1,
               beta = 0.0004,             # transmission rate
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
               c_start = baseline_2017$mean,
               c_lockdown = april_2020$mean,
               c_relaxed = september_2020$mean,
               c_very_relaxed = june_2020$mean,
               c_normal = baseline_2017$mean,
               keep_cm_fixed = TRUE,
               vac_inputs = NULL,
               use_cases = TRUE,                           # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 35.7/100000 * sum(n_vec),        # 35.7  # 20 for IC admissions
               no_vac = FALSE,
               t_calendar_start = yday(as.Date("2020-01-01"))   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
)

# initial values
empty_state <- c(rep(0, 9))
init <- c(
  t = 0,
  S = n_vec - 1,
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

# --------------------------------------------------
# likelihood function
likelihood_func <- function(x,
                            contact_matrix,
                            t,
                            data,
                            params,
                            init,
                            stochastic = FALSE) {
  #params$beta <- x[1] # pars["beta"]
  
  # params$beta1 <- beta1 #pars["beta1"]
  r0 <- x[1]
  S_diag <- diag(init[c(2:10)])
  rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  params$beta <- (r0 / rho) * params$gamma
  
  if (stochastic){
    seir_out <- stochastic_age_struct_seir_ode(times = t,init = init, params = params)
    out <- apply(seir_out, 3, rowSums)
    daily_cases <- params$sigma * (out[,"E"] + out[,"Ev_1d"] + out[,"Ev_2d"]) * params$p_report
    daily_cases <- ifelse(daily_cases == 0, 0.0001, daily_cases) # prevent likelihood function function from being Inf
  } else {
    seir_out <- lsoda(init,t,age_struct_seir_ode,params) #
    seir_out <- as.data.frame(seir_out)
    out <- postprocess_age_struct_model_output(seir_out)
    daily_cases <- (params$sigma * rowSums(out$E + out$Ev_1d + out$Ev_2d)) * params$p_report
    daily_cases <- ifelse(daily_cases == 0, 0.0001, daily_cases) # prevent likelihood function function from being Inf
  }
  
  inc_obs <- data$inc

  # lik <- sum(dpois(x = inc_obs,lambda = incidence,log=TRUE))
  alpha <- x[2]
  size <- daily_cases * (alpha/(1-alpha))
  lik <- -sum(dnbinom(x = inc_obs, mu = daily_cases, size = size, log = TRUE))
  #print(lik)
  lik
}
# ---------------------------------------------------
# Determine MLE using optim

res <- optim(par = c(2.3, 0.01), 
             fn = likelihood_func,
             method = "L-BFGS-B",
             lower = c(1,0.0001),
             upper = c(10,1),
             t = time_vec,
             data = osiris2,
             params = params,
             init = init,
             stochastic = FALSE,
             hessian = TRUE
             )

res$par

# plot to check fit --------------------------------
S_diag <- diag(init[c(2:10)])
rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
params$beta <- (res$par[1] / rho) * params$gamma
#params$beta <- res$par[1]
seir_out <- lsoda(init, time_vec, age_struct_seir_ode, params) #
seir_out <- as.data.frame(seir_out)
out_mle <- postprocess_age_struct_model_output(seir_out)
daily_cases_mle <- params$sigma * rowSums(out_mle$E + out_mle$Ev_1d + out_mle$Ev_2d) * params$p_report

plot(daily_cases_mle ~ time_vec, type = "l") # ylim = c(0, 8000)
points(osiris2$inc ~ time_vec, col = "red", pch = 16)

# Hessian Magic ------------------------------------
# estres_exp <- optim(previous_estimates$par, minloglik_exp,
#                     ts =  tvector, Ns = Nvector, obs = observeds, method = "BFGS", hessian = TRUE)
parameter_draws <- mvtnorm::rmvnorm(200, res$par, solve(res$hessian))
betas <- data.frame(beta = (parameter_draws[,1] / rho) * params$gamma) %>%
  mutate(index = 1:200)
# --------------------------------------------------
# run simulation over many parameter values
function_wrapper <- function(x, contact_matrix, init, t){
  params$beta <- x[1]
  params$c_start <- contact_matrix[[x[2]]]
  seir_out <- lsoda(init, t, age_struct_seir_ode, params) #
  seir_out <- as.data.frame(seir_out)
  out_mle <- postprocess_age_struct_model_output(seir_out)
  daily_cases <- params$sigma * rowSums(out_mle$E + out_mle$Ev_1d + out_mle$Ev_2d) * params$p_report
  return(daily_cases)
}
# run model for each combination of parameters
out <- apply(betas, 1, function_wrapper, contact_matrix = baseline_2017, 
             init = init, t = time_vec) # rows are time points, columns are different simulations
# --------------------------------------------------

# plot with confidence bounds
bounds <- apply(out,1,function(x) quantile(x, c(0.025,0.975)))

plot(osiris2$inc ~ time_vec, col = "red", pch = 16, 
     xlab = "Time (days)",ylab = "Daily Cases"
     #, ylim = c(0,700)
     ) 
lines(daily_cases_mle,col = "blue")
lines(bounds[1,], col = "blue", lty = 2, lwd = 0.5)
lines(bounds[2,], col = "blue", lty = 2, lwd = 0.5)
legend("topright", c("Data","Model","95% credible intervals"),
       col = c("red","blue","blue"), lty = c(0,1,2), pch = c(16,NA,NA))

# --------------------------------------------------
# re-run mle for each set of time points
# --------------------------------------------------

breakpoints_2020 <- list( 
  date = c( 
    as.Date("2020-03-16"),  # closure of hospitality, schools, daycares
    as.Date("2020-03-31"),  # rest of march until we use April 2020 contact matrices
    as.Date("2020-05-11"),  # reopening of primary schools and daycare @ 50% capacity
    as.Date("2020-06-01"),  # cafes/restaurants reopen with max 30 visitor,
                            # primary schools and daycares open at 100% 
    as.Date("2020-07-01"),  # secondary schools reopen @ 100% capacity
    as.Date("2020-09-29"),  # max 3 visitors at home, max group size 30, face masks in public areas
    as.Date("2020-10-20"),  # cafes/restaurants close, no team sports, 3 visitors at home per day
    as.Date("2020-12-01"),  # masks mandatory in all public and indoor areas
    as.Date("2020-12-15"),  # non-essential shops close
    as.Date("2020-12-31")), # end of year
    
  contact_matrix = list( baseline_2017, baseline_2017, 
                         april_2020, april_2020, june_2020, 
                         june_2020, september_2020,september_2020,
                         september_2020)
)

n_bp <- length(breakpoints_2020$date)
mles <- matrix(rep(NA, 2*n_bp), nrow = n_bp)
colnames(mles) <- c("beta", "alpha")
out_mle <- list()
bounds_list <- list()
daily_cases_mle <- list()

for (j in 1:n_bp) {
  
  # set contact matrix for time window
  if (j == n_bp){
    params$c_start <- breakpoints_2020$contact_matrix[[j-1]]$mean
  } else {
    params$c_start <- breakpoints_2020$contact_matrix[[j]]$mean
  }
  
  if (j == 1) {
    # if first time window, start time at 0
    end_day <- yday(breakpoints_2020$date[j]) - 1
    times <- seq(0, end_day, by = 1)
    # set initial conditions to those specified earlier in script
    init_update <- init
    pars <- c(2.3, 0.01)
    S_diag <- diag(init_update[c(2:10)])
    rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  } else {
    start_day <- yday(breakpoints_2020$date[j-1]) - 1 # shift days by 1 because we start time at 0 (not 1)
    end_day <- yday(breakpoints_2020$date[j]) - 1
    times <- seq(start_day, end_day, by = 1)
    # update initial conditions based on last time window
    init_update <- c(
      t = times[1],
      S = as.numeric(tail(out_mle[[j - 1]]$S, 1)),
      Shold_1d = as.numeric(tail(out_mle[[j - 1]]$Shold_1d, 1)),
      Sv_1d = as.numeric(tail(out_mle[[j - 1]]$Sv_1d, 1)),
      Shold_2d = as.numeric(tail(out_mle[[j - 1]]$Shold_2d, 1)),
      Sv_2d = as.numeric(tail(out_mle[[j - 1]]$Sv_2d, 1)),
      E = as.numeric(tail(out_mle[[j - 1]]$E, 1)),
      Ev_1d = as.numeric(tail(out_mle[[j - 1]]$Ev_1d, 1)),
      Ev_2d = as.numeric(tail(out_mle[[j - 1]]$Ev_2d, 1)),
      I = as.numeric(tail(out_mle[[j - 1]]$I, 1)),
      Iv_1d = as.numeric(tail(out_mle[[j - 1]]$Iv_1d, 1)),
      Iv_2d = as.numeric(tail(out_mle[[j - 1]]$Iv_2d, 1)),
      H = as.numeric(tail(out_mle[[j - 1]]$H, 1)),
      Hv_1d = as.numeric(tail(out_mle[[j - 1]]$Hv_1d, 1)),
      Hv_2d = as.numeric(tail(out_mle[[j - 1]]$Hv_2d, 1)),
      H_IC = as.numeric(tail(out_mle[[j - 1]]$H_IC, 1)),
      H_ICv_1d = as.numeric(tail(out_mle[[j - 1]]$H_ICv_1d, 1)),
      H_ICv_2d = as.numeric(tail(out_mle[[j - 1]]$H_ICv_2d, 1)),
      IC = as.numeric(tail(out_mle[[j - 1]]$IC, 1)),
      ICv_1d = as.numeric(tail(out_mle[[j - 1]]$ICv_1d, 1)),
      ICv_2d = as.numeric(tail(out_mle[[j - 1]]$ICv_2d, 1)),
      D = as.numeric(tail(out_mle[[j - 1]]$D, 1)),
      R = as.numeric(tail(out_mle[[j - 1]]$R, 1)),
      Rv_1d = as.numeric(tail(out_mle[[j - 1]]$Rv_1d, 1)),
      Rv_2d = as.numeric(tail(out_mle[[j - 1]]$Rv_2d, 1))
    )
    pars <- c((mles[j-1,1]/params$gamma)*rho, mles[j-1,2])
    S_diag <- diag(init_update[c(2:10)])
    rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  }
  
  # subset data for time window
  osiris_sub <- osiris1[times + 1, ]
  
  # optimize
  res <- optim(par = pars, 
               fn = likelihood_func,
               method = "L-BFGS-B",
               lower = c(0,0.0001),
               upper = c(10,1),
               t = times,
               data = osiris_sub,
               params = params,
               init = init_update,
               stochastic = FALSE,
               hessian = TRUE
  )
  
  # store MLE
  mles[j,1] <- (res$par[1] / rho) * params$gamma
  mles[j,2] <- res$par[2]
  
  # draw 200 parameter values
  parameter_draws_mle <- mvtnorm::rmvnorm(200, res$par, solve(res$hessian))
  betas <- data.frame(beta = (parameter_draws_mle[,1] / rho) * params$gamma) %>%
    mutate(index = 1:200)
  # --------------------------------------------------
  # run model for each combination of parameters
  out <- apply(betas, 1, function_wrapper, 
               contact_matrix = breakpoints_2020$contact_matrix[[j]],
               init = init_update, t = times) # rows are time points, columns are different simulations
  
  # get confidence bounds of model runs
  bounds_list[[j]] <- apply(out,1,function(x) quantile(x, c(0.025,0.975)))
  # --------------------------------------------------
  # re-run simulation with MLE parameter values
  # --------------------------------------------------
  # transform R0 to beta
  params$beta <- mles[j,1]
  if (j == n_bp){
    params$c_start <- breakpoints_2020$contact_matrix[[j-1]]$mean
  } else {
    params$c_start <- breakpoints_2020$contact_matrix[[j]]$mean
  }
  # run simulation
  seir_out <- lsoda(init_update, times, age_struct_seir_ode, params)
  seir_out <- as.data.frame(seir_out)
  out_mle[[j]] <- postprocess_age_struct_model_output(seir_out)
  daily_cases_mle[[j]] <- params$sigma * rowSums(out_mle[[j]]$E + out_mle[[j]]$Ev_1d + out_mle[[j]]$Ev_2d) * params$p_report
  
  # --------------------------------------------------
  # plot
  plot(osiris_sub$inc ~ times, col = "red", pch = 16, 
       xlab = "Time (days)",ylab = "Daily Cases" #, ylim = c(0,5000)
  ) 
  lines(times, daily_cases_mle[[j]],col = "blue")
  lines(times, bounds_list[[j]][1,], col = "blue", lty = 2, lwd = 0.5)
  lines(times, bounds_list[[j]][2,], col = "blue", lty = 2, lwd = 0.5)
  legend("topright", c("Data","Model","95% credible intervals"),
         col = c("red","blue","blue"), lty = c(0,1,2), pch = c(16,NA,NA))
  
} # end of for loop over breakpoints

last_date_in_osiris <- "2021-05-25"
# saveRDS(out_mle, file = paste0("output_from_fits_", last_date_in_osiris, ".rds"))
# saveRDS(daily_cases_mle, file = paste0("cases_from_fits_", last_date_in_osiris, ".rds"))

saveRDS(out_mle, file = "output_2020_from_model_fit.rds")
saveRDS(daily_cases_mle, file = "cases_2020_from_model_fit.rds")

# --------------------------------------------------
#  combine all piecewise results to plot together
cases_all <- unlist(daily_cases_mle)
cases_all1 <- unique(cases_all) # remove duplicate time points
times_all <- 1:length(cases_all1)
model_fit <- data.frame(time = times_all, cases = cases_all1)
saveRDS(model_fit, file = paste0("model_fit_df_", last_date_in_osiris, ".rds"))

plot(osiris1$inc[times_all] ~ times_all,
  col = "red", pch = 16,
  xlab = "Time (days)", ylab = "Incidence"#, ylim = c(0, 7600)
)
lines(times_all, cases_all1, col = "blue")
legend("bottomright", c("Osiris Data", "Model Fit"),
  col = c("red", "blue"), lty = c(0, 1), pch = c(16, NA)
)
# --------------------------------------------------

# --------------------------------------------------
# save initial conditions for forward simulation
init_forward <- c(
  t = tail(times_all, 1) - 1,
  S = as.numeric(tail(out_mle[[8]]$S, 1)),
  Shold_1d = as.numeric(tail(out_mle[[8]]$Shold_1d, 1)),
  Sv_1d = as.numeric(tail(out_mle[[8]]$Sv_1d, 1)),
  Shold_2d = as.numeric(tail(out_mle[[8]]$Shold_2d, 1)),
  Sv_2d = as.numeric(tail(out_mle[[8]]$Sv_2d, 1)),
  E = as.numeric(tail(out_mle[[8]]$E, 1)),
  Ev_1d = as.numeric(tail(out_mle[[8]]$Ev_1d, 1)),
  Ev_2d = as.numeric(tail(out_mle[[8]]$Ev_2d, 1)),
  I = as.numeric(tail(out_mle[[8]]$I, 1)),
  Iv_1d = as.numeric(tail(out_mle[[8]]$Iv_1d, 1)),
  Iv_2d = as.numeric(tail(out_mle[[8]]$Iv_2d, 1)),
  H = as.numeric(tail(out_mle[[8]]$H, 1)),
  Hv_1d = as.numeric(tail(out_mle[[8]]$Hv_1d, 1)),
  Hv_2d = as.numeric(tail(out_mle[[8]]$Hv_2d, 1)),
  H_IC = as.numeric(tail(out_mle[[8]]$H_IC, 1)),
  H_ICv_1d = as.numeric(tail(out_mle[[8]]$H_ICv_1d, 1)),
  H_ICv_2d = as.numeric(tail(out_mle[[8]]$H_ICv_2d, 1)),
  IC = as.numeric(tail(out_mle[[8]]$IC, 1)),
  ICv_1d = as.numeric(tail(out_mle[[8]]$ICv_1d, 1)),
  ICv_2d = as.numeric(tail(out_mle[[8]]$ICv_2d, 1)),
  D = as.numeric(tail(out_mle[[8]]$D, 1)),
  R = as.numeric(tail(out_mle[[8]]$R, 1)),
  Rv_1d = as.numeric(tail(out_mle[[8]]$Rv_1d, 1)),
  Rv_2d = as.numeric(tail(out_mle[[8]]$Rv_2d, 1))
)
saveRDS(init_forward, file = paste0("init_conditions_", last_date_in_osiris, ".rds"))

