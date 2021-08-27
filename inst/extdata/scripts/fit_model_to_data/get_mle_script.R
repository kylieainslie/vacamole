# find mle of beta from SEIR model fit to real data

# --------------------------------------------------
# load packages
#library(segmented)
library(ggplot2)
library(tidyverse)
library(lubridate)

# Read in OSIRIS data ------------------------------
source("inst/extdata/scripts/model_run_helper.R")
source("~/vacamole/R/model_run_wrapper.R")
# read in OSIRIS data
osiris <- readRDS("inst/extdata/data/real_data/Osiris_Data_20210730_1043.rds")

osiris1 <- osiris %>%
  #filter(!is.na(date)) %>%
  #complete(date = seq.Date(min(date), max(date), by="day"), fill = list(inc = 0)) %>%
  mutate(roll_avg = zoo::rollmean(inc, k = 7, fill = 0)) %>%
  filter(date < max(date)-2) # remove last 3 days due to reporting delay

# plot data
p <- ggplot(osiris1, aes(x = date, y = inc)) +
  geom_line() +
  geom_line(aes(x = date, y = roll_avg, color = "red")) +
  theme(panel.background = element_blank())
p

#time_vec <- seq(0, 76, by = 1)
# --------------------------------------------------
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


# specify model parameters
params <- list(dt = 1,
               beta = 0.0003934816,       # transmission rate
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
               vac_inputs = basis1,
               use_cases = TRUE,                           # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 35.7/100000 * sum(n_vec),        # 35.7  # 20 for IC admissions
               no_vac = FALSE,
               t_calendar_start = yday(as.Date("2020-01-01"))   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
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
  lik
}
# ---------------------------------------------------
# Determine MLE using optim

# res <- optim(par = c(2.3, 0.01), 
#              fn = likelihood_func,
#              method = "L-BFGS-B",
#              lower = c(1,0.001),
#              upper = c(10,1),
#              t = time_vec,
#              data = osiris2,
#              params = params,
#              init = init,
#              stochastic = FALSE,
#              hessian = TRUE
#              )
# 
# res$par

# plot to check fit --------------------------------
# S_diag <- diag(init[c(2:10)])
# rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
# params$beta <- (res$par[1] / rho) * params$gamma
# #params$beta <- res$par[1]
# seir_out <- lsoda(init, time_vec, age_struct_seir_ode, params) #
# seir_out <- as.data.frame(seir_out)
# out_mle <- postprocess_age_struct_model_output(seir_out)
# daily_cases_mle <- params$sigma * rowSums(out_mle$E + out_mle$Ev_1d + out_mle$Ev_2d) * params$p_report
# 
# plot(daily_cases_mle ~ time_vec, type = "l") # ylim = c(0, 8000)
# points(osiris2$inc ~ time_vec, col = "red", pch = 16)

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

breakpoints <- list( 
  date = c( 
    as.Date("2020-03-16"),  # closure of hospitality, schools, daycares
    as.Date("2020-03-24"),  # casinos subject to same measures as food/beverage outlets
    as.Date("2020-04-29"),  # Children and young people can participate in outdoor activities
    as.Date("2020-05-11"),  # reopening of primary schools and daycare @ 50% capacity
    as.Date("2020-06-01"),  # cafes/restaurants reopen with max 30 visitor,
                            # primary schools and daycares open at 100% 
    as.Date("2020-07-01"),  # secondary schools reopen @ 100% capacity
    as.Date("2020-08-06"),  # Small groups allowed at universities
    as.Date("2020-08-18"),  # maximum 6 people in household
    as.Date("2020-09-20"),  # groups <= 50 people, some hospitality closes at midnight or 1 AM 
    as.Date("2020-09-29"),  # max 3 visitors at home, max group size 30, face masks in public areas
    as.Date("2020-10-14"),  # cafes/restaurants close, no team sports, 3 visitors at home per day
    as.Date("2020-10-23"),  # hotels can't sell alcohol after 20:00
    as.Date("2020-11-04"),  # 2 visitors per day, groups <= 20
    as.Date("2020-11-19"),  # 3 visitors per day, groups <= 30
    as.Date("2020-12-01"),  # masks mandatory in all public and indoor areas
    as.Date("2020-12-15"),  # non-essential shops close
    as.Date("2021-01-01"),  # end of year
    as.Date("2021-01-20"),  # 1 visitor per day/curfew (23/01/2021)
    as.Date("2021-02-08"),  # Primary schools, child care, special ed reopen
    #as.Date("2021-02-15"),  # Alpha becomes dominant variant
    as.Date("2021-03-01"),  # secondary schools partially reopen, contact professions reopen (3/3/2021)
    as.Date("2021-03-16"),  # some retail reopens
    as.Date("2021-03-31"),  # curfew to start at 22:00 instead of 21:00
    as.Date("2021-04-19"),  # out of school care fully reopens
    as.Date("2021-04-28"),  # curfew canceled
    as.Date("2021-05-19"),  # <27 can play outdoor sports, groups <= 30, non-essential travel allowed within NL 
    as.Date("2021-06-05"),  # 4 visitors per day, museums reopen, group <= 50, restaurants reopen
    as.Date("2021-06-22"),  # vaccination of 12-17 year olds starts, Delta becomes dominant strain
    as.Date("2021-06-26"),  # all restrictions relaxed, except masks on public transport, nightclubs reopen
    as.Date("2021-07-10"),  # catering industry reopens, test for entry with large events, nightclubs close
    as.Date("2021-07-19"),  # work from home advisory re-instated
    as.Date("2021-07-27")   # last date in osiris
  ),  
  contact_matrix = list( baseline_2017, 
                         baseline_2017, 
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
                         #february_2021,
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
                         june_2021
                         ),
  indicator_2021 = c(rep(0,16), rep(1,16)) #,
  #p_report = c(rep(0.1, 6), rep(0.33, 24)) # case ascertainment lower in first wave
)

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
    if (breakpoints$indicator_2021[j] == 1){
      if(breakpoints$indicator_2021[j-1] == 1){ # wait for two consecutive dates in 2021
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
  osiris_sub <- osiris1[times + 1, ]
  
  # optimize
  res <- optim(par = pars, 
               fn = likelihood_func,
               method = "L-BFGS-B",
               lower = c(0,0.005),
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
  parameter_draws[[j]] <- mvtnorm::rmvnorm(200, res$par, solve(res$hessian))
  beta_draws[[j]] <- data.frame(beta = (parameter_draws[[j]][,1] / rho) * params$gamma) %>%
     mutate(index = 1:200)
  # --------------------------------------------------

} # end of for loop over breakpoints

# save outputs
saveRDS(mles, file = paste0(path, "mles_from_fits_", todays_date, ".rds"))
saveRDS(beta_draws, file = paste0(path, "beta_draws_from_fits_", todays_date, ".rds"))
# name list elements for easier indexing
# names(out_mle) <- paste0("end_date_", breakpoints$date)
# names(daily_cases_mle) <- paste0("end_date_", breakpoints$date)

# ----------------------------------------------------
# run simulations for mle, lower, and upper bounds 
# of beta
# ----------------------------------------------------
beta_mles <- readRDS("inst/extdata/results/model_fits/mles_from_fits_2021-08-25.rds")
beta_mles_list <- split(beta_mles, seq(nrow(beta_mles)))
beta_draws <- readRDS("inst/extdata/results/model_fits/beta_draws_from_fits_2021-08-25.rds")
# beta_dat <- bind_rows(lapply(beta_draws, function(x)quantile(x[,1], probs = c(0.025, 0.975)))) %>%
#   mutate(mle = beta_mles[,1]) %>%
#   rename(lower = `2.5%`,
#          upper = `97.5%`)


mle_run <- model_run_wrapper(breakpoints = breakpoints, beta_values = beta_mles_list, init_conditions = init, params = params)
ci_run  <- model_run_wrapper(breakpoints = breakpoints, beta_values = beta_draws, init_conditions = init, params = params, mle = FALSE)
ci_out <- list()
for (i in 1:31){
  ci_out[[i]] <- do.call("rbind", ci_run[[i]])
}
ci_out_wide <- do.call("cbind", ci_out)
matplot(t(ci_out_wide), type = "l")

# save outputs -------------------------------------
path <- "inst/extdata/results/model_fits/"
todays_date <- Sys.Date()

# --------------------------------------------------
#  combine all piecewise results to plot together
cases_mle <- unique(unlist(mle_run))
cases_lower <- unique(unlist(lower_run))
cases_upper <- unique(unlist(upper_run))
times_all <- 1:length(cases_mle)

model_fit <- data.frame(time = times_all, date = osiris1$date, real = osiris1$inc, mle = cases_mle, lower = cases_lower, upper = cases_upper)
# saveRDS(model_fit, file = paste0("model_fit_df_", todays_date, ".rds"))

p <- ggplot(data = model_fit, aes(x = date, y = upper)) +
  geom_line() +
  #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_point(data = model_fit, aes(x = date, y = real), color = "red")

p
# --------------------------------------------------
# old code
# init_update <- c(t = times[1],
                 # S = as.numeric(tail(out_mle[[j - 1]]$S, 1)),
                 # Shold_1d = as.numeric(tail(out_mle[[j - 1]]$Shold_1d, 1)),
                 # Sv_1d = as.numeric(tail(out_mle[[j - 1]]$Sv_1d, 1)),
                 # Shold_2d = as.numeric(tail(out_mle[[j - 1]]$Shold_2d, 1)),
                 # Sv_2d = as.numeric(tail(out_mle[[j - 1]]$Sv_2d, 1)),
                 # E = as.numeric(tail(out_mle[[j - 1]]$E, 1)),
                 # Ev_1d = as.numeric(tail(out_mle[[j - 1]]$Ev_1d, 1)),
                 # Ev_2d = as.numeric(tail(out_mle[[j - 1]]$Ev_2d, 1)),
                 # I = as.numeric(tail(out_mle[[j - 1]]$I, 1)),
                 # Iv_1d = as.numeric(tail(out_mle[[j - 1]]$Iv_1d, 1)),
                 # Iv_2d = as.numeric(tail(out_mle[[j - 1]]$Iv_2d, 1)),
                 # H = as.numeric(tail(out_mle[[j - 1]]$H, 1)),
                 # Hv_1d = as.numeric(tail(out_mle[[j - 1]]$Hv_1d, 1)),
                 # Hv_2d = as.numeric(tail(out_mle[[j - 1]]$Hv_2d, 1)),
                 # H_IC = as.numeric(tail(out_mle[[j - 1]]$H_IC, 1)),
                 # H_ICv_1d = as.numeric(tail(out_mle[[j - 1]]$H_ICv_1d, 1)),
                 # H_ICv_2d = as.numeric(tail(out_mle[[j - 1]]$H_ICv_2d, 1)),
                 # IC = as.numeric(tail(out_mle[[j - 1]]$IC, 1)),
                 # ICv_1d = as.numeric(tail(out_mle[[j - 1]]$ICv_1d, 1)),
                 # ICv_2d = as.numeric(tail(out_mle[[j - 1]]$ICv_2d, 1)),
                 # D = as.numeric(tail(out_mle[[j - 1]]$D, 1)),
                 # R = as.numeric(tail(out_mle[[j - 1]]$R, 1)),
                 # Rv_1d = as.numeric(tail(out_mle[[j - 1]]$Rv_1d, 1)),
                 # Rv_2d = as.numeric(tail(out_mle[[j - 1]]$Rv_2d, 1))
# )

