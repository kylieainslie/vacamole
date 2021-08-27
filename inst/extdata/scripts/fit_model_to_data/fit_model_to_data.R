# Fit model to data to get appropriate initial conditions ---------------

# To run this script the pre-amble from the run_age_struct_seir.R 
# MUST be run! ----------------------------------------------------------
# Load packages ---------------------------------------------------------
library(lubridate)
library(FSA)
library(zoo)

# read in hospital admission counts -------------------------------------
nice <- readRDS("inst/extdata/data/nice_admissions_20210407.rds")
nice_old <- readRDS("inst/extdata/data/nice_admissions_20210212.rds") %>%
  complete(., AdmissionDate = full_seq(AdmissionDate, 1)) %>%
  mutate(n = ifelse(is.na(n), 0, n),
         moving_avg3 = zoo::rollmean(n, k = 3, fill = NA))
  
nice1 <- nice %>%
  group_by(AdmissionDate, IsIc) %>%
  summarise_at(.vars = "n", .funs= sum) %>%
  filter(IsIc == 0) %>%
  # make dates consecutive
  complete(., AdmissionDate = full_seq(AdmissionDate, 1)) %>%
  mutate(n = ifelse(is.na(n), 0, n),
         n = as.double(n))
nice1$moving_avg3 <- zoo::rollmean(nice1$n, k = 3, fill = NA)
nice1$moving_avg7 <- zoo::rollmean(nice1$n, k = 7, fill = NA)

plot(nice1$AdmissionDate, nice1$moving_avg7, type = "l")

# initial states ---------------------------------------------------
# Jacco's suggested way to determine initial conditions
init_states_dat <- data.frame(age_group = c("0-9", "10-19", "20-29", "30-39", "40-49", 
                                            "50-59", "60-69", "70-79","80+"),
                              n = n_vec,
                              # from Scott
                              n_recovered = c(30720, 397100, 642600, 419000, 412200, 505900,
                                              349100, 206800, 115900 + 33200),
                              # from sitrep for 26 januari tot 2 februari: 
                              # https://www.rivm.nl/coronavirus-covid-19/actueel/wekelijkse-update-epidemiologische-situatie-covid-19-in-nederland)
                              n_cases = c(835, 2851, 4591, 3854, 3925, 5191, 3216, 1819, 
                                          1376 + 485),
                              # from NICE data (n_hosp/n_ic refers to occupancy on 1 Feb 2021)
                              n_hosp = c(2, 1, 8, 19, 29, 76, 142, 159, 186),
                              n_ic = c(0, 2, 6, 9, 25, 83, 181, 150, 12),
                              # from NICE data: people with length of stay >= 9 days
                              n_hosp_after_ic = c(2, 2, 9, 15, 45, 158, 321, 392, 266) 
) %>%
  mutate(n_infections = n_cases * 3, #p_reported_by_age,
         init_E = n_infections * (2/7),
         init_I = n_infections * (2/7),
         init_S = n - n_recovered - init_E - init_I  - n_hosp - n_ic - n_hosp_after_ic)

# Specify initial values -------------------------------------------
empty_state <- c(rep(0,9))
init <- c(t = 0,                  
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

# Input parameters -------------------------------------------------
params <- list(beta = 0.0006,                  # transmission rate
               gamma = 0.5,                    # 1/gamma = infectious period
               sigma = 0.5,                    # 1/sigma = latent period
               delta = 1,                      # scaling constant 
               N = n_vec,                      # Population (no need to change)
               h = h,                          # Rate from infection to hospital admission/ time from infection to hosp admission
               i1 = i1,
               i2 = i2,
               d = d, 
               d_ic = d_ic,
               d_hic = d_hic,
               r = r,
               r_ic = r_ic,
               p_report = 1/3, #p_reported_by_age,
               #c_start = B_prime,
               c_lockdown = t5,
               c_relaxed = t4,
               c_very_relaxed = t3,
               c_normal = t1,
               vac_schedule = basis1,
               ve = ve,
               delay = delays,
               hosp_multiplier = h_multiplier,
               use_cases = TRUE,              # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 35.7/100000 * sum(n_vec),        # 20 for IC admissions
               no_vac = TRUE
               )


# write wrapper function for seir model code that outputs only hospital
# admission counts ------------------------------------------------------
fit_to_data_wrapper <- function(x, params, model_params, init_states, init_states_dat){
 
  reff_jan <- params[1]
  delta_jan <- params[2]
  reff_feb <- params[3]
  delta_feb <- params[4]
  reff_march <- params[5]
  delta_march <- params[6]
  
  # list of dates at different intervention points
  t_vec <- lubridate::yday(as.Date(c("2021-01-05",
                                     "2021-02-01"#,
                                     #"2021-03-01",
                                     #"2021-04-01"
                                     )))
  

  # loop over time periods -------------------------------------------
  rtn <- list() # store output for each iteration here
  for (i in 1:(length(t_vec)-1)){
    times <- seq(t_vec[i], t_vec[i+1], by = 1)
    cat("Time window", i, "\n")
    delta <- (i == 1) * delta_jan + (i == 2) * delta_feb + (i == 3) * delta_march
    reff <- (i == 1) * reff_jan + (i == 2) * reff_feb + (i == 3) * reff_march
    
    # determine transmission rate
    S2 = diag(init_states_dat$init_S)
    rho2 = as.numeric(eigs(S2 %*% model_params$c_lockdown,1)$values)
    beta2 = reff/rho2 * model_params$gamma
    # check
    B <- model_params$c_lockdown
    K2 = beta2 * (1/model_params$gamma) * S2 %*% B
    as.numeric(eigs(K2,1)$values) # this should be reff
    
    # callibrate distribution of cases across age groups ----------------
    w <- eigen(K2)$vectors[,1]/sum(eigen(K2)$vectors[,1]) # should match dist_cases (I think!)
    x <- init_states_dat$n_cases / sum(init_states_dat$n_cases)
    A <- diag(x/w)
    B_prime <- A %*% B
    rho2_prime = as.numeric(eigs(S2 %*% B_prime,1)$values)
    beta2_prime = reff/rho2_prime * model_params$gamma
    K2_prime <- beta2_prime * (1/model_params$gamma) * S2 %*% B_prime
    dom_eig_vec <- eigen(K2_prime)$vectors[,1]
    w_prime <- dom_eig_vec/sum(dom_eig_vec)
    as.numeric(eigs(K2_prime,1)$values)
    
    model_params$beta <- beta2_prime
    model_params$c_lockdown <- B_prime
    
  # Specify initial values -------------------------------------------
  if (i == 1){
    init <- init_states
  } else {init <- init_new}
    
  # Solve model ------------------------------------------------------
  cat("Solving model...","\n")
  seir_out <- lsoda(init,times,age_struct_seir_ode,model_params)
  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output(seir_out)
  # update initial states for next time window -----------------------
  init_new <- c(t = dim(out$S)[1],
                S = as.numeric(tail(out$S,1)),
                Shold_1d = as.numeric(tail(out$Shold_1d,1)),
                Sv_1d = as.numeric(tail(out$Sv_1d,1)),
                Shold_2d = as.numeric(tail(out$Shold_2d,1)),
                Sv_2d = as.numeric(tail(out$Sv_2d,1)),
                E = as.numeric(tail(out$E,1)),
                Ev_1d = as.numeric(tail(out$Ev_1d,1)),
                Ev_2d = as.numeric(tail(out$Ev_2d,1)),
                I = as.numeric(tail(out$I,1)),
                Iv_1d = as.numeric(tail(out$Iv_1d,1)),
                Iv_2d = as.numeric(tail(out$Iv_2d,1)),
                H = as.numeric(tail(out$H,1)),
                Hv_1d = as.numeric(tail(out$Hv_1d,1)),
                Hv_2d = as.numeric(tail(out$Hv_2d,1)),
                H_IC = as.numeric(tail(out$H_IC,1)),
                H_ICv_1d = as.numeric(tail(out$H_ICv_1d,1)),
                H_ICv_2d = as.numeric(tail(out$H_ICv_2d,1)),
                IC = as.numeric(tail(out$IC,1)),
                ICv_1d = as.numeric(tail(out$ICv_1d,1)),
                ICv_2d = as.numeric(tail(out$ICv_2d,1)),
                D = as.numeric(tail(out$D,1)),
                R = as.numeric(tail(out$R,1)),
                Rv_1d = as.numeric(tail(out$Rv_1d,1)),
                Rv_2d = as.numeric(tail(out$Rv_2d,1))
  )
  #print(init_i)
  hosp_by_age_group <- sweep((out$I + out$Iv_1d + out$Iv_2d), 2, h, "*")
  out_vec <- rowSums(hosp_by_age_group)
  rtn[[i]] <- out_vec
  }
  
  rtn2 <- unlist(rtn)
  
  return(rtn2)
}

# use nls to fit non-linear least squares ----------------------------
# parameter initial values
init_reff <- c(1, 1, 1)
init_delta <- c(0.9, 1, 1.1)
init_params <- c(init_reff[1], init_delta[1],
                 init_reff[2], init_delta[2],
                 init_reff[3], init_delta[3])

lower_reff <- c(0.75, 0.75, 0.75)
lower_delta <- c(0, 0, 0)
lower_params <- c(lower_reff[1], lower_delta[1],
                 lower_reff[2], lower_delta[2],
                 lower_reff[3], lower_delta[3])

upper_reff <- c(2, 2, 2)
upper_delta <- c(2, 2, 2)
upper_params <- c(upper_reff[1], upper_delta[1],
                  upper_reff[2], upper_delta[2],
                  upper_reff[3], upper_delta[3])

# fit model
nice1_sub <- nice1 %>%
  filter(!is.na(moving_avg7),
         AdmissionDate < as.Date("2021-02-02"))

model_fit <- nls(moving_avg7 ~ fit_to_data_wrapper(AdmissionDate, params, model_params, init, init_states_dat), 
                 data = nice1_sub, # where x and y are
                 start = list(params = init_params), #initial conditions
                 lower = lower_params, # lower bound of parameter values
                 upper = upper_params, # upper bound of parameter values
                 trace = T,
                 algorithm = "port", # must be used to specify bounds
                 control = list(#maxiter = 100000, 
                                warnOnly=T)
                 )
# plot ---------------------------------------------------------------
# loess_fit <- loess(n ~ time, data=nice1, span=0.5)
# smoothed <- predict(loess_fit) 
ymax <- max(max(model_fit$m$fitted()),nice_sub$moving_avg3)
plot(moving_avg3~AdmissionDate,data=nice_sub,type = "l", ylim = c(0,ymax),
     xlab = "Date", ylab = "Hospital Admissions") 
lines(x = nice_sub$AdmissionDate, y = model_fit$m$fitted(), 
      type = "l", col = "red")
# lines(smoothed, x = nice1$AdmissionDate, col = "blue")

