
# Create list of parameter values for input into model solver
params <- list(dt = 1,
               beta =  0.0004,            # transmission rate
               beta1 = 0.14,              # amplitude of seasonal forcing
               gamma = g,                 # 1/gamma = infectious period
               sigma = s,                 # 1/sigma = latent period
               delta = NULL,              # scaling constant for beta (if NULL it is excluded)
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
               c_start = t1,
               c_lockdown = t1,
               c_relaxed = t1,
               c_very_relaxed = t1,
               c_normal = t1,
               keep_cm_fixed = FALSE,
               vac_inputs = basis1_no_wane,
               use_cases = TRUE,                           # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 35.7/100000 * sum(n_vec),      #35.7  # 20 for IC admissions
               no_vac = FALSE,
               t_calendar_start = yday(as.Date("2021-01-31")),   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
               breakpoints = NULL  # breakpoints - start_date    # time points when parameters can change (if NULL, then beta is constant over time)
)

times <- seq(0, 40, by = 1)

# Specify initial values -------------------------------------------
empty_state <- c(rep(0, 9))
init <- c(
  t = 0,
  S = c(rep(100000, 9)),
  Shold_1d = empty_state,
  Sv_1d = empty_state,
  Shold_2d = empty_state,
  Sv_2d = empty_state,
  E = empty_state,
  Ev_1d = empty_state,
  Ev_2d = empty_state,
  I = c(rep(10,9)),
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

# single simulation run --------------------------------------------
# if time doesn't start at 0 we need to initialise the contact 
# matrices flags
if(times[1] != 0){
  flag_relaxed <- 0 # start with relaxed contact matrix
  flag_very_relaxed <- 0
  flag_normal <- 0
}
# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times,stochastic_age_struct_seir_ode,params) #
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)
