# Test script for stochastic model version

# Create list of parameter values for input into model solver
params <- list(dt = 1/6,                  # units, 1 = day
               beta =  0.0004,            # transmission rate, units: per day
               beta1 = 0.14,              # amplitude of seasonal forcing
               gamma = g,                 # rate of becoming infectious, units: per day
               sigma = s,                 # rate of becoming infected, units: per day
               epsilon = 0.01,            # import case rate, units: per day
               #mu = 0.00003653,          # birth rate, units: per day
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
               keep_cm_fixed = TRUE,
               vac_inputs = NULL,
               use_cases = TRUE,                           # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 35.7/100000 * sum(n_vec),      #35.7  # 20 for IC admissions
               no_vac = FALSE,
               t_calendar_start = yday(as.Date("2021-01-31")),   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
               breakpoints = NULL  # breakpoints - start_date    # time points when parameters can change (if NULL, then beta is constant over time)
)

times <- seq(0, 200, by = 1)

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

# determine transmission rate (beta) for r0 ------------------------
r0 <- 2.3
S_diag <- diag(c(rep(100000, 9)))
rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
beta <- (r0 / rho) * params$gamma
# check
K <- (1 / params$gamma) * beta * S_diag %*% params$c_start
as.numeric(eigs(K, 1)$values) # this should be r0
params$beta <- beta

# single simulation run --------------------------------------------
# if time doesn't start at 0 we need to initialise the contact 
# matrices flags
if(times[1] != 0){
  flag_relaxed <- 0 # start with relaxed contact matrix
  flag_very_relaxed <- 0
  flag_normal <- 0
}
# Solve model ------------------------------------------------------
seir_out <- stochastic_age_struct_seir_ode(times,init,params) #
#seir_out <- as.data.frame(seir_out)
#out <- postprocess_age_struct_model_output(seir_out)

# Plot -------------------------------------------------------------
#seirds_col <- c("#8c8cd9", "#e67300", "#d279a6", "#ff4d4d", "#999966", "#660000")
my_cols <- viridis_pal(option = "D")(7)
x_res <- apply(seir_out, 3, rowSums)
x_times <- as.numeric(str_remove(rownames(x_res), "[t]"))
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(x_times, x_res[,c("S","E","I","H", "IC", "D", "R")], xlab = "Time", ylab = "Number of individuals",
        type = "l", col = my_cols, lty = 1)
legend("left", lwd = 1, col = my_cols, legend = ,c("S","E","I","H", "IC", "D", "R"), bty = "n")

