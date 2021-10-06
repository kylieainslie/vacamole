# Test script for stochastic model version

# Create list of parameter values for input into model solver
params <- list(dt = 1,                  # units, 1 = day
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
               t_calendar_start = yday(as.Date("2020-01-01")),   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
               breakpoints = NULL  # breakpoints - start_date    # time points when parameters can change (if NULL, then beta is constant over time)
)

times <- seq(0, 100, by = 1)

# Specify initial values -------------------------------------------
empty_state <- c(rep(0, 9))
init <- c(
  t = 0,
  S = c(rep(10000, 9)),
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

# determine transmission rate (beta) for r0 ------------------------
r0 <- 2.3
S_diag <- diag(c(rep(1000, 9)))
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
my_cols <- viridis_pal(option = "D")(4)
x_res <- apply(seir_out, 3, rowSums)
x_times <- as.numeric(str_remove(rownames(x_res), "[t]"))
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(x_times, x_res[,c("S","E","I","R")], xlab = "Time", ylab = "Number of individuals",
        type = "l", col = my_cols, lty = 1)
legend("left", lwd = 1, col = my_cols, legend = ,c("S","E","I","R"), bty = "n")

# Multiple runs ----------------------------------------------------
library(foreach)
library(doParallel)
library(parallel)

n_cores <- detectCores()
registerDoParallel(n_cores)
n_sims <- 100

mult_sim_out <- foreach (i = 1:n_sims) %dopar% {
  stochastic_age_struct_seir_ode(times,init,params)
}

# a little data post-processing
x_res_mult <- lapply(mult_sim_out ,FUN = function(x) as.data.frame(apply(x, 3, rowSums))) 
x_times <- as.numeric(str_remove(rownames(x_res_mult[[1]]), "[t]"))
x_res_mult1 <- x_res_mult %>%
  bind_rows(., .id = "sim") %>%
  mutate(time = rep(x_times,n_sims)) %>%
  select(sim, time, S:Rv_1d)
names(x_res_mult1) <- gsub("_", ".", names(x_res_mult1)) # need to do this for the pivot longer statement later

x_plot <- x_res_mult1 %>%
  select(-sim) %>%
  group_by(time) %>%
  summarise_all(list(mean = mean, 
                lower = ~quantile(., probs = 0.025),
                upper = ~quantile(., probs = 0.975))) %>%
  pivot_longer(!time,
               names_to = c("state", "stat"),
               names_sep = c("_"),
               values_to = "value") %>%
  pivot_wider(names_from = "stat",
              values_from = "value")

# plot
p <- ggplot(x_plot %>%
              filter(state %in% c("S", "E", "I", "R")), aes(x = time, y = mean, color = state, fill = state)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, colour = NA) + 
  ylab("Number of Individuals") +
  xlab("Time") +
  theme(legend.position = "bottom",
        panel.background = element_blank()#,
        #axis.text.x = element_text(angle = 45, hjust = 1)
        )
p
my_cols <- viridis_pal(option = "D")(24)
my_cols_transp <- paste0(my_cols, "1A")

par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(x_times, x_res_mult1, xlab = "Time", ylab = "Number of individuals",
        type = "l", col = c(rep(my_cols_transp, 10)), lty = 1)
legend("left", lwd = 1, col = my_cols, legend = ,c("S","E","I","R"), bty = "n")


