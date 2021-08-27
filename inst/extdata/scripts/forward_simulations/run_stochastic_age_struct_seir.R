# Simulate age-structured SEIR compartmental model with 2 vaccine doses

# Data and model parameters are loaded/defined in the script inst/extdata/scripts/model_run_helper.R
source("inst/extdata/scripts/model_run_helper.R")

beta_mle <- 0.00061
start_date <- lubridate::yday(as.Date("2021-05-25"))
end_date <- lubridate::yday(as.Date("2022-03-30")) + 365
# Create list of parameter values for input into model solver
params <- list(dt = 1/6,
               beta = beta_mle,           # transmission rate
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
               c_start = t2,
               c_lockdown = t3,
               c_relaxed = t4,
               c_very_relaxed = t3,
               c_normal = t1,
               keep_cm_fixed = FALSE,
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

times <- seq(start_date, end_date, by = 1) - yday(as.Date("2021-01-31"))
# set initial_conditions start time the first time point
initial_conditions["t"] <- times[1]
# single simulation run --------------------------------------------
# if time doesn't start at 0 we need to initialise the contact 
# matrices flags
if(times[1] != 0){
  flag_relaxed <- 0 # start with relaxed contact matrix
  flag_very_relaxed <- 0
  flag_normal <- 0
}
# Solve model ------------------------------------------------------
seir_out <- lsoda(initial_conditions,times,stochastic_age_struct_seir_ode,params) #
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)

# parallel simulations ---------------------------------------------
# simulations to run
sims_vec <- 1:10

# set up parallel
library(parallel)
n_cores <- detectCores()

# wrapper function for parallel runs
# can't use this on Windows!
parallel_wrapper <- function(n_sim){
  
  # if time doesn't start at 0 we need to initialise the contact matrices flags
  if(times[1] != 0){
    flag_relaxed <- 0 # start with relaxed contact matrix
    flag_very_relaxed <- 0
    flag_normal <- 0
  }
  # Solve model ----------------------------------------------------
  seir_out <- lsoda(initial_conditions,times,stochastic_age_struct_seir_ode,params) #
  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output(seir_out)
  return(out)
}

out_parallel <- mclapply(sims_vec, parallel_wrapper, mc.cores = n_cores)

# quick check ------------------------------------------------------
# cases_from_fit <- (params$sigma * (init_from_fit1$E + init_from_fit1$Ev_1d + init_from_fit1$Ev_2d)) / 3
cases <- (params$sigma * (out$E + out$Ev_1d + out$Ev_2d)) * params$p_report
plot(seq(1, dim(cases)[1], by = 1), rowSums(cases), type = "l", col = "blue",
     xlab="Time (days)",ylab="Daily Cases")
points(osiris2$inc~times,col="red",pch=16) 
#lines(seq(1, dim(cases)[1], by = 1), rowSums(cases), col = "black")

# Summarise results ------------------------------------------------
tag <- "basis_no_wane_beta1_0_15_16june"
results <- summarise_results(out, params, start_date = "2021-01-31", 
                             times = times, vac_inputs = params$vac_inputs)
saveRDS(results, paste0("inst/extdata/results/",tag,".rds"))

# Make summary table ------------------------------------------------
summary_tab <- results$df_summary %>%
  # filter(date >= as.Date("2021-04-01"),
  #        date < as.Date("2021-09-01")) %>%
  group_by(outcome) %>%
  summarise_at(.vars = "value", .funs = sum) 
summary_tab

# Make plot ---------------------------------------------------------
# summary over all age groups
p <- ggplot(results$df_summary %>%
              filter(outcome == "new_cases"), aes(x = date, y = value, color = outcome)) +
  geom_line() +
  labs(y = "Value", x = "Time (days)", color = "Outcome") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  #geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~outcome, scales = "free")
p

# combine with model fit
basis_cases <- results$df %>%
  filter(outcome == "new_cases") %>%
  select(time, age_group, outcome, value) %>%
  group_by(time) %>%
  summarise_at(.vars = "value", .funs = "sum") %>%
  filter(time != 79) %>% # remove first time point because it's a repeat of the last time point of the data fit
  rename(cases = value)

cases_fit_and_model <- bind_rows(daily_cases_from_fit, basis_cases) %>%
  mutate(date = time + as.Date("2021-01-31"))

p_fit <- ggplot(cases_fit_and_model, aes(x = date, y = cases)) +
  geom_line() +
  geom_point(data = osiris1 %>% filter(date < as.Date("2021-05-26")), aes(x = date, y = inc, color = "Osiris data")) +
  labs(y = "Daily Cases", x = "Time (days)") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  #geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_fit

# save plot to file
ggsave(paste0("inst/extdata/results/plot_",tag,".jpg"),
       plot = p_fit,
       height = 6,
       width = 12,
       dpi = 300)


