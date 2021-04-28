# Simulate age-structured SEIR compartmental model with 2 vaccine doses

# Data and model parameters are loaded/defined in the script inst/extdata/scripts/model_run_helper.R

# Create list of parameter values for input into model solver
params <- list(beta = beta2_prime,           # transmission rate
               gamma = g,                      # 1/gamma = infectious period
               sigma = s,                      # 1/sigma = latent period
               delta = 1,                      # scaling constant for beta
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
               c_start = B_prime,
               c_lockdown = B_prime,
               c_relaxed = t4,
               c_very_relaxed = t3,
               c_normal = t1,
               vac_inputs = basis1,
               use_cases = TRUE,              # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 35.7/100000 * sum(n_vec),        # 20 for IC admissions
               no_vac = FALSE,
               t_calendar_start = yday(as.Date("2021-01-31")),   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
               breakpoints = c() # time points when parameters can change (if NULL, then beta is constant over time)
)

# Determine time points over which to solve model
t_max <- lubridate::yday(as.Date("2021-04-20")) # last day of osiris data that we fit
times_fit <- seq(0,t_max - params$t_calendar_start, by = 1)  


end_date <- lubridate::yday(as.Date("2021-08-31"))
times_forward <- seq(0, end_date - params$t_calandar_start, by = 1)
# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times_forward,age_struct_seir_ode,params) #
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)

# quick check ------------------------------------------------------
cases_from_fit <- (params$sigma * (init_from_fit1$E + init_from_fit1$Ev_1d + init_from_fit1$Ev_2d)) / 3
cases <- (params$sigma * (out$E + out$Ev_1d + out$Ev_2d)) / 3
cases_all <- bind_rows(cases_from_fit, cases)
plot(seq(1, dim(cases_all)[1], by = 1), rowSums(cases_all), type = "l", col = "blue",
     xlab="Time (days)",ylab="Daily Cases", ylim = c(0,7000))
points(osiris2$inc~seq(1,weeks*7,by=1),col="red",pch=16) 


# Summarise results ------------------------------------------------
tag <- "basis_19April"
results <- summarise_results(out, params, start_date = "2021-04-04", 
                             times = times_forward, vac_inputs = basis1)
saveRDS(results, paste0("inst/extdata/results/res_",tag,".rds"))

# Make summary table ------------------------------------------------
summary_tab <- results$df_summary %>%
  filter(date >= as.Date("2021-04-01"),
         date < as.Date("2021-09-01")) %>%
  group_by(outcome) %>%
  summarise_at(.vars = "value", .funs = sum) 
summary_tab

# Make plot ---------------------------------------------------------
# summary over all age groups
p <- ggplot(results$df_summary %>%
              filter(outcome == "new_infections"), aes(x = date, y = value, color = outcome)) +
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

# by age group
p2 <- ggplot(results$df, aes(x = date, y = value, color = age_group)) +
  geom_line() +
  labs(y = "Value", x = "Time (days)", color = "Outcome") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~outcome, scales = "free")
p2

# plot vaccination
ggplot(results$df_vac, aes(x = date, y = value, color = dose)) +
  geom_line() +
  labs(y = "Value", x = "Time (days)", color = "Dose") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  #geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) #+
  #facet_wrap(~outcome, scales = "free")

