# COVID scenarios project
# Use different Pico contact matrices to simulate different measures in first wave

# Set-up -------------------------------------------
# --------------------------------------------------
# load packages
devtools::load_all()
library(vacamole)

source("inst/extdata/scripts/helpers/model_run_helper.R")

# Read in OSIRIS data ------------------------------
last_date_in_osiris <- "2021-11-28"
case_data <- readRDS(paste0("inst/extdata/data/case_data_upto_", last_date_in_osiris,".rds")) 
first_wave_data <- case_data %>%
  filter(date <= as.Date("2020-07-01"))

# plot data for first wave
p <- ggplot(first_wave_data, aes(x = date, y = inc)) +
  geom_line() +
  geom_line(aes(x = date, y = roll_avg, color = "red")) +
  theme(panel.background = element_blank())
p

# specify model parameters --------------------------
params <- list(beta = 0.0003934816,       # transmission rate
               beta1 = 0.14,              # amplitude of seasonal forcing
               gamma = g,                 # 1/gamma = infectious period
               sigma = s,                 # 1/sigma = latent period
               epsilon = 0.01,            # import case rate
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
               vac_inputs = NULL,
               use_cases = TRUE,                           # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 35.7/100000 * sum(n_vec),        # 35.7  # 20 for IC admissions
               no_vac = TRUE,
               t_calendar_start = yday(as.Date("2020-01-01")),   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
               beta_change = NULL
)
# ---------------------------------------------------
# define wrapper function ---------------------------
scenarios_wrapper <- function(contact_matrix, data, beta_mles, beta_draws, params){

# split list of mles
beta_mles_list <- split(beta_mles, seq(nrow(beta_mles)))
  
# specify which contact matrix to use after exponential growth
if(contact_matrix == "april_2017"){cm <- april_2017
} else if (contact_matrix == "april_2020"){cm <- april_2020
} else if (contact_matrix == "june_2020"){cm <- june_2020
} else if (contact_matrix == "september_2020"){cm <- september_2020
} else if (contact_matrix == "february_2021"){cm <- february_2021
} else {cm <- june_2021}

# define time window breakpoints
breakpoints <- list(
  date = c(as.Date("2020-03-16"),
           as.Date("2020-03-24"),
           as.Date("2020-04-29"),
           as.Date("2020-05-11"),
           as.Date("2020-06-01"),
           as.Date("2020-07-01")),
  contact_matrix = list(april_2017, cm, cm, cm, cm, cm))

# run for 200 contact matrices
n_bp <- length(breakpoints$date)

mle_run <- model_run_wrapper(breakpoints = breakpoints, beta_values = beta_mles_list, init_conditions = init, params = params)
ci_run  <- model_run_wrapper(breakpoints = breakpoints, beta_values = beta_draws, init_conditions = init, params = params, mle = FALSE)
ci_out <- list()
for (i in 1:n_bp){
  ci_out[[i]] <- do.call("rbind", ci_run[[i]])
}

ci_out_wide <- do.call("cbind", ci_out)
#matplot(t(ci_out_wide), type = "l")

# get quantiles and plot
bounds <- apply(ci_out_wide, 2, quantile, probs = c(0.025, 0.975))
#matplot(t(bounds), type = "l")

# save outputs -------------------------------------
# --------------------------------------------------
#  combine all piecewise results to plot together
cases_mle <- unique(unlist(mle_run))
cases_lower <- unique(bounds[1,])
cases_upper <- unique(bounds[2,])
times_all <- 1:length(cases_mle)

model_fit <- data.frame(time = times_all, date = first_wave_data$date, real = first_wave_data$inc, mle = cases_mle, lower = cases_lower, upper = cases_upper)
saveRDS(model_fit, file = paste0(path_out, "scenario_fit_df_", contact_matrix, ".rds"))
# --------------------------------------------------
return(model_fit)
}

# ----------------------------------------------------
# run simulations for mle, lower, and upper bounds 
# of beta
# ----------------------------------------------------
fit_date <- "2021-10-01"
path_out <- "/rivm/s/ainsliek/code/vacamole/inst/extdata/results/model_fits/"

# read in beta estimates
b_mles <- readRDS(paste0(path_out,"mles_from_fits_",fit_date,".rds"))
b_draws <- readRDS(paste0(path_out,"beta_draws_from_fits_",fit_date,".rds"))

# run scenarios wrapper function
res_april_2017     <- scenarios_wrapper(contact_matrix = "april_2017", data = first_wave_data, beta_mles = b_mles, beta_draws = b_draws, params = params)
res_april_2020     <- scenarios_wrapper(contact_matrix = "april_2020", data = first_wave_data, beta_mles = b_mles, beta_draws = b_draws, params = params)
res_june_2020      <- scenarios_wrapper(contact_matrix = "june_2020", data = first_wave_data, beta_mles = b_mles, beta_draws = b_draws, params = params)
res_september_2020 <- scenarios_wrapper(contact_matrix = "september_2020", data = first_wave_data, beta_mles = b_mles, beta_draws = b_draws, params = params)
res_february_2021  <- scenarios_wrapper(contact_matrix = "february_2021", data = first_wave_data, beta_mles = b_mles, beta_draws = b_draws, params = params)
res_june_2021      <- scenarios_wrapper(contact_matrix = "june_2021", data = first_wave_data, beta_mles = b_mles, beta_draws = b_draws, params = params)

# --------------------------------------------------
# combine all results together and plot together
# uncomment to read in saved model fits, so they don't have to be run again
# res_april_2017     <- readRDS(file = paste0(path_out, "scenario_fit_df_april_2017.rds"))
# res_april_2020     <- readRDS(file = paste0(path_out, "scenario_fit_df_april_2020.rds"))
# res_june_2020      <- readRDS(file = paste0(path_out, "scenario_fit_df_june_2020.rds"))
# res_september_2020 <- readRDS(file = paste0(path_out, "scenario_fit_df_september_2020.rds"))
# res_february_2021  <- readRDS(file = paste0(path_out, "scenario_fit_df_february_2021.rds"))
# res_june_2021      <- readRDS(file = paste0(path_out, "scenario_fit_df_june_2021.rds"))

results <- bind_rows(res_april_2017, res_april_2020, res_june_2020, res_september_2020, res_february_2021, res_june_2021,
                     .id = "contact_matrix") %>%
  mutate(scenario = case_when(
    contact_matrix == 1 ~ "april_2017",
    contact_matrix == 2 ~ "april_2020",
    contact_matrix == 3 ~ "june_2020",
    contact_matrix == 4 ~ "september_2020",
    contact_matrix == 5 ~ "february_2021",
    contact_matrix == 6 ~ "june_2021"
  ))

# linear scale ------------------------------------
p_lin <- ggplot(data = results, #%>%
               #filter(scenario == "april_2017"), 
             aes(x = date, y = mle, color = scenario)) + # linetype="solid", 
  #geom_point(data = results, aes(x = date, y = real, color = "Osiris notifications")) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = scenario), alpha = 0.05, colour = NA) +
  #scale_color_manual(values = c("red"),
  #                   labels = c("Osiris notifications")) +
  #scale_fill_manual(values = c("grey70")) +
  #scale_linetype_manual(values=c(1), labels = c("Model Fit")) +
  #scale_shape_manual(values=c(NA,20)) +
  labs(y = "Daily Cases", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom", #"none"
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title=element_text(size=14)) #+
  #facet_wrap(.~scenario, scales = "free_y")
p_lin

# log scale ----------------------------------------
p_log <- ggplot(data = results, #%>%
                #filter(scenario == "april_2017"), 
                aes(x = date, y = mle, color = scenario)) + # linetype="solid", 
  #geom_point(data = results, aes(x = date, y = real, color = "Osiris notifications")) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = scenario), alpha = 0.05, colour = NA) +
  scale_y_continuous(trans='log10') + 
  #scale_color_manual(values = c("red"),
  #                   labels = c("Osiris notifications")) +
  #scale_fill_manual(values = c("grey70")) +
  #scale_linetype_manual(values=c(1), labels = c("Model Fit")) +
  #scale_shape_manual(values=c(NA,20)) +
  labs(y = "LogDaily Cases", x = "Date") +
  #ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom", #"none"
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title=element_text(size=14)) #+
#facet_wrap(.~scenario, scales = "free_y")
p_log

# save figures -------------------------------------
ggsave(filename = "/rivm/s/ainsliek/results/covid_scenarios/pico_scenarios_linear_plot_facet_wrap.pdf", plot = p_lin,
       units = "in", height = 8, width = 10, dpi = 300)

ggsave(filename = "/rivm/s/ainsliek/results/covid_scenarios/pico_scenarios_log_plot.pdf", plot = p_log,
       units = "in", height = 8, width = 10, dpi = 300)

# --------------------------------------------------
# old code -----------------------------------------
# --------------------------------------------------
# initial values
# empty_state <- c(rep(0, 9))
# init <- c(
#   t = 0,
#   S = c(n_vec[1:4], n_vec[5]-1, n_vec[6:9]),
#   Shold_1d = empty_state,
#   Sv_1d = empty_state,
#   Shold_2d = empty_state,
#   Sv_2d = empty_state,
#   E = empty_state,
#   Ev_1d = empty_state,
#   Ev_2d = empty_state,
#   I = c(rep(0,4),1,rep(0,4)),
#   Iv_1d = empty_state,
#   Iv_2d = empty_state,
#   H = empty_state,
#   Hv_1d = empty_state,
#   Hv_2d = empty_state,
#   H_IC = empty_state,
#   H_ICv_1d = empty_state,
#   H_ICv_2d = empty_state,
#   IC = empty_state,
#   ICv_1d = empty_state,
#   ICv_2d = empty_state,
#   D = empty_state,
#   R = empty_state,
#   Rv_1d = empty_state,
#   Rv_2d = empty_state
# )
# 

# --------------------------------------------------
# fit model to data for exponential growth period
# --------------------------------------------------
# first we need to determine transmission rates
# during the first wave
# 1 January 2020 to 15 March 2020
# 15 March 2020 to 1 July 2020

# define time window breakpoints
# breakpoints <- list( 
#   date = c(as.Date("2020-03-16"), 
#            as.Date("2020-03-24"),
#            as.Date("2020-04-29"),
#            as.Date("2020-05-11"),
#            as.Date("2020-06-01"),
#            as.Date("2020-07-01")),  
#   contact_matrix = list(april_2017, cm))
# 
# n_bp <- length(breakpoints$date)
# 
# # create empty objects to store outputs
# mles <- matrix(rep(NA, 2*n_bp), nrow = n_bp)
# colnames(mles) <- c("beta", "alpha")
# out_mle <- list()
# parameter_draws <- list()
# beta_draws <- list()
# daily_cases_mle <- list()
# 
# # loop over time windows
# for (j in 1:n_bp) {
#   # set contact matrix for time window
#   if (j == n_bp){
#     params$c_start <- breakpoints$contact_matrix[[j-1]]$mean
#   } else {
#     params$c_start <- breakpoints$contact_matrix[[j]]$mean
#   }
#   
#   if (j == 1) {
#     # if first time window, start time at 0
#     end_day <- yday(breakpoints$date[j]) - 1
#     
#     times <- seq(0, end_day, by = 1)
#     # set initial conditions to those specified earlier in script
#     init_update <- init
#     pars <- c(2.3, 0.01)
#     S_diag <- diag(init_update[c(2:10)])
#     rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
#   } else {
#     start_day <- yday(breakpoints$date[j-1]) - 1 
#     end_day <- yday(breakpoints$date[j]) - 1
#     times <- seq(start_day, end_day, by = 1)
#     
#     # update initial conditions based on last time window
#     init_update <- c(t = times[1], unlist(lapply(unname(out_mle[[j-1]]), tail,1)))
#     
#     pars <- c((mles[j-1,1]/params$gamma)*rho, mles[j-1,2])
#     S_diag <- diag(init_update[c(2:10)])
#     rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
#   }
#   
#   # subset data for time window
#   case_data_sub <- case_data[times + 1, ]
#   
#   # optimize
#   res <- optim(par = pars, 
#                fn = likelihood_func,
#                method = "L-BFGS-B",
#                lower = c(0,0.005),
#                upper = c(10,1),
#                t = times,
#                data = case_data_sub,
#                params = params,
#                init = init_update,
#                stochastic = FALSE,
#                hessian = TRUE
#   )
#   
#   # store MLEs
#   mles[j,1] <- (res$par[1] / rho) * params$gamma
#   mles[j,2] <- res$par[2]
#   #mles[j,3] <- res$par[3]
#   
#   # draw 200 parameter values
#   hess_mat <- res$hessian
#   parameter_draws[[j]] <- mvtnorm::rmvnorm(200, res$par, solve(hess_mat))
#   beta_draws[[j]] <- data.frame(beta = (parameter_draws[[j]][,1] / rho) * params$gamma) %>%
#     mutate(index = 1:200)
# # --------------------------------------------------
# # run for mle to get initial conditions for next timepoint
#   params$beta <- mles[j,1]
#   seir_out <- lsoda(init_update, times, age_struct_seir_ode, params)
#   seir_out <- as.data.frame(seir_out)
#   out_mle[[j]] <- postprocess_age_struct_model_output(seir_out)
#   cases <- params$sigma * rowSums(out_mle[[j]]$E + out_mle[[j]]$Ev_1d + out_mle[[j]]$Ev_2d) * params$p_report
#   
#   # plot for quick check of fit
#   plot(cases~times, type = "l")
#   points(times, case_data_sub$inc, pch = 16, col = "red")
# } # end loop
# 
# # ---------------------------------------------------
# # output
# todays_date <- Sys.Date()
# path_out <- "/rivm/s/ainsliek/code/vacamole/inst/extdata/results/model_fits/"
# # save outputs
# saveRDS(mles, file = paste0(path_out, "mles_exp_growth_", todays_date, ".rds"))
# saveRDS(beta_draws, file = paste0(path_out, "beta_draws_exp_growth_", todays_date, ".rds"))
# saveRDS(out_mle, file = paste0(path_out, "output_from_exp_growth_", todays_date, ".rds"))
# 


