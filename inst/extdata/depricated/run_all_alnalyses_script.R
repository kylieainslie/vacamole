# main analysis
for (i in 1:4){
  scenario_names <- c("o2y", "y2o", "alt", "no_vac_healthy")
  if (i == 1){ vac_strat <- prop_o2y}
  if (i == 2){ vac_strat <- prop_y2o}
  if (i == 3){ vac_strat <- prop_alt}
  if (i == 4){ vac_strat <- prop_no_vac_healthy}
  
  t_max <- dim(vac_strat)[1] - 1
  times <- seq(0,t_max, by = 1)     # Vector of times 242 = Sept 30, 2021
  
  
params <- list(beta = beta2_prime,           # transmission rate
               gamma = g,                      # 1/gamma = infectious period
               sigma = s,                      # 1/sigma = latent period
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
               c_start = B_prime,
               c_lockdown = B_prime,
               c_relaxed = t4,
               c_very_relaxed = t3,
               c_normal = t1,
               vac_schedule = vac_strat,
               ve = ve,
               delay = delays,
               use_cases = TRUE,              # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 35.7/100000 * sum(n_vec),        # 20 for IC admissions
               no_vac = FALSE,
               #force_relax = NULL,                          # time step when measures are forced to relax regardless of criteria
               t_start_end = 0#,                           # time step when starting contact matrix ends and criteria are used to decide contact matrix
               #init_lambda = beta2_prime * B_prime %*% init_states_dat$init_I
)

# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times,age_struct_seir_ode,params) #
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)

# Summarise results ------------------------------------------------
tag <- paste0(scenario_names[i],"_15March")
results <- summarise_results(out, params, start_date = "2021-01-31", times)
saveRDS(results, paste0("inst/extdata/results/res_",tag,".rds"))
}


# sensitivity analyses
# sa20
for (i in 1:4){
  scenario_names <- c("o2y", "y2o", "alt", "no_vac_healthy")
  if (i == 1){ vac_strat <- prop_o2y}
  if (i == 2){ vac_strat <- prop_y2o}
  if (i == 3){ vac_strat <- prop_alt}
  if (i == 4){ vac_strat <- prop_no_vac_healthy}
  
  params <- list(beta = beta2_prime,           # transmission rate
                 gamma = g,                      # 1/gamma = infectious period
                 sigma = s,                      # 1/sigma = latent period
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
                 c_start = B_prime,
                 c_lockdown = B_prime,
                 c_relaxed = t4,
                 c_very_relaxed = t3,
                 c_normal = t1,
                 vac_schedule = vac_strat,
                 ve = ve_sa20,
                 delay = delays,
                 use_cases = TRUE,              # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
                 thresh_n = 0.5/100000 * sum(n_vec),
                 thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
                 thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
                 thresh_u = 35.7/100000 * sum(n_vec),        # 20 for IC admissions
                 no_vac = FALSE,
                 #force_relax = NULL,                          # time step when measures are forced to relax regardless of criteria
                 t_start_end = 0#,                           # time step when starting contact matrix ends and criteria are used to decide contact matrix
                 #init_lambda = beta2_prime * B_prime %*% init_states_dat$init_I
  )
  
  # Solve model ------------------------------------------------------
  seir_out <- lsoda(init,times,age_struct_seir_ode,params) #
  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output(seir_out)
  
  # Summarise results ------------------------------------------------
  tag <- paste0(scenario_names[i],"_sa20_15March")
  results <- summarise_results(out, params, start_date = "2021-01-31", times)
  saveRDS(results, paste0("inst/extdata/results/res_",tag,".rds"))
}

# sa50
for (i in 1:4){
  scenario_names <- c("o2y", "y2o", "alt", "no_vac_healthy")
  if (i == 1){ vac_strat <- prop_o2y}
  if (i == 2){ vac_strat <- prop_y2o}
  if (i == 3){ vac_strat <- prop_alt}
  if (i == 4){ vac_strat <- prop_no_vac_healthy}
  
  params <- list(beta = beta2_prime,           # transmission rate
                 gamma = g,                      # 1/gamma = infectious period
                 sigma = s,                      # 1/sigma = latent period
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
                 c_start = B_prime,
                 c_lockdown = B_prime,
                 c_relaxed = t4,
                 c_very_relaxed = t3,
                 c_normal = t1,
                 vac_schedule = vac_strat,
                 ve = ve_sa50,
                 delay = delays,
                 use_cases = TRUE,              # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
                 thresh_n = 0.5/100000 * sum(n_vec),
                 thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
                 thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
                 thresh_u = 35.7/100000 * sum(n_vec),        # 20 for IC admissions
                 no_vac = FALSE,
                 #force_relax = NULL,                          # time step when measures are forced to relax regardless of criteria
                 t_start_end = 0#,                           # time step when starting contact matrix ends and criteria are used to decide contact matrix
                 #init_lambda = beta2_prime * B_prime %*% init_states_dat$init_I
  )
  
  # Solve model ------------------------------------------------------
  seir_out <- lsoda(init,times,age_struct_seir_ode,params) #
  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output(seir_out)
  
  # Summarise results ------------------------------------------------
  tag <- paste0(scenario_names[i],"_sa50_15March")
  results <- summarise_results(out, params, start_date = "2021-01-31", times)
  saveRDS(results, paste0("inst/extdata/results/res_",tag,".rds"))
}

