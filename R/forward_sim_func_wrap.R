# ------------------------------------------------------------------
# Helper data wrangling functions ----------------------------------
# ------------------------------------------------------------------

wide_to_long <- function(x, times) {
  
  x %>% 
    mutate(time = times) %>%
    pivot_longer(cols = !time,
                 names_to = c("state", "age_group"),
                 names_sep = -1,
                 values_to = "value")
}

wrangle_results <- function(x){
  
  rtn_mle <- x %>%
    filter(sim == 0)
  
  rtn_bounds <- x %>%
    filter(sim != 0) %>%
    group_by(time, state, age_group) %>%
    summarise(lower = quantile(value, probs = 0.025),
              upper = quantile(value, probs = 0.975)) 
  
  rtn <- left_join(rtn_mle, rtn_bounds, by = c("time", "state", "age_group")) %>%
    rename(mean = value) %>%
    select(-sim)
  
  return(rtn)
}

# ------------------------------------------------------------------
# Forward simulation function wrapper ------------------------------
# ------------------------------------------------------------------

forward_sim_func_wrap <- function(start_date,
                                  end_date,
                                  init_cond,
                                  beta_m,
                                  vac_inputs,
                                  beta_c,
                                  beta_d,
                                  t_normal,
                                  contact_matrices,
                                  tag){
  
  # empty list for output 
  out <- list()
  # specify time points ----------------------------------------------
  start_date <- lubridate::yday(as.Date(start_date) + 365) + 365
  end_date <- lubridate::yday(as.Date(end_date)) + (365*2) 
  times <- seq(start_date, end_date, by = 1)

  initial_conditions <- c(t=times[1], init_cond)

  # Create list of parameter values for input into model solver ------
  params <- list(beta = beta_m,             # transmission rate
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
                 c_start = contact_matrices$june_2021$mean,
                 c_lockdown = contact_matrices$february_2021$mean,
                 c_relaxed = contact_matrices$june_2020$mean,
                 c_very_relaxed = contact_matrices$june_2021$mean,
                 c_normal = contact_matrices$baseline_2017$mean,
                 keep_cm_fixed = FALSE,
                 vac_inputs = vac_inputs,
                 use_cases = TRUE,                           # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
                 thresh_n = 0.5/100000 * sum(n_vec),         # somewhat arbitrary cut-off ***need to check if realistic
                 thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
                 thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
                 thresh_u = 100000/100000 * sum(n_vec),        #35.7  # 20 for IC admissions
                 no_vac = FALSE,
                 t_calendar_start = yday(as.Date("2020-01-01")),   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
                 beta_change =  beta_c[1],
                 t_normal = t_normal
  )


  # if time doesn't start at 0 we need to initialise the contact 
  # matrices flags ---------------------------------------------------
  assign("flag_relaxed", 0, envir = .GlobalEnv) 
  assign("flag_very_relaxed", 0, envir = .GlobalEnv)
  assign("flag_normal", 0, envir = .GlobalEnv)
  
  #  Solve model ------------------------------------------------------
  # mle
  print("mle")
  seir_out <- lsoda(initial_conditions,times,age_struct_seir_ode,params) 
  seir_out <- as.data.frame(seir_out)
  out_mle <- postprocess_age_struct_model_output(seir_out)
  tmp_mle <- lapply(out_mle, wide_to_long, times) %>%
    bind_rows() %>%
    mutate(sim = 0)

  # run for beta draws ------------------------------------------------
  n_beta <- length(beta_d)
  rtn_out <- list()

  for(i in 1:n_beta){
    print(i)
  
    # reset flags
    # need to make sure they are saved to the global environment
    assign("flag_relaxed", 0, envir = .GlobalEnv) 
    assign("flag_very_relaxed", 0, envir = .GlobalEnv)
    assign("flag_normal", 0, envir = .GlobalEnv)
    
    # change parameters
    params$beta <- beta_d[i]
    params$beta_c <- beta_c[i+1]
    params$c_start <- contact_matrices$june_2021[[i]]
    params$normal <- contact_matrices$baseline_2017[[i]]
    
    # run model
    seir_out <- lsoda(initial_conditions, times, age_struct_seir_ode, params)
    seir_out <- as.data.frame(seir_out)
    out <- postprocess_age_struct_model_output(seir_out)
    
    # convert output from each simulation into long data frame
    tmp <- lapply(out, wide_to_long, times) %>%
        bind_rows() %>%
        mutate(sim = i)
  
    rtn_out[[i]] <- tmp
  }

  
  # return outouts
  rtn_out$mle <- tmp_mle
  rtn <- bind_rows(rtn_out)
  
  saveRDS(rtn, paste0("inst/extdata/results/",tag,".rds"))
  
  return(rtn)

} # end of function

# test plot
# ggplot(data = test %>%
#          filter(state == "E",
#                 time < 600), aes(x = time, y = value, color = as.factor(sim))) +
#   geom_line() +
#   facet_wrap(~age_group)


