# ------------------------------------------------------------------
# Hleper data wrangling function -----------------------------------
# ------------------------------------------------------------------

wide_to_long <- function(x, times) {
  
  x %>% 
    mutate(time = times) %>%
    pivot_longer(cols = !time,
                 names_to = c("state", "age_group"),
                 names_sep = -1,
                 values_to = "value")
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
                 beta_change =  beta_c 
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
  
  # get outcomes ------------------------------------------------------
  # cases_mle <- params$sigma * rowSums(out_mle$E + out_mle$Ev_1d + out_mle$Ev_2d) * params$p_report
  # infectious_mle <- (out_mle$I + out_mle$Iv_1d + out_mle$Iv_2d)
  # hosp_admissions_mle <- rowSums(sweep(infectious_mle, 2, h, "*"))
  # hosp_occ_mle <- (out_mle$H + out_mle$Hv_1d + out_mle$Hv_2d)
  # ic_mle <- rowSums(sweep(hosp_occ_mle, 2, i1, "*"))
  # ic_occ_mle <- (out_mle$IC + out_mle$ICv_1d + out_mle$ICv_2d)
  # hosp_after_ic_mle <- sweep(ic_occ_mle, 2, i2, "*")
  # hosp_after_ic_occ_mle <- (out_mle$H_IC + out_mle$H_ICv_1d + out_mle$H_ICv_2d)
  # deaths_mle <- rowSums(sweep(ic_occ_mle, 2, d_ic, "*") + sweep(hosp_occ_mle, 2, d, "*") + sweep(hosp_after_ic_occ_mle, 2, d_hic, "*"))


  # run for beta draws ------------------------------------------------
  n_beta <- length(beta_d)
  # rtn_cases <- matrix(,nrow = length(times), ncol = n_beta)
  # rtn_hosp <- rtn_cases
  # rtn_ic <- rtn_cases
  # rtn_deaths <- rtn_cases

  rtn_out <- list()

  for(i in 1:n_beta){
    print(i)
  
    # reset flags
    # flag_relaxed <- 0
    # flag_very_relaxed <- 0
    # flag_normal <- 0
    
    assign("flag_relaxed", 0, envir = .GlobalEnv) 
    assign("flag_very_relaxed", 0, envir = .GlobalEnv)
    assign("flag_normal", 0, envir = .GlobalEnv)
    
    # change parameters
    params$beta <- beta_d[i]
    params$c_start <- contact_matrices$june_2021[[i]]
    params$normal <- contact_matrices$baseline_2017[[i]]
    
    # print(params$beta)
    # print(params$c_start)
    # print(params$normal)
    # run model
    seir_out <- lsoda(initial_conditions, times, age_struct_seir_ode, params)
    seir_out <- as.data.frame(seir_out)
    out <- postprocess_age_struct_model_output(seir_out)
    
    # convert output from each simulation into long data frame
    tmp <- lapply(out, wide_to_long, times) %>%
        bind_rows() %>%
        mutate(sim = i)
  
    rtn_out[[i]] <- tmp
  
    # get outcomes
    # daily_cases <- params$sigma * rowSums(out$E + out$Ev_1d + out$Ev_2d) * params$p_report
    # infectious <- (out$I + out$Iv_1d + out$Iv_2d)
    # hosp_admissions <- rowSums(sweep(infectious, 2, h, "*"))
    # hosp_occ <- (out$H + out$Hv_1d + out$Hv_2d)
    # ic <- rowSums(sweep(hosp_occ, 2, i1, "*"))
    # ic_occ <- (out$IC + out$ICv_1d + out$ICv_2d)
    # hosp_after_ic <- sweep(ic_occ, 2, i2, "*")
    # hosp_after_ic_occ <- (out$H_IC + out$H_ICv_1d + out$H_ICv_2d)
    # daily_deaths <- rowSums(sweep(ic_occ, 2, d_ic, "*") + sweep(hosp_occ, 2, d, "*") + sweep(hosp_after_ic_occ, 2, d_hic, "*"))
    # 
    # rtn_cases[,i]  <- daily_cases
    # rtn_hosp[,i]   <- hosp_admissions
    # rtn_ic[,i]     <- ic
    # rtn_deaths[,i] <- daily_deaths
  }

  # get confidence bounds of model runs for total pop
  # bounds_cases  <- apply(rtn_cases,1,function(x) quantile(x, c(0.025,0.975)))
  # bounds_hosp   <- apply(rtn_hosp,1,function(x) quantile(x, c(0.025,0.975)))
  # bounds_ic     <- apply(rtn_ic,1,function(x) quantile(x, c(0.025,0.975)))
  # bounds_deaths <- apply(rtn_deaths,1,function(x) quantile(x, c(0.025,0.975)))
  # 
  # # output ---------------------------------------------------------
  # df_total <- data.frame(time = times,
  #                        cases_mle = cases_mle,
  #                        cases_lower = bounds_cases[1,],
  #                        cases_upper = bounds_cases[2,],
  #                        hosp_mle = hosp_admissions_mle,
  #                        hosp_lower = bounds_hosp[1,],
  #                        hosp_upper = bounds_hosp[2,],
  #                        ic_mle = ic_mle,
  #                        ic_lower = bounds_ic[1,],
  #                        ic_upper = bounds_ic[2,],
  #                        deaths_mle = deaths_mle,
  #                        deaths_lower = bounds_deaths[1,],
  #                        deaths_upper = bounds_deaths[2,]
  # )
  
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


