#' Summarise output from SEIR model for plotting
#' @param seir_output list of dataframes from post-processing
#' @param params list of parameter values
#' @param start_date calendar date of start of simulation
#' @param times vector of time points
#' @param vac_inputs vaccination inputs (from the output of convert_vac_schedule())
#' @return List of summary results
#' @keywords vacamole
#' @export

summarise_results <- function(seir_output, params, t_vec) {
  
  # get force of infection (lambda) --------------------------------------------
  calendar_day <- lubridate::yday(as.Date(t_vec, origin = params$calendar_start_date))
  beta_t <- params$beta * (1 + params$beta1 * cos(2 * pi * calendar_day / 365.24)) 
  lambda <- get_foi(x  = seir_output, 
                    y1 = params$eta_trans1[t_vec,-1], 
                    y2 = params$eta_trans2[t_vec,-1], 
                    y3 = params$eta_trans3[t_vec,-1], 
                    y4 = params$eta_trans4[t_vec,-1], 
                    y5 = params$eta_trans5[t_vec,-1],
                    beta = beta_t, 
                    contact_mat = params$contact_mat,
                    times = t_vec)
  
  # calculate infections -------------------------------------------------------
  new_infections <- lambda * (seir_output$S + seir_output$Shold_1d +
    (params$eta1[t_vec,-1] * (seir_output$Sv_1d + seir_output$Shold_2d)) +
    (params$eta2[t_vec,-1] * (seir_output$Sv_2d + seir_output$Shold_3d)) +
    (params$eta3[t_vec,-1] * (seir_output$Sv_3d + seir_output$Shold_4d)) +
    (params$eta4[t_vec,-1] * (seir_output$Sv_4d + seir_output$Shold_5d)) +
    (params$eta5[t_vec,-1] * seir_output$Sv_5d)
  ) %>%
    rename_with(., ~ paste0("age_group",1:9))
  
  new_infections <- new_infections %>%
    mutate(target_variable = "inc infection",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  # calculate cases ------------------------------------------------------------
  new_cases <- sweep(params$sigma * (seir_output$E + seir_output$Ev_1d + 
    seir_output$Ev_2d + seir_output$Ev_3d + seir_output$Ev_4d + seir_output$Ev_5d),
    2, params$p_report, "*") %>%
    rename_with(., ~ paste0("age_group",1:9)) 
  
  new_cases <- new_cases %>%
    mutate(target_variable = "inc case",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  
  # calculate hospital admissions ----------------------------------------------
  hosp_admissions <- sweep(seir_output$I + 
    params$eta_hosp1[t_vec,-1] * seir_output$Iv_1d +
    params$eta_hosp2[t_vec,-1] * seir_output$Iv_2d +
    params$eta_hosp3[t_vec,-1] * seir_output$Iv_3d +
    params$eta_hosp4[t_vec,-1] * seir_output$Iv_4d +
    params$eta_hosp5[t_vec,-1] * seir_output$Iv_5d,
    2, params$h, "*") %>%
    rename_with(., ~ paste0("age_group",1:9)) 
  
  hosp_admissions <- hosp_admissions %>%
    mutate(target_variable = "inc hosp",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  
  # calculate IC admissions ----------------------------------------------------
  ic_admissions <- sweep(seir_output$H + seir_output$Hv_1d + seir_output$Hv_2d + 
    seir_output$Hv_3d + seir_output$Hv_4d + seir_output$Hv_5d, 2, params$i1, "*") %>%
    rename_with(., ~ paste0("age_group",1:9)) 
  
  ic_admissions <- ic_admissions %>%
    mutate(target_variable = "inc icu",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  
  # calculate deaths -----------------------------------------------------------
  new_deaths <- 
    # from hospital compartment
    sweep(seir_output$H + seir_output$Hv_1d + seir_output$Hv_2d + 
          seir_output$Hv_3d + seir_output$Hv_4d + seir_output$Hv_5d, 2, params$d,
          "*") +
    # from IC compartment
    sweep(seir_output$IC + seir_output$ICv_1d + seir_output$ICv_2d + 
          seir_output$ICv_3d + seir_output$ICv_4d + seir_output$ICv_5d, 2, 
          params$d_ic, "*") + 
    # from hospital after IC
    sweep(seir_output$H_IC + seir_output$H_ICv_1d + seir_output$H_ICv_2d + 
         seir_output$H_ICv_3d + seir_output$H_ICv_4d + seir_output$H_ICv_5d, 2,
         params$d_hic, "*")
   
  
  new_deaths <- new_deaths %>%
    rename_with(., ~ paste0("age_group",1:9)) %>%
    mutate(target_variable = "inc death",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  
  # Create object into format for scenario hub ---------------------------------
  # bind all outcome data frames together and wrangle
  rtn <- bind_rows(new_infections, new_cases, hosp_admissions, ic_admissions, 
                   new_deaths) %>%
    pivot_longer(cols = age_group1:age_group9,
                 names_to = "age_group",
                 names_prefix = "age_group",
                 values_to = "value") %>%
    select(time, date, epiweek, year, horizon, target_variable, age_group, value)

  return(rtn)
}
