#' Summarise output from SEIR model for plotting
#' @param seir_output list of dataframes from post-processing
#' @param params list of parameter values
#' @param start_date calendar date of start of simulation
#' @param times vector of time points
#' @param vac_inputs vaccination inputs (from the output of convert_vac_schedule())
#' @return List of summary results
#' @keywords vacamole
#' @export

summarise_results <- function(seir_output, params) {
  
  # get time vector ------------------------------------------------------------
  times <- seir_out$time
  # get force of infection (lambda) --------------------------------------------
  calendar_day <- lubridate::yday(as.Date(times, origin = params$calendar_start_date))
  beta_t <- params$beta * (1 + params$beta1 * cos(2 * pi * calendar_day / 365.24)) 
  lambda <- get_foi(x  = seir_output, 
                    y1 = params$eta_trans1[times,-1], 
                    y2 = params$eta_trans2[times,-1], 
                    y3 = params$eta_trans3[times,-1], 
                    y4 = params$eta_trans4[times,-1], 
                    y5 = params$eta_trans5[times,-1],
                    beta = beta_t, 
                    contact_mat = params$contact_mat,
                    times = times)
  
  # calculate infections -------------------------------------------------------
  new_infections <- lambda * (seir_output$S + seir_output$Shold_1d +
    (params$eta1[times,-1] * (seir_output$Sv_1d + seir_output$Shold_2d)) +
    (params$eta2[times,-1] * (seir_output$Sv_2d + seir_output$Shold_3d)) +
    (params$eta3[times,-1] * (seir_output$Sv_3d + seir_output$Shold_4d)) +
    (params$eta4[times,-1] * (seir_output$Sv_4d + seir_output$Shold_5d)) +
    (params$eta5[times,-1] * seir_output$Sv_5d)
  )
  
  # calculate cases ------------------------------------------------------------
  new_cases <- sweep(params$sigma * (seir_output$E + seir_output$Ev_1d + 
    seir_output$Ev_2d + seir_output$Ev_3d + seir_output$Ev_4d + seir_output$Ev_5d),
    2, params$p_report, "*")

  # calculate hospital admissions ----------------------------------------------
  hosp_admissions <- sweep(seir_output$I + 
    params$eta_hosp1[times,-1] * seir_output$Iv_1d +
    params$eta_hosp2[times,-1] * seir_output$Iv_2d +
    params$eta_hosp3[times,-1] * seir_output$Iv_3d +
    params$eta_hosp4[times,-1] * seir_output$Iv_4d +
    params$eta_hosp5[times,-1] * seir_output$Iv_5d,
    2, params$h, "*")
  
  
  hosp_admissions <- sweep(infectious, 2, params$h, "*")
  hosp_occ <- (seir_output$H + seir_output$Hv_1d + seir_output$Hv_2d)
  ic <- sweep(hosp_occ, 2, params$i1, "*")
  ic_occ <- (seir_output$IC + seir_output$ICv_1d + seir_output$ICv_2d)
  hosp_after_ic <- sweep(ic, 2, params$i2, "*")
  hosp_after_ic_occ <- (seir_output$H_IC + seir_output$H_ICv_1d + seir_output$H_ICv_2d)
  total_hosp_occ <- (seir_output$H + seir_output$Hv_1d + seir_output$Hv_2d) + (seir_output$H_IC + seir_output$H_ICv_1d + seir_output$H_ICv_2d)

  # deaths
  daily_deaths <- sweep(ic_occ, 2, params$d_ic, "*") + sweep(hosp_occ, 2, params$d, "*") + sweep(hosp_after_ic_occ, 2, params$d_hic, "*")

  # Create object for plotting ---------------------------------------
  # convert from wide to long format

  inf_long <- infections %>%
    mutate(time = times) %>%
    pivot_longer(
      cols = starts_with("E"),
      names_to = "age_group",
      names_prefix = "E",
      values_to = "infections"
    )

  new_inf_long <- new_infections %>%
    mutate(time = times) %>%
    pivot_longer(
      cols = starts_with("S"),
      names_to = "age_group",
      names_prefix = "S",
      values_to = "new_infections"
    )

  cases_long <- cases %>%
    mutate(time = times) %>%
    pivot_longer(
      cols = starts_with("I"),
      names_to = "age_group",
      names_prefix = "I",
      values_to = "cases"
    )

  new_cases_long <- new_cases %>%
    mutate(time = times) %>%
    pivot_longer(
      cols = starts_with("E"),
      names_to = "age_group",
      names_prefix = "E",
      values_to = "new_cases"
    )

  hosp_long <- total_hosp_occ %>%
    mutate(time = times) %>%
    pivot_longer(
      cols = starts_with("H"),
      names_to = "age_group",
      names_prefix = "H",
      values_to = "hospitalisations"
    )

  hosp_admin_long <- hosp_admissions %>%
    mutate(time = times) %>%
    pivot_longer(
      cols = starts_with("I"),
      names_to = "age_group",
      names_prefix = "I",
      values_to = "hospital_admissions"
    )

  ic_long <- ic %>%
    mutate(time = times) %>%
    pivot_longer(
      cols = starts_with("H"),
      names_to = "age_group",
      names_prefix = "H",
      values_to = "ic_admissions"
    )

  deaths_long <- daily_deaths %>%
    mutate(time = times) %>%
    pivot_longer(
      cols = starts_with("IC"),
      names_to = "age_group",
      names_prefix = "IC",
      values_to = "new_deaths"
    )


  df <- left_join(inf_long, new_inf_long, by = c("time", "age_group")) %>%
    left_join(.data, cases_long, by = c("time", "age_group")) %>%
    left_join(.data, new_cases_long, by = c("time", "age_group")) %>%
    left_join(.data, hosp_long, by = c("time", "age_group")) %>%
    left_join(.data, hosp_admin_long, by = c("time", "age_group")) %>%
    left_join(.data, ic_long, by = c("time", "age_group")) %>%
    left_join(.data, deaths_long, by = c("time", "age_group")) %>%
    pivot_longer(
      cols = c(
        "infections", "new_infections", "cases", "new_cases",
        "hospitalisations", "hospital_admissions", "ic_admissions",
        "new_deaths"
      ),
      names_to = "outcome",
      values_to = "value"
    ) %>%
    mutate(date = .data$time + as.Date(start_date)) %>%
    select(.data$time, .data$date, .data$age_group, .data$outcome, .data$value)


  df_summary <- df %>%
    group_by(.data$time, .data$date, .data$outcome) %>%
    summarise_at(.vars = "value", .funs = sum)

  rtn <- list(
    lambda = lambda_est,
    df = df,
    df_summary = df_summary
    # df_vac = df_vac
  )
  return(rtn)
}
