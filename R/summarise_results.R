#' Summarise output from SEIR model for plotting
#' @param seir_output list of dataframes from postprocessing
#' @param params list of parameter values
#' @param start_date calendar date of start of simulation
#' @return List of summary results
#' @keywords vacamole
#' @export

summarise_results <- function(seir_output, params, start_date){

res <- get_vac_rate_2(times, params$vac_schedule, params$ve, params$delay)
eta_dose1 <- res %>%
  select(time, age_group, eta_dose1) %>%
  pivot_wider(names_from = age_group, names_prefix = "eta_dose1_",
              values_from = eta_dose1)
eta_dose2 <- res %>%
  select(time, age_group, eta_dose2) %>%
  pivot_wider(names_from = age_group, names_prefix = "eta_dose2_",
              values_from = eta_dose2)

lambda_est <- get_foi(dat = seir_output, 
                      beta = params$beta, 
                      sigma = params$sigma,
                      i1 = params$i1,
                      p_report = params$p_report,
                      N = params$N, 
                      c_lockdown = params$c_lockdown, 
                      c_relaxed = params$c_relaxed,
                      c_very_relaxed = params$c_very_relaxed,
                      c_normal = params$c_normal,
                      thresh_l = params$thresh_l,
                      thresh_m = params$thresh_m,
                      thresh_u = params$thresh_u,
                      use_cases = params$use_cases,
                      force_relax = params$force_relax)

lambda_est1 <- lambda_est$lambda %>%
  pivot_wider(names_from = age_group, names_prefix = "age_group_", values_from = foi)

# infections
new_infections <- (out$S + out$Shold_1d + (eta_dose1[,-1] * (out$Sv_1d + out$Shold_2d)) + (eta_dose2[,-1] * out$Sv_2d)) * lambda_est1[,-1]
infections <- (out$E + out$Ev_1d + out$Ev_2d)

#infectious/cases
new_infectious <- params$sigma * (out$E + out$Ev_1d + out$Ev_2d)
new_cases <- sweep(new_infectious, 2, p_reported_by_age, "*")
infectious <- (out$I + out$Iv_1d + out$Iv_2d)
cases <- sweep(infectious, 2, p_reported_by_age, "*")

# hospitalisations/ic admissions
hosp_admissions <- sweep(infectious, 2, h, "*")
hosp_occ <- (out$H + out$Hv_1d + out$Hv_2d)
ic <- sweep(hosp_occ, 2, i1, "*")
ic_occ <- (out$IC + out$ICv_1d + out$ICv_2d)
hosp_after_ic <- sweep(ic, 2, i2, "*")
total_hosp_occ <- (out$H + out$Hv_1d + out$Hv_2d) + (out$H_IC + out$H_ICv_1d + out$H_ICv_2d)

# deaths
daily_deaths <- sweep(ic, 2, d_ic, "*") + sweep(hosp_admissions, 2, d, "*") + sweep(hosp_after_ic, 2, d_hic, "*")

# Create object for plotting ---------------------------------------
# convert from wide to long format
inf_long <- infections %>% 
  mutate(time = times) %>%
  pivot_longer(cols = starts_with("E"), 
               names_to = "age_group", 
               names_prefix = "E",
               values_to = "infections")

new_inf_long <- new_infections %>% 
  mutate(time = times) %>%
  pivot_longer(cols = starts_with("S"), 
               names_to = "age_group", 
               names_prefix = "S",
               values_to = "new_infections")

cases_long <- cases %>% 
  mutate(time = times) %>%
  pivot_longer(cols = starts_with("I"), 
               names_to = "age_group", 
               names_prefix = "I",
               values_to = "cases")

new_cases_long <- new_cases %>% 
  mutate(time = times) %>%
  pivot_longer(cols = starts_with("E"), 
               names_to = "age_group", 
               names_prefix = "E",
               values_to = "new_cases")

hosp_long <- total_hosp_occ %>%
  mutate(time = times) %>%
  pivot_longer(cols = starts_with("H"), 
               names_to = "age_group", 
               names_prefix = "H",
               values_to = "hospitalisations")

hosp_admin_long <- hosp_admissions %>%
  mutate(time = times) %>%
  pivot_longer(cols = starts_with("I"), 
               names_to = "age_group", 
               names_prefix = "I",
               values_to = "hospital_admissions")

ic_long <- ic %>%
  mutate(time = times) %>%
  pivot_longer(cols = starts_with("H"), 
               names_to = "age_group", 
               names_prefix = "H",
               values_to = "ic_admissions")

deaths_long <- daily_deaths %>%
  mutate(time = times) %>%
  pivot_longer(cols = starts_with("H"),
               names_to = "age_group",
               names_prefix = "H",
               values_to = "new_deaths")


df <- left_join(inf_long, new_inf_long, by = c("time", "age_group")) %>%
  left_join(.,cases_long, by = c("time", "age_group")) %>%
  left_join(., new_cases_long, by = c("time", "age_group")) %>%
  left_join(.,hosp_long, by = c("time", "age_group")) %>%
  left_join(., hosp_admin_long, by = c("time", "age_group")) %>%
  left_join(.,ic_long, by = c("time", "age_group")) %>%
  left_join(.,deaths_long, by = c("time", "age_group")) %>%
  pivot_longer(cols = c("infections", "new_infections", "cases", "new_cases",
                        "hospitalisations", "hospital_admissions", "ic_admissions",
                        "new_deaths"),
               names_to = "outcome",
               values_to = "value") %>%
  mutate(date = time + as.Date(start_date)) %>%
  select(time, date, age_group, outcome, value)


df_summary <- df %>%
  group_by(time, date, outcome) %>%
  summarise_at(.vars = "value", .funs = sum) %>%
  mutate(scenario = tag)

  rtn <- list(lambda = lambda_est,
              df = df,
              df_summary = df_summary)
  return(rtn)
}






