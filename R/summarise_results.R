#' Summarise output from SEIR model for plotting
#' @param seir_output list of dataframes from postprocessing
#' @param params list of parameter values
#' @param start_date calendar date of start of simulation
#' @param times vector of time points
#' @param vac_inputs vaccination inputs (from the output of convert_vac_schedule())
#' @return List of summary results
#' @keywords vacamole
#' @export

summarise_results <- function(seir_output, params, start_date, times, vac_inputs){

alpha_dose1 <- vac_inputs$alpha_dose1[times,]
alpha_dose2 <- vac_inputs$alpha_dose2[times,]
eta_dose1 <- vac_inputs$eta_dose1[times,]
eta_dose2 <- vac_inputs$eta_dose2[times,]
delay_dose1 <- vac_inputs$delay_dose1[times,]
delay_dose2 <- vac_inputs$delay_dose2[times,]
eta_hosp_dose1 <- vac_inputs$eta_hosp_dose1[times,]
eta_hosp_dose2 <- vac_inputs$eta_hosp_dose2[times,]
eta_trans_dose1 <- vac_inputs$eta_trans_dose1[times,]
eta_trans_dose2 <- vac_inputs$eta_trans_dose2[times,]

vac_inputs_new = list(alpha_dose1 = alpha_dose1, 
                      alpha_dose2 = alpha_dose2, 
                      eta_dose1 = eta_dose1, 
                      eta_dose2 = eta_dose2,
                      delay_dose1 = delay_dose1,
                      delay_dose2 = delay_dose2,
                      eta_hosp_dose1 = eta_hosp_dose1, 
                      eta_hosp_dose2 = eta_hosp_dose2,
                      eta_trans_dose1 = eta_trans_dose1, 
                      eta_trans_dose2 = eta_trans_dose2)

lambda_est <- get_foi(dat = seir_output, params = params, vac_inputs = vac_inputs_new)

lambda_est1 <- lambda_est$lambda %>%
  pivot_wider(names_from = .data$age_group, names_prefix = "age_group_", values_from = .data$foi)

# infections
new_infections <- (seir_output$S + seir_output$Shold_1d + 
                     (eta_dose1[,-1] * (seir_output$Sv_1d + seir_output$Shold_2d)) + 
                     (eta_dose2[,-1] * seir_output$Sv_2d)) * lambda_est1[,-1]
infections <- (seir_output$E + seir_output$Ev_1d + seir_output$Ev_2d)

#infectious/cases
new_infectious <- params$sigma * (seir_output$E + seir_output$Ev_1d + seir_output$Ev_2d)
new_cases <- sweep(new_infectious, 2, params$p_report, "*")
infectious <- (seir_output$I + seir_output$Iv_1d + seir_output$Iv_2d)
cases <- sweep(infectious, 2, params$p_report, "*")

# hospitalisations/ic admissions
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
  pivot_longer(cols = starts_with("IC"),
               names_to = "age_group",
               names_prefix = "IC",
               values_to = "new_deaths")


df <- left_join(inf_long, new_inf_long, by = c("time", "age_group")) %>%
  left_join(.data,cases_long, by = c("time", "age_group")) %>%
  left_join(.data, new_cases_long, by = c("time", "age_group")) %>%
  left_join(.data,hosp_long, by = c("time", "age_group")) %>%
  left_join(.data, hosp_admin_long, by = c("time", "age_group")) %>%
  left_join(.data,ic_long, by = c("time", "age_group")) %>%
  left_join(.data,deaths_long, by = c("time", "age_group")) %>%
  pivot_longer(cols = c("infections", "new_infections", "cases", "new_cases",
                        "hospitalisations", "hospital_admissions", "ic_admissions",
                        "new_deaths"),
               names_to = "outcome",
               values_to = "value") %>%
  mutate(date = .data$time + as.Date(start_date)) %>%
  select(.data$time, .data$date, .data$age_group, .data$outcome, .data$value)


df_summary <- df %>%
  group_by(.data$time, .data$date, .data$outcome) %>%
  summarise_at(.vars = "value", .funs = sum)

  rtn <- list(lambda = lambda_est,
              df = df,
              df_summary = df_summary
              # df_vac = df_vac
              )
  return(rtn)
}






