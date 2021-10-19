# -----------------------------------------------------------
# Data wrangling for manuscript figures script
# -----------------------------------------------------------

# load packages ---------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
# -----------------------------------------------------------
source("inst/extdata/scripts/helpers/model_run_helper.R")
source("R/forward_sim_func_wrap.R")

# read in simulation results --------------------------------
file_date <- "2021-10-09"
#file_path <- "C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/vaccination_modelling/vacamole_files/code/vacamole/inst/extdata/results/main_analysis/"
file_path <- "/rivm/s/ainsliek/results/impact_vac/"

# delta, no wane
delta_5plus  <- readRDS(paste0(file_path, "results_5plus_delta_", file_date, ".rds"))
delta_12plus <- readRDS(paste0(file_path, "results_12plus_delta_", file_date, ".rds"))
delta_18plus <- readRDS(paste0(file_path, "results_18plus_delta_", file_date, ".rds"))

# delta, wane
delta_5plus_wane  <- readRDS(paste0(file_path, "results_5plus_delta_wane_", file_date, ".rds"))
delta_12plus_wane <- readRDS(paste0(file_path, "results_12plus_wane_delta_", file_date, ".rds"))
delta_18plus_wane <- readRDS(paste0(file_path, "results_18plus_wane_delta_", file_date, ".rds"))

# wrangle raw results ----------------------------------------
delta_5plus_wrangled  <- wrangle_results(delta_5plus) %>%
  mutate(Immunity = "No Waning", Scenario = "Vaccination of 5+")
delta_12plus_wrangled <- wrangle_results(delta_12plus) %>%
  mutate(Immunity = "No Waning", Scenario = "Vaccination of 12+")
delta_18plus_wrangled <- wrangle_results(delta_18plus) %>%
  mutate(Immunity = "No Waning", Scenario = "Vaccination of 18+")
delta_5plus_wane_wrangled  <- wrangle_results(delta_5plus_wane) %>%
  mutate(Immunity = "Waning", Scenario = "Vaccination of 5+")
delta_12plus_wane_wrangled <- wrangle_results(delta_12plus_wane) %>%
  mutate(Immunity = "Waning", Scenario = "Vaccination of 12+")
delta_18plus_wane_wrangled <- wrangle_results(delta_18plus_wane) %>%
  mutate(Immunity = "Waning", Scenario = "Vaccination of 18+")

# combine wrangled data into single data set ------------------
data_combined <- bind_rows(delta_5plus_wrangled,
                           delta_12plus_wrangled, 
                           delta_18plus_wrangled,
                           delta_5plus_wane_wrangled,
                           delta_12plus_wane_wrangled, 
                           delta_18plus_wane_wrangled) %>%
  mutate(Scenario = factor(Scenario, levels = c("Vaccination of 5+", "Vaccination of 12+", "Vaccination of 18+"))) %>%
  pivot_wider(names_from = state, values_from = mean:upper) %>%
  group_by(Immunity, Scenario, time) %>%
  mutate(h = params$h,
         i1 = params$i1,
         i2 = params$i2,
         d = params$d,
         d_ic = params$d_ic,
         d_hic = params$d_hic,
         cases_unvac_mle = params$sigma * mean_E * params$p_report,
         cases_unvac_lower = params$sigma * lower_E * params$p_report,
         cases_unvac_upper = params$sigma * upper_E * params$p_report,
         cases_partvac_mle = params$sigma * mean_Ev_1d * params$p_report,
         cases_partvac_lower = params$sigma * lower_Ev_1d * params$p_report,
         cases_partvac_upper = params$sigma * upper_Ev_1d * params$p_report,
         cases_fullvac_mle = params$sigma * mean_Ev_2d * params$p_report,
         cases_fullvac_lower = params$sigma * lower_Ev_2d * params$p_report,
         cases_fullvac_upper = params$sigma * upper_Ev_2d * params$p_report,
         hosp_unvac_mle = h * mean_I,
         hosp_unvac_lower = h * lower_I,
         hosp_unvac_upper = h * upper_I,
         hosp_partvac_mle = h * mean_Iv_1d,
         hosp_partvac_lower = h * lower_Iv_1d,
         hosp_partvac_upper = h * upper_Iv_1d,
         hosp_fullvac_mle = h * mean_Iv_2d,
         hosp_fullvac_lower = h * lower_Iv_2d,
         hosp_fullvac_upper = h * upper_Iv_2d,
         ic_unvac_mle = i1 * mean_H,
         ic_unvac_lower = i1 * lower_H,
         ic_unvac_upper = i1 * upper_H,
         ic_partvac_mle = i1 * mean_Hv_1d,
         ic_partvac_lower = i1 * lower_Hv_1d,
         ic_partvac_upper = i1 * upper_Hv_1d,
         ic_fullvac_mle = i1 * mean_Hv_2d,
         ic_fullvac_lower = i1 * lower_Hv_2d,
         ic_fullvac_upper = i1 * upper_Hv_2d,
         deaths_unvac_mle = d * mean_H + d_ic * mean_IC + d_hic * mean_H_IC,
         deaths_unvac_lower = d * lower_H + d_ic * lower_IC + d_hic * lower_H_IC,
         deaths_unvac_upper = d * upper_H + d_ic * upper_I + d_hic * upper_H_IC,
         deaths_partvac_mle = d * mean_Hv_1d + d_ic * mean_ICv_1d + d_hic * mean_H_ICv_1d,
         deaths_partvac_lower = d * lower_Hv_1d + d_ic * lower_ICv_1d + d_hic * lower_H_ICv_1d,
         deaths_partvac_upper = d * upper_Hv_1d + d_ic * upper_ICv_1d + d_hic * upper_H_ICv_1d,
         deaths_fullvac_mle = d * mean_Hv_2d + d_ic * mean_ICv_2d + d_hic * mean_H_ICv_2d,
         deaths_fullvac_lower = d * lower_Hv_2d + d_ic * lower_ICv_2d + d_hic * lower_H_ICv_2d,
         deaths_fullvac_upper = d * upper_Hv_2d + d_ic * upper_ICv_2d + d_hic * upper_H_ICv_2d
  ) %>%
  select(Immunity, Scenario, time, age_group, cases_unvac_mle:deaths_fullvac_upper) %>%
  pivot_longer(cases_unvac_mle:deaths_fullvac_upper, names_to = c("outcome", "vac_status", "estimate"), names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = "estimate", values_from = "value")


# make plots --------------------------------------------------
all_res_for_plot <- data_combined %>%
  ungroup() %>%
  mutate(date = time + as.Date("2020-01-01"),
         outcome = factor(case_when(
           outcome == "cases" ~ "Daily Cases",
           outcome == "hosp" ~ "Hospital Admissions",
           outcome == "ic" ~ "IC Admissions",
           outcome == "deaths" ~ "Daily Deaths"
         ), levels = c("Daily Cases","Hospital Admissions","IC Admissions","Daily Deaths")),
         age_group2 = case_when(
           age_group == 1 ~ "0-9 years",
           age_group == 2 ~ "10-19 years",
           age_group %in% c(3:9) ~ ">19 years"),
         age_group2 = factor(age_group2, levels = c("0-9 years", "10-19 years", ">19 years")),
         age_group = case_when(
           age_group == 1 ~ "0-9",
           age_group == 2 ~ "10-19",
           age_group == 3 ~ "20-29",
           age_group == 4 ~ "30-39",
           age_group == 5 ~ "40-49",
           age_group == 6 ~ "50-59",
           age_group == 7 ~ "60-69",
           age_group == 8 ~ "70-79",
           age_group == 9 ~ "80+"
         )
  )

# save output ---------------------------------------------------
#saveRDS(all_res_for_plot, file = "inst/extdata/results/all_res_for_plot.rds")
