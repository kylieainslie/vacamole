# -----------------------------------------------------------
# Figure SA script - sensitivity analysis
# simulated outcomes in whole pop w/o waning 12+ vs. 18+
# for the Alpha variant
# -----------------------------------------------------------

# load packages ---------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
# source files ----------------------------------------------
source("inst/extdata/scripts/helpers/model_run_helper.R")
source("R/forward_sim_func_wrap.R")

# read in simulation results --------------------------------
file_date <- "2021-10-01"
file_path <- "/rivm/s/ainsliek/results/impact_vac/"
# no waning
alpha_12plus <- readRDS(paste0(file_path, "results_12plus_alpha_", file_date, ".rds"))
alpha_18plus <- readRDS(paste0(file_path, "results_18plus_alpha_", file_date, ".rds"))

# waning
alpha_12plus_wane <- readRDS(paste0(file_path, "results_12plus_alpha_wane_", file_date, ".rds"))
alpha_18plus_wane <- readRDS(paste0(file_path, "results_18plus_alpha_wane_", file_date, ".rds"))

# wrangle raw results ----------------------------------------
alpha_12plus_wrangled <- wrangle_results(alpha_12plus) %>%
  mutate(Scenario = "12+", Immunity = "No Waning")
alpha_18plus_wrangled <- wrangle_results(alpha_18plus) %>%
  mutate(Scenario = "18+", Immunity = "No Waning")
alpha_12plus_wane_wrangled <- wrangle_results(alpha_12plus_wane) %>%
  mutate(Scenario = "12+", Immunity = "Waning")
alpha_18plus_wane_wrangled <- wrangle_results(alpha_18plus_wane) %>%
  mutate(Scenario = "18+", Immunity = "Waning")

# combine results into single data frame ---------------------
dat_combined_s4 <- bind_rows(alpha_12plus_wrangled,
                             alpha_18plus_wrangled,
                             alpha_12plus_wane_wrangled,
                             alpha_18plus_wane_wrangled) %>%
  pivot_wider(names_from = state, values_from = mean:upper) %>%
  group_by(Scenario, Immunity, time) %>%
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
  select(Scenario, Immunity, time, age_group, cases_unvac_mle:deaths_fullvac_upper) %>%
  pivot_longer(cases_unvac_mle:deaths_fullvac_upper, names_to = c("outcome", "vac_status", "estimate"), names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = "estimate", values_from = "value")

# for plot ----------------------------------------------------
dat_for_plot <- dat_combined_s4 %>%
  ungroup() %>%
  mutate(date = time + as.Date("2020-01-01"),
         outcome = factor(case_when(
            outcome == "cases" ~ "Daily Cases",
            outcome == "hosp" ~ "Hospital Admissions",
            outcome == "ic" ~ "IC Admissions",
            outcome == "deaths" ~ "Daily Deaths"
         ), levels = c("Daily Cases","Hospital Admissions","IC Admissions","Daily Deaths")),
         age_group2 = case_when(
           age_group == 2 ~ "10-19 year olds",
           age_group %in% c(1,3:9) ~ "0-9 & >19 year olds"
         ),
         Scenario = case_when(
           Scenario == "12+" ~ "Vaccination of 12+", 
           Scenario == "18+" ~ "Vaccination of 18+"
         ),
         Scenario = factor(Scenario, levels = c("Vaccination of 12+", "Vaccination of 18+")),
         Immunity = case_when(
           Immunity == "No Waning" ~ "No waning immunity",
           Immunity == "Waning" ~ "Waning immunity"
         )
         ) %>%
  select(-time)

dat_figS4 <- dat_for_plot %>%
  group_by(Scenario, Immunity, age_group2, date, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum") %>%
  ungroup()

# figure S4 --------------------------------------------------
# 12+ vs. 18+, 10-19 age group, no waning --------------------

figS4 <- ggplot(data = dat_figS4 %>%
                   filter(outcome == "Daily Cases"),
                aes(x = date, y = mle, fill = Immunity, linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Immunity), alpha = 0.3) +
  geom_line(aes(color = Immunity), size = 0.7) +
  labs(y = "Daily Cases", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  guides(fill = guide_legend(""), color = guide_legend(""),
         linetype = guide_legend("Strategy")) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) +
  facet_grid(age_group2~., scales = "free_y")
figS4

ggsave(filename = paste0(file_path,"figure S4.jpg"), plot = figS4,
       units = "in", height = 10, width = 12, dpi = 300)

# figure S5 --------------------------------------------------
# 12+ vs. 18+, 10-19 age group, no waning --------------------
figS5 <- ggplot(data = dat_figS4 %>%
                   filter(outcome %in% c("Hospital Admissions", "IC Admissions")), 
                aes(x = date, y = mle, fill = Immunity,linetype = Scenario)) +
  geom_line(aes(color = Immunity)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Immunity), alpha = 0.3) +
  labs(y = "Value", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  guides(fill = guide_legend(""), color = guide_legend(""),
         linetype = guide_legend("Strategy")) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) +
  facet_grid(age_group2~outcome, scales = "free_y")
figS5

ggsave(filename = paste0(file_path,"figure S5.jpg"), plot = figS5,
       units = "in", height = 10, width = 12, dpi = 300)

# table S4 ----------------------------------------------------
tableS4_10_19 <- all_res_for_plot %>%
  filter(age_group == 2,
         outcome != "Daily Deaths") %>%
  group_by(Scenario, R0, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

tableS4_not_10_19 <- all_res_for_plot %>%
  filter(age_group != 2,
         outcome != "Daily Deaths") %>%
  group_by(Scenario, R0, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

# calculate percent differnce
tableS4_not_10_19_12plus <- table1_not_10_19 %>% filter(Scenario == "12+")
tableS4_not_10_19_18plus <- table1_not_10_19 %>% filter(Scenario == "18+")
perc_diff <- (tableS4_not_10_19_12plus[,4:6] * 100)/tableS4_not_10_19_18plus[,4:6] - 100
