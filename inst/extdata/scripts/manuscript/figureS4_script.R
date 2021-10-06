# -----------------------------------------------------------
# Figure SA script - sensitivity analysis
# simulated outcomes in whole pop w/o waning 12+ vs. 18+
# when vac of 12-17 starts in Oct or Nov
# -----------------------------------------------------------

# load packages ---------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
# source files ----------------------------------------------
source("inst/extdata/scripts/helpers/model_run_helper.R")
source("R/forward_sim_func_wrap.R")

# set working directory -------------------------------------
setwd("inst/extdata/results/sensitivity_analysis/")

# read in simulation results --------------------------------
file_date <- "2021-09-16"

# 12+ - October
mle_res_12plus_oct <- readRDS(paste0("results_12plus_oct_mle_beta_", file_date, ".rds"))
lower_res_12plus_oct <- readRDS(paste0("results_12plus_oct_lower_beta_", file_date, ".rds"))
upper_res_12plus_oct <- readRDS(paste0("results_12plus_oct_upper_beta_", file_date, ".rds"))

# 12+ - November
mle_res_12plus_nov <- readRDS(paste0("results_12plus_nov_mle_beta_", file_date, ".rds"))
lower_res_12plus_nov <- readRDS(paste0("results_12plus_nov_lower_beta_", file_date, ".rds"))
upper_res_12plus_nov <- readRDS(paste0("results_12plus_nov_upper_beta_", file_date, ".rds"))

# 12+ - June
file_date2 <- "2021-09-22"
mle_res_12plus_jun <- readRDS(paste0("results_12plus_jun_mle_beta_", file_date2, ".rds"))
lower_res_12plus_jun <- readRDS(paste0("results_12plus_jun_lower_beta_", file_date2, ".rds"))
upper_res_12plus_jun <- readRDS(paste0("results_12plus_jun_upper_beta_", file_date2, ".rds"))

# 18+
mle_res_18plus <- readRDS(paste0("results_18plus_nov_mle_beta_", file_date, ".rds"))
lower_res_18plus <- readRDS(paste0("results_18plus_nov_lower_beta_", file_date, ".rds"))
upper_res_18plus <- readRDS(paste0("results_18plus_nov_upper_beta_", file_date, ".rds"))

# wrangle raw results ----------------------------------------
# 12 + - June
mle_12plus_jun_all <- wrangle_results(mle_res_12plus_jun)
lower_12plus_jun_all <- wrangle_results(lower_res_12plus_jun)
upper_12plus_jun_all <- wrangle_results(upper_res_12plus_jun)

all_12plus_jun <- bind_rows(mle_12plus_jun_all, lower_12plus_jun_all, upper_12plus_jun_all, .id = "R0") %>%
  mutate(
    R0 = case_when(
      R0 == 1 ~ "4.6",
      R0 == 2 ~ "3.45",
      R0 == 3 ~ "5.75"
    ), Scenario = "12+ (June)") %>%
  pivot_wider(names_from = state, values_from = mean:upper) %>%
  group_by(Scenario, R0, time) %>%
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
  select(Scenario, R0, time, age_group, cases_unvac_mle:deaths_fullvac_upper) %>%
  pivot_longer(cases_unvac_mle:deaths_fullvac_upper, names_to = c("outcome", "vac_status", "estimate"), names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = "estimate", values_from = "value")

# 12 + - October
mle_12plus_oct_all <- wrangle_results(mle_res_12plus_oct)
lower_12plus_oct_all <- wrangle_results(lower_res_12plus_oct)
upper_12plus_oct_all <- wrangle_results(upper_res_12plus_oct)

all_12plus_oct <- bind_rows(mle_12plus_oct_all, lower_12plus_oct_all, upper_12plus_oct_all, .id = "R0") %>%
  mutate(
    R0 = case_when(
    R0 == 1 ~ "4.6",
    R0 == 2 ~ "3.45",
    R0 == 3 ~ "5.75"
  ), Scenario = "12+ (October)") %>%
  pivot_wider(names_from = state, values_from = mean:upper) %>%
  group_by(Scenario, R0, time) %>%
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
  select(Scenario, R0, time, age_group, cases_unvac_mle:deaths_fullvac_upper) %>%
  pivot_longer(cases_unvac_mle:deaths_fullvac_upper, names_to = c("outcome", "vac_status", "estimate"), names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = "estimate", values_from = "value")

# 12+ - November
mle_12plus_nov_all <- wrangle_results(mle_res_12plus_oct)
lower_12plus_nov_all <- wrangle_results(lower_res_12plus_oct)
upper_12plus_nov_all <- wrangle_results(upper_res_12plus_oct)

all_12plus_nov <- bind_rows(mle_12plus_nov_all, lower_12plus_nov_all, upper_12plus_nov_all, .id = "R0") %>%
  mutate(
    R0 = case_when(
      R0 == 1 ~ "4.6",
      R0 == 2 ~ "3.45",
      R0 == 3 ~ "5.75"
    ), Scenario = "12+ (November)") %>%
  pivot_wider(names_from = state, values_from = mean:upper) %>%
  group_by(Scenario, R0, time) %>%
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
  select(Scenario, R0, time, age_group, cases_unvac_mle:deaths_fullvac_upper) %>%
  pivot_longer(cases_unvac_mle:deaths_fullvac_upper, names_to = c("outcome", "vac_status", "estimate"), names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = "estimate", values_from = "value")

# 18+
mle_18plus_all <- wrangle_results(mle_res_18plus)
lower_18plus_all <- wrangle_results(lower_res_18plus)
upper_18plus_all <- wrangle_results(upper_res_18plus)

all_18plus <- bind_rows(mle_18plus_all, lower_18plus_all, upper_18plus_all, .id = "R0") %>%
  mutate(R0 = case_when(
    R0 == 1 ~ "4.6",
    R0 == 2 ~ "3.45",
    R0 == 3 ~ "5.75"
  ), Scenario = "18+") %>%
  pivot_wider(names_from = state, values_from = mean:upper) %>%
  group_by(Scenario, R0, time) %>%
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
  select(Scenario, R0, time, age_group, cases_unvac_mle:deaths_fullvac_upper) %>%
  pivot_longer(cases_unvac_mle:deaths_fullvac_upper, names_to = c("outcome", "vac_status", "estimate"), names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = "estimate", values_from = "value")

# make plots --------------------------------------------------
all_res_for_plot <- bind_rows(all_12plus_jun, all_12plus_oct, all_12plus_nov, all_18plus) %>%
  ungroup() %>%
  mutate(date = time + as.Date("2020-01-01"),
         outcome = factor(case_when(
            outcome == "cases" ~ "Daily Cases",
            outcome == "hosp" ~ "Hospital Admissions",
            outcome == "ic" ~ "IC Admissions",
            outcome == "deaths" ~ "Daily Deaths"
         ), levels = c("Daily Cases","Hospital Admissions","IC Admissions","Daily Deaths"))) %>%
  select(-time)

# figure 1a - 12+ vs. 18+, 10-19 age group, no waning ---------
figS4a <- ggplot(data = all_res_for_plot %>%
                        filter(age_group == 2,
                               outcome != "Daily Deaths"
                               #outcome == "Daily Cases"
                               ) %>%
                        group_by(Scenario, R0, date, outcome) %>%
                        summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum"), 
                aes(x = date, y = mle, fill = R0, linetype = Scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = R0), alpha = 0.3) +
  labs(y = "Value", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) +
  facet_wrap(~outcome, scales = "free_y", nrow = 4)
figS4a

# figure 1b - 12+ vs. 18+, !10-19 age group, no waning --------
figS4b <- ggplot(data = all_res_for_plot %>%
                        filter(age_group != 2,
                               outcome != "Daily Deaths") %>%
                        group_by(Scenario, R0, date, outcome) %>%
                        summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum"), 
                aes(x = date, y = mle, fill = R0, 
                                              linetype = Scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = R0), alpha = 0.3) +
  labs(y = "Value", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) +
  facet_wrap(~outcome, scales = "free_y", nrow = 4)
figS4b

figS4_no_legend <- plot_grid(figS4a + theme(legend.position = "none"), 
                  figS4b + theme(legend.position = "none"), 
                  labels = "AUTO", nrow = 1)

legend <- get_legend(
  figS4a + theme(legend.box.margin = margin(0, 0, 0, 12))
)

figS4 <- plot_grid(figS4_no_legend, legend, rel_heights = c(3, .4), nrow = 2)
figS4

ggsave(filename = "../figures/figure S4.jpg", plot = figS4,
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
