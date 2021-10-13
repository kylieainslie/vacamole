# -----------------------------------------------------------
# Figure 2 script
# simulated outcomes w waning 5+ vs. 12+ vs. 18+
# -----------------------------------------------------------

# load packages ---------------------------------------------
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

source("R/forward_sim_func_wrap.R")
# read in simulation results --------------------------------
file_date <- "2021-10-09"
file_path <- "C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/vaccination_modelling/vacamole_files/code/vacamole/inst/extdata/results/sensitivity_analysis/"

#file_path <- "inst/extdata/results/sensitivity_analysis/"

# delta, wane
delta_5plus_wane  <- readRDS(paste0(file_path, "results_5plus_delta_wane_", file_date, ".rds"))
delta_12plus_wane <- readRDS(paste0(file_path, "results_12plus_wane_delta_", file_date, ".rds"))
delta_18plus_wane <- readRDS(paste0(file_path, "results_18plus_wane_delta_", file_date, ".rds"))

# wrangle raw results ----------------------------------------
delta_5plus_wane_wrangled  <- wrangle_results(delta_5plus_wane) 
delta_12plus_wane_wrangled <- wrangle_results(delta_12plus_wane) 
delta_18plus_wane_wrangled <- wrangle_results(delta_18plus_wane) 

data_combined_wane <- bind_rows(delta_5plus_wane_wrangled,
                           delta_12plus_wane_wrangled, 
                           delta_18plus_wane_wrangled,
                           .id = "Scenario"
) %>%
  mutate(Scenario = case_when(
    Scenario == 1 ~ "5+",
    Scenario == 2 ~ "12+",
    Scenario == 3 ~ "18+"), 
    Scenario = factor(Scenario, levels = c("5+", "12+", "18+"))) %>%
  pivot_wider(names_from = state, values_from = mean:upper) %>%
  group_by(Scenario, time) %>%
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
  select(Scenario, time, age_group, cases_unvac_mle:deaths_fullvac_upper) %>%
  pivot_longer(cases_unvac_mle:deaths_fullvac_upper, names_to = c("outcome", "vac_status", "estimate"), names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = "estimate", values_from = "value")

# make plots --------------------------------------------------
all_res_for_plot_wane <- data_combined_wane %>%
  ungroup() %>%
  mutate(date = time + as.Date("2020-01-01"),
         outcome = factor(case_when(
           outcome == "cases" ~ "Daily Cases",
           outcome == "hosp" ~ "Hospital Admissions",
           outcome == "ic" ~ "IC Admissions",
           outcome == "deaths" ~ "Daily Deaths"
         ), levels = c("Daily Cases","Hospital Admissions","IC Admissions","Daily Deaths")),
         Immunity = "Waning") %>%
  select(-time)

all_res_for_plot <- all_res_for_plot %>%
  mutate(Immunity = "No Waning")

all_res <- bind_rows(all_res_for_plot, all_res_for_plot_wane)

# figure 2a - 12+ vs. 18+, no waning vs. waning, 10-19 
# age group ---------------------------------------------------
dat_fig2 <- all_res %>%
  filter(outcome == "Daily Cases",
         date >= as.Date("2021-11-01")) %>%
  mutate(age_group2 = case_when(
    age_group == 1 ~ "0-9",
    age_group == 2 ~ "10-19",
    age_group %in% c(3:9) ~ ">19"),
    age_group2 = factor(age_group2, levels = c("0-9", "10-19", ">19"))) %>%
  group_by(Immunity, Scenario, age_group2, date, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

fig2a <- ggplot(data = dat_fig2 %>%
                  filter(Immunity == "No Waning"), 
                aes(x = date, y = mle, fill = age_group2,linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group2), alpha = 0.3) +
  geom_line(aes(color = age_group2), size = 1) +
  geom_line(data = dat_fig2 %>%
              filter(Immunity == "Waning") %>%
              mutate(age_group3 = age_group2),
            aes(x = date, y = mle, linetype = Scenario, color = age_group2), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
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
        axis.title=element_text(size=14,face="bold")) #+
#facet_wrap(~age_group2, scales = "free_y", nrow = 1)
fig2a

ggsave(filename = "inst/extdata/results/figure 2.jpg", plot = fig2a,
       units = "in", height = 10, width = 12, dpi = 300)

# data wrangling --------------------------------------------
table2 <- all_res_for_plot_wane %>%
  filter(#age_group == 2,
    outcome != "Daily Deaths") %>%
  mutate(age_group2 = case_when(
    age_group == 1 ~ "0-9",
    age_group == 2 ~ "10-19",
    age_group %in% c(3:9) ~ ">19"),
    age_group2 = factor(age_group2, levels = c("0-9", "10-19", ">19"))) %>%
  group_by(Scenario, age_group2, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum") %>%
  ungroup()

# calculate percent difference
table2a <- table2 %>%
  group_by(age_group2, outcome) %>%
  mutate(abs_diff = mle - mle[Scenario == "18+"],
         abs_diff_lower = lower - lower[Scenario == "18+"],
         abs_diff_upper = upper - upper[Scenario == "18+"],
         perc_diff = (mle * 100)/mle[Scenario == "18+"] - 100,
         perc_diff_lower = (lower * 100)/lower[Scenario == "18+"] - 100,
         perc_diff_upper = (upper * 100)/upper[Scenario == "18+"] - 100) %>%
  mutate_if(is.numeric, round, 1) %>%
  as.data.frame()

save_path <- "C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/vaccination_modelling/vacamole_files/results/main_analysis/"
write.csv(table2a, file = paste0(save_path, "table2.csv"))

# all age groups together
table2 %>%
  group_by(Scenario, outcome) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  group_by(outcome) %>%
  mutate(abs_diff = mle - mle[Scenario == "18+"],
         abs_diff_lower = lower - lower[Scenario == "18+"],
         abs_diff_upper = upper - upper[Scenario == "18+"],
         perc_diff = (mle * 100)/mle[Scenario == "18+"] - 100,
         perc_diff_lower = (lower * 100)/lower[Scenario == "18+"] - 100,
         perc_diff_upper = (upper * 100)/upper[Scenario == "18+"] - 100) %>%
  mutate_if(is.numeric, round, 1) %>%
  as.data.frame()
# cases by age group -----------------------------------
dat_figS_waning <- all_res %>%
  filter(outcome == "Daily Cases",
         date >= as.Date("2021-11-01")) %>%
  group_by(Immunity, Scenario, age_group, date, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

figS_waning <- ggplot(data = dat_figS_waning %>%
                  filter(Immunity == "Waning"), 
                aes(x = date, y = mle, fill = age_group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group), alpha = 0.3) +
  geom_line(aes(color = age_group), size = 1) +
  #scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(y = "Daily Cases", x = "Date") +
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
  facet_grid(Scenario~., scales = "free_y")
figS_waning

save_path_fig <- "C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/vaccination_modelling/vacamole_files/results/figures/"
ggsave(filename = paste0(save_path_fig, "figure_waning_by_age_group.jpg"), plot = figS_waning,
       units = "in", height = 10, width = 12, dpi = 300)


