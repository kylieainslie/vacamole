# -----------------------------------------------------------
# Figure 1 script
# simulated outcomes in whole pop w/o waning 12+ vs. 18+
# -----------------------------------------------------------

# load packages ---------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
# -----------------------------------------------------------
source("R/forward_sim_func_wrap.R")
# read in simulation results --------------------------------
file_date <- "2021-10-01"
file_path <- "inst/extdata/results/main_analysis/"
# alpha, no wane
alpha_12plus <- readRDS(paste0(file_path, "results_12plus_alpha_", file_date, ".rds"))
alpha_18plus <- readRDS(paste0(file_path, "results_18plus_alpha_", file_date, ".rds"))
#upper_res_12plus <- readRDS(paste0("inst/extdata/results/results_12plus_upper_beta_", file_date, ".rds"))

# delta, no wane
delta_12plus <- readRDS(paste0(file_path, "results_12plus_delta_", file_date, ".rds"))
delta_18plus <- readRDS(paste0(file_path, "results_18plus_delta_", file_date, ".rds"))
#upper_res_18plus <- readRDS(paste0("inst/extdata/results/results_18plus_upper_beta_", file_date, ".rds"))

# wrangle raw results ----------------------------------------
# 12 +
alpha_12plus_wrangled <- wrangle_results(alpha_12plus) %>%
  mutate(Scenario = "12+", Variant = "Alpha")
alpha_18plus_wrangled <- wrangle_results(alpha_18plus) %>%
  mutate(Scenario = "18+", Variant = "Alpha")
delta_12plus_wrangled <- wrangle_results(delta_12plus) %>%
  mutate(Scenario = "12+", Variant = "Delta")
delta_18plus_wrangled <- wrangle_results(delta_18plus) %>%
  mutate(Scenario = "18+", Variant = "Delta")

data_combined <- bind_rows(alpha_12plus_wrangled, alpha_18plus_wrangled
                           #delta_12plus_wrangled, delta_18plus_wrangled
                           ) %>%
  pivot_wider(names_from = state, values_from = mean:upper) %>%
  group_by(Variant, Scenario, time) %>%
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
  select(Variant, Scenario, time, age_group, cases_unvac_mle:deaths_fullvac_upper) %>%
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
         ), levels = c("Daily Cases","Hospital Admissions","IC Admissions","Daily Deaths"))) %>%
  select(-time)

# figure 1a - 12+ vs. 18+, 10-19 age group, no waning ---------
fig1a <- ggplot(data = all_res_for_plot %>%
                        filter(age_group == 2,
                               outcome != "Daily Deaths") %>%
                        group_by(Variant, Scenario, date, outcome) %>%
                        summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum"), 
                aes(x = date, y = mle, fill = Variant, linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Variant), alpha = 0.3) +
  geom_line(aes(color = Variant)) +
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
fig1a

# figure 1b - 12+ vs. 18+, !10-19 age group, no waning --------
fig1b <- ggplot(data = all_res_for_plot %>%
                        filter(age_group != 2,
                               outcome != "Daily Deaths") %>%
                        group_by(Variant, Scenario, date, outcome) %>%
                        summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum"), 
                aes(x = date, y = mle, fill = Variant,linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Variant), alpha = 0.3) +
  geom_line(aes(color = Variant)) +
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
fig1b

fig1_no_legend <- plot_grid(fig1a + theme(legend.position = "none"), 
                  fig1b + theme(legend.position = "none"), 
                  labels = "AUTO", nrow = 1)

legend <- get_legend(
  fig1a + theme(legend.box.margin = margin(0, 0, 0, 12))
)

fig1 <- plot_grid(fig1_no_legend, legend, rel_heights = c(3, .4), nrow = 2)
fig1

ggsave(filename = "inst/extdata/results/figure 1.jpg", plot = fig1,
       units = "in", height = 10, width = 12, dpi = 300)

# table 1 ----------------------------------------------------
table1_10_19 <- all_res_for_plot %>%
  filter(age_group == 2,
         outcome != "Daily Deaths") %>%
  group_by(Scenario, R0, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

table1_not_10_19 <- all_res_for_plot %>%
  filter(age_group != 2,
         outcome != "Daily Deaths") %>%
  group_by(Scenario, R0, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

# calculate percent differnce
table1_not_10_19_12plus <- table1_not_10_19 %>% filter(Scenario == "12+")
table1_not_10_19_18plus <- table1_not_10_19 %>% filter(Scenario == "18+")
perc_diff <- (table1_not_10_19_12plus[,4:6] * 100)/table1_not_10_19_18plus[,4:6] - 100
