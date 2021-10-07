# -----------------------------------------------------------
# Figure 2 script
# simulated outcomes in whole pop w waning 12+ vs. 18+
# -----------------------------------------------------------

# load packages ---------------------------------------------
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

source("R/forward_sim_func_wrap.R")
# read in simulation results --------------------------------
file_date <- "2021-10-01"
file_path <- "inst/extdata/results/main_analysis/"
# alpha
alpha_12plus_wane   <- readRDS(paste0(file_path,"results_12plus_wane_alpha_", file_date, ".rds"))
alpha_18plus_wane <- readRDS(paste0(file_path,"results_18plus_wane_alpha_", file_date, ".rds"))
#upper_res_12plus_wane <- readRDS(paste0("inst/extdata/results/results_12plus_wane_upper_beta_", file_date, ".rds"))

# delta
delta_12plus_wane   <- readRDS(paste0(file_path, "results_12plus_wane_delta_", file_date, ".rds"))
delta_18plus_wane <- readRDS(paste0(file_path, "results_18plus_wane_delta_", file_date, ".rds"))
#upper_res_18plus_wane <- readRDS(paste0("inst/extdata/results/results_18plus_wane_upper_beta_", file_date, ".rds"))

# wrangle raw results ----------------------------------------
# 12 +
alpha_12plus_wane_wrangled <- wrangle_results(alpha_12plus_wane) %>%
  mutate(Scenario = "12+", Variant = "Alpha")
alpha_18plus_wane_wrangled <- wrangle_results(alpha_18plus_wane) %>%
  mutate(Scenario = "18+", Variant = "Alpha")
delta_12plus_wane_wrangled <- wrangle_results(delta_12plus_wane) %>%
  mutate(Scenario = "12+", Variant = "Delta")
delta_18plus_wane_wrangled <- wrangle_results(delta_18plus_wane) %>%
  mutate(Scenario = "18+", Variant = "Delta")

data_wane_combined <- bind_rows(alpha_12plus_wane_wrangled, 
                                alpha_18plus_wane_wrangled, 
                                delta_12plus_wane_wrangled,
                                delta_18plus_wane_wrangled) %>%
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
all_res_for_plot_wane <- data_wane_combined %>%
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
dat_fig2a <- all_res %>%
  filter(age_group == 2,
         outcome == "Daily Cases") %>%
  group_by(Variant, Scenario, Immunity, date, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

fig2a <- ggplot(data = dat_fig2a, 
                aes(x = date, y = mle, fill = Immunity, linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Immunity), alpha = 0.3) +
  geom_line(aes(color = Immunity)) +
  labs(y = "Daily Cases", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_blank(),#lement_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) +
  facet_wrap(~Variant, scales = "free_y", nrow = 1)
fig2a

# figure 2b - 12+ vs. 18+, no waning vs. waning, !10-19 
# age group ---------------------------------------------------
dat_fig2b <- all_res %>%
  filter(age_group != 2,
         outcome == "Daily Cases") %>%
  group_by(Variant, Scenario, Immunity, date, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

fig2b <- ggplot(data = dat_fig2b, 
                aes(x = date, y = mle, fill = Immunity, 
                    linetype = Scenario)) +
  geom_line(aes(color = Immunity)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Immunity), alpha = 0.3) +
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
  facet_wrap(~Variant, scales = "free_y", nrow = 1)
fig2b

fig2_no_legend <- plot_grid(fig2a + theme(legend.position = "none"), 
                            fig2b + theme(legend.position = "none"), 
                            labels = "AUTO", nrow = 2, rel_heights = c(0.65,1))

legend2 <- get_legend(
  fig2a + theme(legend.box.margin = margin(0, 0, 0, 12))
)

fig2ab <- plot_grid(fig2_no_legend, legend2, rel_heights = c(3, .4), nrow = 2)
fig2ab

ggsave(filename = "inst/extdata/results/figure 2 new.jpg", plot = fig2ab,
       units = "in", height = 10, width = 12, dpi = 300)

# table 2 -----------------------------------------------------
table2_10_19_wane <- all_res_for_plot_wane %>%
  filter(age_group == 2,
         outcome != "Daily Deaths") %>%
  group_by(Variant, Scenario, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

table2_not_10_19_wane <- all_res_for_plot_wane %>%
  filter(age_group != 2,
         outcome != "Daily Deaths") %>%
  group_by(Variant, Scenario, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

# calculate percent difference_wane
# 10-19
table2_10_19_12plus_wane <- table2_10_19_wane %>% filter(Scenario == "12+")
table2_10_19_18plus_wane <- table2_10_19_wane %>% filter(Scenario == "18+")
(table2_10_19_12plus_wane[,4:6] * 100)/table2_10_19_18plus_wane[,4:6] - 100

# not 10-19
table2_not_10_19_12plus_wane <- table2_not_10_19_wane %>% filter(Scenario == "12+")
table2_not_10_19_18plus_wane <- table2_not_10_19_wane %>% filter(Scenario == "18+")
(table2_not_10_19_12plus_wane[,4:6] * 100)/table2_not_10_19_18plus_wane[,4:6] - 100




