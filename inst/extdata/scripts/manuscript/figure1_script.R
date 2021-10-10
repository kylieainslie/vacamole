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
source("inst/extdata/scripts/helpers/model_run_helper.R")
source("R/forward_sim_func_wrap.R")
# read in simulation results --------------------------------
file_date <- "2021-10-09"
file_path <- "inst/extdata/results/main_analysis/"
# alpha, no wane
# alpha_12plus <- readRDS(paste0(file_path, "results_12plus_alpha_", file_date, ".rds"))
# alpha_18plus <- readRDS(paste0(file_path, "results_18plus_alpha_", file_date, ".rds"))

# delta, no wane
delta_5plus  <- readRDS(paste0(file_path, "results_5plus_delta_", file_date, ".rds"))
delta_12plus <- readRDS(paste0(file_path, "results_12plus_delta_", file_date, ".rds"))
delta_18plus <- readRDS(paste0(file_path, "results_18plus_delta_", file_date, ".rds"))

# wrangle raw results ----------------------------------------
# alpha_12plus_wrangled <- wrangle_results(alpha_12plus) %>%
#   mutate(Scenario = "12+", Variant = "Alpha")
# alpha_18plus_wrangled <- wrangle_results(alpha_18plus) %>%
#   mutate(Scenario = "18+", Variant = "Alpha")
delta_5plus_wrangled  <- wrangle_results(delta_5plus) 
delta_12plus_wrangled <- wrangle_results(delta_12plus) 
delta_18plus_wrangled <- wrangle_results(delta_18plus) 

data_combined <- bind_rows(delta_5plus_wrangled,
                           delta_12plus_wrangled, 
                           delta_18plus_wrangled,
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

# figure 1a: 5+ vs. 12+ vs. 18+ by age group ----------------------
fig1a <- ggplot(data = all_res_for_plot %>%
                        filter(#age_group == 2,
                               outcome == "Daily Cases",
                               date >= as.Date("2021-11-01")) %>%
                        mutate(age_group2 = case_when(
                          age_group == 1 ~ "0-9",
                          age_group == 2 ~ "10-19",
                          age_group %in% c(3:9) ~ ">19"),
                          age_group2 = factor(age_group2, levels = c("0-9", "10-19", ">19"))) %>%
                        group_by(Scenario, age_group2, date, outcome) %>%
                        summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum"), 
                aes(x = date, y = mle, fill = age_group2,linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group2), alpha = 0.3) +
  geom_line(aes(color = age_group2), size = 1) +
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
  #facet_wrap(~outcome, scales = "free_y", nrow = 1)
fig1a

# data wrangling --------------------------------------------
table1 <- all_res_for_plot %>%
  filter(#age_group == 2,
    outcome != "Daily Deaths",
    date >= as.Date("2021-11-01")) %>%
  mutate(age_group2 = case_when(
    age_group == 1 ~ "0-9",
    age_group == 2 ~ "10-19",
    age_group %in% c(3:9) ~ ">19"),
    age_group2 = factor(age_group2, levels = c("0-9", "10-19", ">19"))) %>%
  group_by(Scenario, age_group2, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

# calculate percent differnce
table1 <- table1 %>%
  group_by(age_group2, outcome) %>%
  mutate(abs_diff = mle - mle[Scenario == "18+"],
         abs_diff_lower = lower - lower[Scenario == "18+"],
         abs_diff_upper = upper - upper[Scenario == "18+"],
         perc_diff = (mle * 100)/mle[Scenario == "18+"] - 100,
         perc_diff_lower = (lower * 100)/lower[Scenario == "18+"] - 100,
         perc_diff_upper = (upper * 100)/upper[Scenario == "18+"] - 100)


# bar plot of percent difference ----------------------------
fig1_inset <- ggplot(data = table1 %>%
                       filter(Scenario != "18+"), 
                     aes(x = outcome, y = abs(perc_diff), fill = age_group2)) +
  geom_bar(stat = "Identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = abs(perc_diff_upper), ymax = abs(perc_diff_lower), width = 0.2),
                position = position_dodge(0.9)) +
  labs(x = "Outcome", y = "Percent Reduction (%)", fill = "Age Group") +
  facet_wrap(~Scenario, nrow = 2) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face = "bold"))
fig1_inset

fig1 <- fig1a + annotation_custom(ggplotGrob(fig1_inset), 
                                  xmin = as.Date("2022-02-01"), xmax = as.Date("2022-04-04"),
                                  ymin = 15000, ymax = 100000)
fig1
ggsave(filename = "inst/extdata/results/figures/figure_1_w_inset.jpg", plot = fig1,
       units = "in", height = 10, width = 12, dpi = 300)


# old -------------------------------------------------------
# fig1_no_legend <- plot_grid(fig1a + theme(legend.position = "none") + ylim(0, 100000), 
#                   fig1b + theme(legend.position = "none") + ylim(0, 100000), 
#                   fig1c + theme(legend.position = "none"),
#                   labels = "AUTO", nrow = 1)
# 
# legend <- get_legend(
#   fig1a + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
# 
# fig1 <- plot_grid(fig1_no_legend, legend, rel_heights = c(3, .4), nrow = 2)
# fig1