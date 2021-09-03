# -----------------------------------------------------------
# Figure 2 script
# simulated outcomes in whole pop w waning 12+ vs. 18+
# -----------------------------------------------------------

# load packages ---------------------------------------------
library(dplyr)
library(ggplot2)
library(cowplot)

# read in simulation results --------------------------------
file_date <- "2021-09-02"

# 12+
mle_res_12plus_wane   <- readRDS(paste0("inst/extdata/results/results_12plus_wane_mle_beta_", file_date, ".rds"))
lower_res_12plus_wane <- readRDS(paste0("inst/extdata/results/results_12plus_wane_lower_beta_", file_date, ".rds"))
upper_res_12plus_wane <- readRDS(paste0("inst/extdata/results/results_12plus_wane_upper_beta_", file_date, ".rds"))

# 18+
mle_res_18plus_wane   <- readRDS(paste0("inst/extdata/results/results_18plus_wane_mle_beta_", file_date, ".rds"))
lower_res_18plus_wane <- readRDS(paste0("inst/extdata/results/results_18plus_wane_lower_beta_", file_date, ".rds"))
upper_res_18plus_wane <- readRDS(paste0("inst/extdata/results/results_18plus_wane_upper_beta_", file_date, ".rds"))

# wrangle raw results ----------------------------------------
# 12 +
mle_12plus_wane_all <- wrangle_results(mle_res_12plus_wane)
lower_12plus_wane_all <- wrangle_results(lower_res_12plus_wane)
upper_12plus_wane_all <- wrangle_results(upper_res_12plus_wane)

all_12plus_wane <- bind_rows(mle_12plus_wane_all, lower_12plus_wane_all, upper_12plus_wane_all, .id = "R0") %>%
  mutate(R0 = case_when(
    R0 == 1 ~ "4.6",
    R0 == 2 ~ "3.45",
    R0 == 3 ~ "5.75"
  ), Scenario = "12+") %>%
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
mle_18plus_wane_all <- wrangle_results(mle_res_18plus_wane)
lower_18plus_wane_all <- wrangle_results(lower_res_18plus_wane)
upper_18plus_wane_all <- wrangle_results(upper_res_18plus_wane)

all_18plus_wane <- bind_rows(mle_18plus_wane_all, lower_18plus_wane_all, upper_18plus_wane_all, .id = "R0") %>%
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
all_res_for_plot_wane <- bind_rows(all_12plus_wane, all_18plus_wane) %>%
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
fig2a <- ggplot(data = all_res %>%
                  filter(age_group == 2,
                         outcome == "Daily Cases",
                         R0 == "4.6") %>%
                  group_by(Scenario, Immunity, date, outcome) %>%
                  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum"), 
                aes(x = date, y = mle, fill = Immunity, linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Immunity), alpha = 0.3) +
  geom_line() +
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
fig2a

# figure 2b - 12+ vs. 18+, no waning vs. waning, !10-19 
# age group ---------------------------------------------------
fig2b <- ggplot(data = all_res %>%
                  filter(age_group != 2,
                         outcome == "Daily Cases",
                         R0 == "4.6") %>%
                  group_by(Scenario, Immunity, date, outcome) %>%
                  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum"), 
                aes(x = date, y = mle, fill = Immunity, 
                    linetype = Scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Immunity), alpha = 0.3) +
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
fig2b

fig2_no_legend <- plot_grid(fig2a + theme(legend.position = "none"), 
                            fig2b + theme(legend.position = "none"), 
                            labels = "AUTO", nrow = 2)

legend2 <- get_legend(
  fig2a + theme(legend.box.margin = margin(0, 0, 0, 12))
)

fig2ab <- plot_grid(fig2_no_legend, legend2, rel_heights = c(3, .4), nrow = 2)

# figure 2c - bar chart of cases by vac status ----------------
totals <- all_res %>%
  filter(outcome == "Daily Cases",
         R0 == "4.6") %>%
  group_by(Scenario, R0, Immunity) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

data_for_bar_plot <- all_res %>%
  filter(outcome == "Daily Cases",
         R0 == "4.6") %>%
  group_by(Scenario, R0, Immunity, vac_status) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum") %>%
  left_join(., totals, by = c("Scenario", "R0", "Immunity")) %>%
  mutate(mle_prop = mle.x/mle.y,
         lower_prop = lower.x/lower.y,
         upper_prop = upper.x/upper.y,
         vac_status = factor(case_when(
           vac_status == "unvac" ~ "Unvaccinated",
           vac_status == "partvac" ~ "Partially Vaccinated",
           vac_status == "fullvac" ~ "Fully Vaccinated"
         ), levels = c("Unvaccinated", "Partially Vaccinated", "Fully Vaccinated")))

fig2c <- ggplot(data = data_for_bar_plot, 
                aes(x=vac_status, y=mle_prop, fill=Immunity)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=lower_prop, ymax=upper_prop), width=.2,
                position=position_dodge(.9)) +
  ylim(0,1) +
  labs(y = "Proportion of Daily Cases", x = "Vaccination Status") +
  facet_wrap(~Scenario, nrow = 2) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold"))
fig2c

fig2 <- plot_grid(fig2ab, fig2c, labels = c("", "C"), rel_widths = c(1, 0.7), ncol = 2)
fig2

ggsave(filename = "inst/extdata/results/figure 2.jpg", plot = fig2,
       units = "in", height = 8, width = 12, dpi = 300)



