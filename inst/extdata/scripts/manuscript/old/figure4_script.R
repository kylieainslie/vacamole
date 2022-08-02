# -----------------------------------------------------------
# Figure 4 script (or supplemental figure 5)
# simulated outcomes in whole pop w/o waning 5+ vs. 18+
# -----------------------------------------------------------

# source files -----------------------------------------------
source("inst/extdata/scripts/manuscript/figure1_script.R")

# read in simulation results --------------------------------
file_date <- "2021-10-07"
file_path <- "inst/extdata/results/sensitivity_analysis/"
# alpha, no wane
alpha_5plus <- readRDS(paste0(file_path, "results_5plus_alpha_", file_date, ".rds"))
delta_5plus <- readRDS(paste0(file_path, "results_5plus_delta_", file_date, ".rds"))

# wrangle raw results ----------------------------------------
# 5 +
alpha_5plus_wrangled <- wrangle_results(alpha_5plus) %>%
  mutate(Scenario = "5+", Variant = "Alpha")
delta_5plus_wrangled <- wrangle_results(delta_5plus) %>%
  mutate(Scenario = "5+", Variant = "Delta")

data_combined4 <- bind_rows(alpha_5plus_wrangled, alpha_18plus_wrangled,
                            delta_5plus_wrangled, delta_18plus_wrangled) %>%
  mutate(Scenario = factor(Scenario, levels = c("5+", "18+"))) %>%
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
all_res_for_plot4 <- data_combined4 %>%
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
fig4a <- ggplot(data = all_res_for_plot4 %>%
                  filter(age_group %in% c(1,2),
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
fig4a

# figure 1b - 12+ vs. 18+, !10-19 age group, no waning --------
fig4b <- ggplot(data = all_res_for_plot4 %>%
                  filter(age_group %in% c(3:9),
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
fig4b

fig4_no_legend <- plot_grid(fig4a + theme(legend.position = "none"), 
                            fig4b + theme(legend.position = "none"), 
                            labels = "AUTO", nrow = 1)

legend <- get_legend(
  fig4a + theme(legend.box.margin = margin(0, 0, 0, 12))
)

fig4 <- plot_grid(fig4_no_legend, legend, rel_heights = c(3, .4), nrow = 2)
fig4

ggsave(filename = "inst/extdata/results/figures/figure 4.jpg", plot = fig4,
       units = "in", height = 10, width = 12, dpi = 300)

# table 1 ----------------------------------------------------
table4_0_19 <- all_res_for_plot4 %>%
  filter(age_group %in% c(1,2),
         outcome != "Daily Deaths") %>%
  group_by(Variant, Scenario, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

table4_not_0_19 <- all_res_for_plot4 %>%
  filter(age_group %in% c(3:9),
         outcome != "Daily Deaths") %>%
  group_by(Variant, Scenario, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

# calculate percent differnce
table4_0_19_5plus <- table4_0_19 %>% filter(Scenario == "5+")
table4_0_19_18plus <- table4_0_19 %>% filter(Scenario == "18+")
abs_diff <- table4_0_19_18plus[,4:6] - table4_0_19_5plus[,4:6]
perc_diff <- (table4_0_19_5plus[,4:6] * 100)/table4_0_19_18plus[,4:6] - 100

table4_not_0_19_5plus <- table4_not_0_19 %>% filter(Scenario == "5+")
table4_not_0_19_18plus <- table4_not_0_19 %>% filter(Scenario == "18+")
abs_diff2 <- table4_not_0_19_18plus[,4:6] - table4_not_0_19_5plus[,4:6]
perc_diff2 <- (table4_not_0_19_12plus[,4:6] * 100)/table4_not_0_19_18plus[,4:6] - 100

# plot by age group -------------------------------------------
# a) no waning
dat_fig5 <- all_res_for_plot4 %>%
  filter(outcome != "Daily Deaths") %>%
  group_by(Variant, Scenario, date, age_group, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

fig5 <- ggplot(data = dat_fig5, 
               aes(x = date, y = mle, fill = age_group, linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group), alpha = 0.3) +
  geom_line(aes(color = age_group)) +
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
  facet_grid(outcome~Variant, scales = "free_y")
fig5

ggsave(filename = "inst/extdata/results/figure 5.jpg", plot = fig5,
       units = "in", height = 8, width = 12, dpi = 300)
