# -------------------------------------
# Figure 2 script
# plot of outcomes by age group
# -------------------------------------
library(ggsci)
library(gridExtra)

# -------------------------------------------------------------
save_path <- "/rivm/s/ainsliek/results/impact_vac/resubmission/"

# read in results df -----------------------------------------------------------
df_all <- readRDS("/rivm/s/ainsliek/results/impact_vac/resubmission/results_all.rds")

# make plots -------------------------------------------------------------------
all_res_for_plot <- df_all %>%
  ungroup() %>%
  mutate(
    outcome = factor(case_when(
      target_variable == "inc infection" ~ "Daily Infections",
      target_variable == "inc case" ~ "Daily Cases",
      target_variable == "inc hosp" ~ "Hospital Admissions",
      target_variable == "inc icu" ~ "IC Admissions",
      target_variable == "inc death" ~ "Daily Deaths"
    ), levels = c("Daily Infections", "Daily Cases", "Hospital Admissions", "IC Admissions", "Daily Deaths")),
    age_group2 = case_when(
      age_group == 1 ~ "0-9 years",
      age_group == 2 ~ "10-19 years",
      age_group %in% c(3:9) ~ ">19 years"
    ),
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
    ),
    scenario_id = factor(scenario_id, levels = c("Vaccination in 5+", 
                                                 "Vaccination in 12+",
                                                 "Vaccination in 18+"))
  )

# Figure 2: summarise by age_group levels --------------------------------------
dat_fig2 <- all_res_for_plot %>%
  group_by(scenario_id, outcome, date, sample, age_group) %>%
  summarise(sum = sum(value)) %>%
  ungroup() %>%
  group_by(scenario_id, outcome, date, age_group) %>%
  summarise(mean  = mean(sum),
            q025 = quantile(sum, probs = 0.025),
            q25  = quantile(sum, probs = 0.25),
            q75  = quantile(sum, probs = 0.75),
            q975 = quantile(sum, probs = 0.975)
  ) %>%
  select(date, scenario_id, age_group, outcome, mean:q975) %>%
  # calculate cases per 100,000 people
  mutate(pop_size = case_when(
    age_group == "0-9" ~ n_vec[1],
    age_group == "10-19" ~ n_vec[2],
    age_group == "20-29" ~ n_vec[3],
    age_group == "30-39" ~ n_vec[4],
    age_group == "40-49" ~ n_vec[5],
    age_group == "50-59" ~ n_vec[6],
    age_group == "60-69" ~ n_vec[7],
    age_group == "70-79" ~ n_vec[8],
    age_group == "80+"   ~ n_vec[9],
  ))

aaas_cols <- pal_aaas("default")(3)
my_colors <- c(aaas_cols[3], aaas_cols[2], gray.colors(7, start = 0.2, end = 0.8))

fig2 <- ggplot(data = dat_fig2 %>%
                 filter(#scenario_id == "Vaccination in 5+",
                   #age_group2 == "10-19 years",
                   date >= as.Date("2021-11-01"),
                   #date <= as.Date("2021-12-31"),
                   outcome == "Daily Cases"), 
                 aes(x = date, y = mean, fill = age_group)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = age_group), alpha = 0.3) +
  geom_line(aes(color = age_group)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  labs(y = "Daily Cases", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  guides(fill = guide_legend("Age Group"), color = guide_legend("Age Group")) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) +
  facet_grid(scenario_id~., scales = "free_y")
fig2

ggsave(filename = paste0(save_path, "figure2_new_color_palette.jpg"), plot = fig3,
       units = "in", height = 8, width = 12, dpi = 300)

# -------------------------------------------------------------
# b) hosp and IC admissions is supplemental figure
dat_figS4 <- all_res_for_plot %>%
  filter(outcome %in% c("Hospital Admissions", "IC Admissions"),
         Immunity == "No Waning",
         date >= as.Date("2021-11-01")) %>%
  group_by(Immunity, Scenario, date, age_group, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

figS4 <- ggplot(data = dat_figS4, 
                 aes(x = date, y = mle, fill = age_group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group), alpha = 0.3) +
  geom_line(aes(color = age_group)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  labs(y = "Value", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  guides(fill = guide_legend("Age Group"), color = guide_legend("Age Group")) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) +
  facet_grid(outcome~Scenario, scales = "free_y")
figS4

# add text annotation
# figS4 + 
#   annotate(
#     geom = "curve", x = as.Date("2022-03-07"), y = 500, xend = as.Date("2022-02-01"), yend = 400, 
#     curvature = .3, arrow = arrow(length = unit(2, "mm"))
#   ) +
#   annotate(geom = "text", x = as.Date("2022-03-08"), y = 500, label = "", hjust = "left")

# ggsave(filename = paste0(save_path, "figureS4.jpg"), plot = figS4,
#        units = "in", height = 8, width = 12, dpi = 300)
