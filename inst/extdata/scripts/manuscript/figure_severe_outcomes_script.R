# grid plot of severe disease outcomes by vaccination scenario

# load required packages -------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

# read in results df -----------------------------------------------------------
df_all <- readRDS("/rivm/s/ainsliek/results/impact_vac/resubmission/results_all.rds")

# make plots --------------------------------------------------
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

# summarise by age_group2 levels
all_res_for_plot_sum <- all_res_for_plot %>%
  group_by(scenario_id, outcome, date, sample, age_group2) %>%
  summarise(sum = sum(value)) %>%
  ungroup() %>%
  group_by(scenario_id, outcome, date, age_group2) %>%
  summarise(mean  = mean(sum),
            q025 = quantile(sum, probs = 0.025),
            q25  = quantile(sum, probs = 0.25),
            q75  = quantile(sum, probs = 0.75),
            q975 = quantile(sum, probs = 0.975)
  ) %>%
  select(date, scenario_id, age_group2, outcome, mean:q975) %>%
  # calculate cases per 100,000 people
  mutate(pop_size = case_when(
    age_group2 == "0-9 years" ~ n_vec[1],
    age_group2 == "10-19 years" ~ n_vec[2],
    age_group2 == ">19 years" ~ sum(n_vec[3:9])
  ))

all_res_total <- all_res_for_plot_sum %>%
  ungroup() %>%
  group_by(scenario_id, outcome, date) %>%
  summarise_at(vars(mean:pop_size), sum)
# ------------------------------------------------------------------------------
# Hospital admissions ----------------------------------------------------------
# proportion of hosp admissions per day (per 100,000 people in each 
# age group) by vac strategy
fig_hosp <- ggplot(data = all_res_for_plot_sum %>%
                  filter(
                    #age_group2 == "10-19 years",
                    date >= as.Date("2021-11-01"),
                    #date <= as.Date("2021-12-31"),
                    outcome %in% c("Hospital Admissions")), #, "IC Admissions", "Daily Deaths"
                aes(x = date, y = (mean/pop_size)*100000, fill = scenario_id, 
                    linetype = scenario_id)) +
  geom_ribbon(aes(ymin = (q025/pop_size)*100000, ymax = (q975/pop_size)*100000, 
                  fill = scenario_id), alpha = 0.1) +
  geom_line(aes(color = scenario_id), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(y = "Cases per 100,000 people", x = "Date of infection") +
  #ylim(0,NA) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA), labels = my_round) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.box="vertical",
        axis.title=element_text(size=14, face="bold")) +
  facet_wrap(~age_group2, scales = "free_y") +
  guides(fill=guide_legend("Scenario"), 
         colour = guide_legend("Scenario"),
         linetype = guide_legend("Scenario"))
fig_hosp

# IC admissions ----------------------------------------------------------
# proportion of IC admissions per day (per 100,000 people in each 
# age group) by vac strategy
fig_ic <- ggplot(data = all_res_for_plot_sum %>%
                     filter(
                       #age_group2 == "10-19 years",
                       date >= as.Date("2021-11-01"),
                       #date <= as.Date("2021-12-31"),
                       outcome %in% c("IC Admissions")), #, "IC Admissions", "Daily Deaths"
                   aes(x = date, y = (mean/pop_size)*100000, fill = scenario_id, 
                       linetype = scenario_id)) +
  geom_ribbon(aes(ymin = (q025/pop_size)*100000, ymax = (q975/pop_size)*100000, 
                  fill = scenario_id), alpha = 0.1) +
  geom_line(aes(color = scenario_id), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(y = "Cases per 100,000 people", x = "Date of infection") +
  #ylim(0,NA) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA), labels = my_round) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.box="vertical",
        axis.title=element_text(size=14, face="bold")) +
  facet_wrap(~age_group2, scales = "free_y") +
  guides(fill=guide_legend("Scenario"), 
         colour = guide_legend("Scenario"),
         linetype = guide_legend("Scenario"))
fig_ic

# Daily Deaths ----------------------------------------------------------
# proportion of deaths per day (per 100,000 people in each 
# age group) by vac strategy
my_round <- function(x){round(x, 3)}

fig_death <- ggplot(data = all_res_for_plot_sum %>%
                   filter(
                     #age_group2 == "10-19 years",
                     date >= as.Date("2021-11-01"),
                     #date <= as.Date("2021-12-31"),
                     outcome %in% c("Daily Deaths")), #, "IC Admissions", "Daily Deaths"
                 aes(x = date, y = (mean/pop_size)*100000, fill = scenario_id, 
                     linetype = scenario_id)) +
  geom_ribbon(aes(ymin = (q025/pop_size)*100000, ymax = (q975/pop_size)*100000, 
                  fill = scenario_id), alpha = 0.1) +
  geom_line(aes(color = scenario_id), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(y = "Cases per 100,000 people", x = "Date of infection") +
  #ylim(0,NA) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA), labels = my_round) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.box="vertical",
        axis.title=element_text(size=14, face="bold")) +
  facet_wrap(~age_group2, scales = "free_y") +
  guides(fill=guide_legend("Scenario"), 
         colour = guide_legend("Scenario"),
         linetype = guide_legend("Scenario"))
fig_death

# create one plot from all three plots
fig_sev <- plot_grid(fig_hosp + theme(legend.position = "none",
                                      axis.text.x = element_blank(),
                                      axis.title.x = element_blank()), 
                     fig_ic + theme(legend.position = "none",
                                    axis.text.x = element_blank(),
                                    axis.title.x = element_blank()),
                     fig_death + theme(legend.position = "none"), 
                     labels = "AUTO", nrow = 3, rel_heights = c(1,1,1.5))

legend <- get_legend(
  fig_hosp + theme(legend.box.margin = margin(0, 0, 0, 12))
)

fig_severe <- plot_grid(fig_sev, legend, rel_heights = c(3,.2), nrow = 2)
fig_severe

ggsave(filename = "/rivm/s/ainsliek/results/impact_vac/resubmission/figure_severe_disease.pdf", 
       plot = fig_severe, units = "in", height = 12, width = 10, dpi = 300)
