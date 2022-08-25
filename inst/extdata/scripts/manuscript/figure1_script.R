# -----------------------------------------------------------
# Figure 1 script
# simulated outcomes in whole pop w/o waning 12+ vs. 18+
# -----------------------------------------------------------

# source files ----------------------------------------------
#source("inst/extdata/scripts/manuscript/data_wrangling_for_figures.R")
#source("inst/extdata/scripts/manuscript/tables_script.R")

# load required packages -------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

# read in results df -----------------------------------------------------------
# df_all <- readRDS("/rivm/s/ainsliek/results/impact_vac/resubmission/results_all.rds")
df_all <- readRDS("S:/R/ainsliek/results/impact_vac/resubmission/results_all.rds")

# population size for determining outcomes per 100,000 people ------------------
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 
              0.13083463,0.14514332, 0.12092904, 0.08807406, 
              0.04622194)
n <- 17407585 # Dutch population size
n_vec <- n * age_dist

# data wrangling for plots -----------------------------------------------------
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
df_fig1 <- all_res_for_plot %>%
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

# total 
df_total <- df_fig1 %>%
  ungroup() %>%
  group_by(scenario_id, outcome, date) %>%
  summarise_at(vars(mean:pop_size), sum)
# ------------------------------------------------------------------------------

# Figure 1 ---------------------------------------------------------------------
# Plot outcomes per day (per 100,000 people in each age group) by vac strategy
my_round <- function(x){round(x, 3)}

# infections -------------------------------------------------------------------
fig1_inf <- ggplot(data = df_fig1 %>%
                  filter(#scenario_id == "Vaccination in 5+",
                    #age_group2 == "10-19 years",
                    #date >= as.Date("2021-11-01"),
                    #date <= as.Date("2021-12-31"),
                    outcome == "Daily Infections"), 
                aes(x = date, y = (mean/pop_size)*100000, fill = scenario_id, linetype = scenario_id)) +
  geom_ribbon(aes(ymin = (q025/pop_size)*100000, ymax = (q975/pop_size)*100000, fill = scenario_id), alpha = 0.1) +
  geom_line(aes(color = scenario_id), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(y = "Infections per 100,000 people", x = "Date") +
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
  facet_wrap(~age_group2, nrow = 3) +
  guides(fill=guide_legend("Scenario"), 
         colour = guide_legend("Scenario"),
         linetype = guide_legend("Scenario"))
#fig1_inf
# ------------------------------------------------------------------------------

# hospital admissions ----------------------------------------------------------
fig1_hosp <- ggplot(data = df_fig1 %>%
                     filter(#scenario_id == "Vaccination in 5+",
                       #age_group2 == "10-19 years",
                       #date >= as.Date("2021-11-01"),
                       #date <= as.Date("2021-12-31"),
                       outcome == "Hospital Admissions"), 
                   aes(x = date, y = (mean/pop_size)*100000, fill = scenario_id, linetype = scenario_id)) +
  geom_ribbon(aes(ymin = (q025/pop_size)*100000, ymax = (q975/pop_size)*100000, fill = scenario_id), alpha = 0.1) +
  geom_line(aes(color = scenario_id), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(y = "Hospital Admissions per 100,000 people", x = "Date") +
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
  facet_wrap(~age_group2, nrow = 3, scales = "free_y"
             ) +
  guides(fill=guide_legend("Scenario"), 
         colour = guide_legend("Scenario"),
         linetype = guide_legend("Scenario"))
#fig1_hosp
# ------------------------------------------------------------------------------

# ICU admissions ---------------------------------------------------------------
fig1_icu <- ggplot(data = df_fig1 %>%
                      filter(#scenario_id == "Vaccination in 5+",
                        #age_group2 == "10-19 years",
                        #date >= as.Date("2021-11-01"),
                        #date <= as.Date("2021-12-31"),
                        outcome == "IC Admissions"), 
                    aes(x = date, y = (mean/pop_size)*100000, fill = scenario_id, linetype = scenario_id)) +
  geom_ribbon(aes(ymin = (q025/pop_size)*100000, ymax = (q975/pop_size)*100000, fill = scenario_id), alpha = 0.1) +
  geom_line(aes(color = scenario_id), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(y = "IC Admissions per 100,000 people", x = "Date") +
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
  facet_wrap(~age_group2, nrow = 3, scales = "free_y"
             ) +
  guides(fill=guide_legend("Scenario"), 
         colour = guide_legend("Scenario"),
         linetype = guide_legend("Scenario"))
#fig1_icu
# ------------------------------------------------------------------------------

# deaths -----------------------------------------------------------------------
fig1_death <- ggplot(data = df_fig1 %>%
                      filter(#scenario_id == "Vaccination in 5+",
                        #age_group2 == "10-19 years",
                        #date >= as.Date("2021-11-01"),
                        #date <= as.Date("2021-12-31"),
                        outcome == "Daily Deaths"), 
                    aes(x = date, y = (mean/pop_size)*100000, fill = scenario_id, 
                        linetype = scenario_id)) +
  geom_ribbon(aes(ymin = (q025/pop_size)*100000, ymax = (q975/pop_size)*100000, 
                  fill = scenario_id), alpha = 0.1) +
  geom_line(aes(color = scenario_id), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(y = "Deaths per 100,000 people", x = "Date") +
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
  facet_wrap(~age_group2, nrow = 3, scales = "free_y"
             ) +
  guides(fill=guide_legend("Scenario"), 
         colour = guide_legend("Scenario"),
         linetype = guide_legend("Scenario"))
#fig1_death
# ------------------------------------------------------------------------------

# combine plots ----------------------------------------------------------------
fig1 <- plot_grid(fig1_inf   + theme(legend.position = "none"), 
                  fig1_hosp  + theme(legend.position = "none"),
                  fig1_icu   + theme(legend.position = "none"),
                  fig1_death + theme(legend.position = "none"), 
                  labels = "AUTO", nrow = 1#, rel_heights = c(1,1,1.5)
                  )

legend <- get_legend(
  fig1_inf + theme(legend.box.margin = margin(0, 0, 0, 12))
)

figure1 <- plot_grid(fig1, legend, rel_heights = c(3,.2), nrow = 2)
figure1

# ggsave(filename = "/rivm/s/ainsliek/results/impact_vac/resubmission/figure1_per100k.pdf",
#        plot = figure1, units = "in", height = 10, width = 26, dpi = 300)

ggsave(filename = "S:/R/ainsliek/results/impact_vac/resubmission/figure1_per100k.pdf",
       plot = figure1, units = "in", height = 10, width = 26, dpi = 300)
