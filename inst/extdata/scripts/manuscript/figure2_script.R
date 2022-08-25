# ------------------------------------------------------------------------------
# Figure 2 script
# Heatmap of % reduction in disease outcomes by vaccination scenario
# ------------------------------------------------------------------------------
# library(ggsci)
# library(gridExtra)
library(tidyr)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(cowplot)
# -------------------------------------------------------------

# read in results df -----------------------------------------------------------
#df_all <- readRDS("/rivm/s/ainsliek/results/impact_vac/resubmission/results_all.rds")
df_all <- readRDS("S:/R/ainsliek/results/impact_vac/resubmission/results_all.rds")

# population size for determining outcomes per 100,000 people ------------------
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 
              0.13083463,0.14514332, 0.12092904, 0.08807406, 
              0.04622194)
n <- 17407585 # Dutch population size
n_vec <- n * age_dist

# data wrangling for plots -----------------------------------------------------
for_plot2 <- df_all %>%
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
    pop_size = case_when(
      age_group == 1 ~ n_vec[1],
      age_group == 2 ~ n_vec[2],
      age_group == 3 ~ n_vec[3],
      age_group == 4 ~ n_vec[4],
      age_group == 5 ~ n_vec[5],
      age_group == 6 ~ n_vec[6],
      age_group == 7 ~ n_vec[7],
      age_group == 8 ~ n_vec[8],
      age_group == 9 ~ n_vec[9]
    ),
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
  ) %>%
  group_by(date, scenario_id, target_variable, age_group, sample) %>%
  mutate(value_per100k = (value/pop_size)*100000)


# calculate absolute and percent difference in disease outcomes (per 100k) -----
table1 <- for_plot2 %>%
  # sum disease outcomes over days
  group_by(scenario_id, outcome, age_group, sample) %>%
  summarise(value_cum = sum(value_per100k)) %>%
  ungroup() %>%
  # determine percent difference for each sample
  group_by(outcome, age_group, sample) %>%
  mutate(
    abs_diff  = value_cum - value_cum[scenario_id == 'Vaccination in 18+'],
    abs_diff = ifelse(is.nan(abs_diff), NA, abs_diff),
    perc_diff = (value_cum - value_cum[scenario_id == 'Vaccination in 18+'])/value_cum[scenario_id == 'Vaccination in 18+'],
    perc_diff = ifelse(is.nan(perc_diff), NA, perc_diff)
    ) %>%
  # calculate summary statistics
  group_by(scenario_id, outcome, age_group) %>%
  summarise(
    abs_diff_mean  = mean(abs_diff, na.rm = TRUE),
    abs_diff_q025 = quantile(abs_diff, probs = 0.025, na.rm = TRUE),
    abs_diff_q975 = quantile(abs_diff, probs = 0.975, na.rm = TRUE),
    perc_diff_mean = mean(perc_diff, na.rm = TRUE),
    perc_diff_q025 = quantile(perc_diff, probs = 0.025, na.rm = TRUE),
    perc_diff_q975 = quantile(perc_diff, probs = 0.975, na.rm = TRUE)#,
    #perc_diff_plot = ifelse(is.nan(abs(perc_diff_mean)*10), NA, abs(perc_diff_mean)*10) #tranforming it so it can be on log
    # scale when plotting the colors
  ) %>%
  ungroup()
table1

# make a tile plot (heatmap) of percent difference by outcome and age group ----
only_12plus <- table1 %>%
  filter(
    outcome != "Daily Cases",
    scenario_id == "Vaccination in 12+"
  )
fig2a <- ggplot(data = only_12plus, 
               aes(x = age_group, y = outcome, fill= perc_diff_mean)) + 
  geom_tile() + 
  scale_fill_viridis_c(
    limits = c(-0.35, 0),
    breaks = c(-0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0),
    labels = c(-0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0),
    #trans = "log",
    direction = 1) +
  theme_light() +
  guides(fill = guide_legend("Percent Difference")) +
  labs(y = "Outcome", x = "Age Group")
fig2a

fig2b <- ggplot(data = table1 %>%
                  filter(
                    outcome != "Daily Cases",
                    scenario_id == "Vaccination in 5+"
                  ), 
                aes(x = age_group, y = outcome, fill= perc_diff_mean)) + 
  geom_tile() + 
  scale_fill_viridis_c(
    limits = c(-0.35, 0),
    breaks = c(-0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0),
    labels = c(-0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0),
    #trans = "log",
    direction = 1) +
  theme_light() +
  guides(fill = guide_legend("Percent Difference")) +
  labs(y = "Outcome", x = "Age Group")
fig2b

fig2 <- plot_grid(fig2a + theme(legend.position = "none"),
                     fig2b + theme(legend.position = "none"), 
                     labels = "AUTO")

legend <- get_legend(
  fig2a + theme(legend.box.margin = margin(t= 0, r= 20, b = 0, l = 0))
)

figure2 <- plot_grid(fig2, legend, rel_widths = c(5,1), nrow = 1)
figure2


# ggsave(filename = "/rivm/s/ainsliek/results/impact_vac/resubmission/figure2.pdf",
#        plot = figure2, units = "in", height = 6, width = 12, dpi = 300)

ggsave(filename = "S:/R/ainsliek/results/impact_vac/resubmission/figure2.pdf",
       plot = figure2, units = "in", height = 6, width = 12, dpi = 300)
