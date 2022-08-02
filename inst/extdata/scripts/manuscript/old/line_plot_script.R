# -----------------------------------------------------------
# Figure S1 script
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
df_all <- readRDS("/rivm/s/ainsliek/results/impact_vac/resubmission/results_all.rds")

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

# plot -------------------------------------------------------------------------
# mean line with ribbon 
figS1 <- ggplot(data = all_res_for_plot %>%
                     filter(outcome %in% c("Daily Infections")#,
                            #date < as.Date("2022-10-01")
                     ), # "inc hosp", "inc icu", "inc death"
                   aes(x = date, y = value, color = scenario_id, group = sample)) +
  geom_line() +
  xlab("Date") + 
  ylab("Mean value") +
  scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.box="vertical",
        axis.title=element_text(size=14,face="bold")) +
  # guides(#fill=guide_legend("Scenario ID"), 
  #        colour = guide_legend("Scenario ID")) +
  facet_grid(.~age_group2)
  # annotate("rect", xmin = as.Date("2022-09-22"), xmax = as.Date("2022-12-15"), ymin = 0, ymax = 200000, 
  #          alpha = .5)
figS1
ggsave(filename = "/rivm/s/ainsliek/results/scenario_hub/round2/case_plot_round2.jpg", 
       plot = p_ribbon,
       units = "in", height = 8, width = 13, dpi = 300)
