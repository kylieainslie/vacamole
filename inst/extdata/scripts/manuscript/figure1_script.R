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

# read in results df -----------------------------------------------------------
df_all <- readRDS("inst/extdata/results/impact_vac/resubmission/results_all.rds")

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
    )
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
  select(date, scenario_id, age_group2, outcome, mean:q975) 

# Figure 1 (main) -------------------------------------------
fig1a <- ggplot(data = all_res_for_plot_sum %>%
                  filter(#scenario_id == "Vaccination in 12+",
                         #age_group2 == "10-19 years",
                         outcome == "Daily Cases"), 
                aes(x = date, y = mean, fill = age_group2, linetype = scenario_id)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = age_group2), alpha = 0.3) +
  geom_line(aes(color = age_group2), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(y = "Cases per day", x = "Date of infection") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.box="vertical",
        axis.title=element_text(size=14,face="bold")) +
  guides(fill=guide_legend("Age Group"), colour = guide_legend("Age Group"),
         linetype = guide_legend("Strategy"))
fig1a

# Figure 1 (inset) ------------------------------------------
# bar plot of percent difference ----------------------------
fig1_inset <- ggplot(data = table1a %>%
                     filter(Scenario != "Vaccination of 18+"), 
                     aes(x = outcome, y = abs(perc_diff), fill = age_group2)) +
  geom_bar(stat = "Identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = abs(perc_diff_upper), ymax = abs(perc_diff_lower), width = 0.2),
                position = position_dodge(0.9)) +
  labs(x = "", y = "Percent Reduction (%)", fill = "Age Group") +
  scale_x_discrete(labels = c("Daily\nCases", "Hospital\nAdmissions", "IC\n Admissions")) +
  facet_wrap(~Scenario, nrow = 2) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12), #angle = 45, hjust = 1, 
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face = "bold")) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
fig1_inset

fig1 <- fig1a + annotation_custom(ggplotGrob(fig1_inset), 
                                  xmin = as.Date("2022-02-01"), xmax = as.Date("2022-04-04"),
                                  ymin = 15000, ymax = 100000)
fig1

# save output -------------------------------------------------
# ggsave(filename = "/rivm/s/ainsliek/results/impact_vac/figure_1_w_inset.pdf", plot = fig1,
#        units = "in", height = 10, width = 12, dpi = 300)
