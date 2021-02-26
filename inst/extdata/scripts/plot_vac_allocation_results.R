# Plot results of different AstraZeneca vaccine allocation strategies--------------
# 1) old to young
# 2) young to old
# 3) alternative
# 4) no vaccination (control)

# load results --------------------------------------------------------------------
# constant contact matrix
# o2y <- readRDS("inst/extdata/results/res_old_to_young.rds") 
# y2o <- readRDS("inst/extdata/results/res_young_to_old.rds")
# alt <- readRDS("inst/extdata/results/res_alternative.rds")
# nov <- readRDS("inst/extdata/results/res_no_vac.rds")

# changing contact matrix
# cases as criteria
o2y_ccm <- readRDS("inst/extdata/results/res_ccm_old_to_young.rds") 
y2o_ccm <- readRDS("inst/extdata/results/res_ccm_young_to_old.rds")
alt_ccm <- readRDS("inst/extdata/results/res_ccm_alternative.rds")
nov_ccm <- readRDS("inst/extdata/results/res_ccm_no_vac.rds")

# ic as criteria
o2y_ccm_ic <- readRDS("inst/extdata/results/res_ccm_old_to_young_w_deaths_ic.rds") 
y2o_ccm_ic <- readRDS("inst/extdata/results/res_ccm_young_to_old_w_deaths_ic.rds")
alt_ccm_ic <- readRDS("inst/extdata/results/res_ccm_alternative_w_deaths_ic.rds")
nov_ccm_ic <- readRDS("inst/extdata/results/res_ccm_no_vac_w_deaths_ic.rds")

# by age group
o2y_ccm_ic_by_ag <- readRDS("inst/extdata/results/res_by_age_group_ccm_old_to_young_w_deaths_ic.rds") %>%
  mutate(scenario = "old_to_young")
y2o_ccm_ic_by_ag <- readRDS("inst/extdata/results/res_by_age_group_ccm_young_to_old_w_deaths_ic.rds") %>%
  mutate(scenario = "young_to_old")
alt_ccm_ic_by_ag <- readRDS("inst/extdata/results/res_by_age_group_ccm_alternative_w_deaths_ic.rds") %>%
  mutate(scenario = "alternative")
nov_ccm_ic_by_ag <- readRDS("inst/extdata/results/res_by_age_group_ccm_no_vac_w_deaths_ic.rds") %>%
  mutate(scenario = "no_vac")

# deferring second dose
orig <- readRDS("inst/extdata/results/res_orig_distribution_schedule_extra_cp.rds") %>%
  mutate(scenario = "original")
no_second_doses <- readRDS("inst/extdata/results/res_no_second_doses_extra_cp.rds") %>%
  mutate(scenario = "second dose deferred")

# data wrangle for plotting/summary table -----------------------------------------
results_cases_thresh <- bind_rows(o2y_ccm, y2o_ccm, alt_ccm, nov_ccm)
results_ic_thresh <- bind_rows(o2y_ccm_ic, y2o_ccm_ic, alt_ccm_ic, nov_ccm_ic)
results_ic_by_age_group <- bind_rows(o2y_ccm_ic_by_ag, y2o_ccm_ic_by_ag, alt_ccm_ic_by_ag, nov_ccm_ic_by_ag)
results_deferral <- bind_rows(orig, no_second_doses)
# calculate life years lost
life_expectancy <- c(77.89, 67.93, 58.08, 48.28, 38.6, 29.22, 20.52, 12.76, 4.35)
deaths_by_ag <- results_ic_by_age_group %>%
  group_by(scenario, outcome, age_group) %>%
  summarise_at(.vars = "value", .funs = sum) %>%
  filter(outcome == "new_deaths")

no_vac_deaths_by_ag <- deaths_by_ag %>%
  filter(scenario == "no_vac")

deaths_prevented <- deaths_by_ag %>%
  mutate(deaths_prevented = no_vac_deaths_by_ag$value - value,
         life_years_lost_prevented = deaths_prevented * life_expectancy,
         life_years_lost = value * life_expectancy)
  
total_life_years_lost <- deaths_prevented %>%
  group_by(scenario) %>%
  summarise_at(.vars = c("value","life_years_lost", "life_years_lost_prevented"), .funs = sum) %>%
  mutate(perc_diff = ((life_years_lost*100)/44505) - 100)

# summary table -------------------------------------------------------------------
results_dat <- results_deferral
# change data frame in summary tab for summary results from ic_thresh
summary_tab <- results_dat %>%
  group_by(scenario, outcome) %>%
  summarise_at(.vars = "value", .funs = sum)

ref_no_vac <- summary_tab %>%
  filter(scenario == "ccm_no_vac_w_deaths_ic") %>%
  group_by(outcome) %>%
  summarise_at(.vars = "value", .funs = sum)

infections_sum <- summary_tab %>%
  filter(outcome == "new_infections") %>%
  mutate(perc_diff = ((value*100)/ref_no_vac$value[8]) - 100)

cases_sum <- summary_tab %>%
  filter(outcome == "new_cases") %>%
  mutate(perc_diff = ((value*100)/ref_no_vac$value[6]) - 100)

hosp_admissions_sum <- summary_tab %>%
  filter(outcome == "hospital_admissions") %>%
  mutate(perc_diff = ((value*100)/ref_no_vac$value[2]) - 100)

ic_admissions_sum <- summary_tab %>%
  filter(outcome == "ic_admissions") %>%
  mutate(perc_diff = ((value*100)/ref_no_vac$value[4]) - 100)

deaths_sum <- summary_tab %>%
  filter(outcome == "new_deaths") %>%
  mutate(perc_diff = ((value*100)/ref_no_vac$value[7]) - 100)

# plot ----------------------------------------------------------------------------
# subset for plotting
for_plot <- results_cases_thresh %>%
  filter(outcome %in% c("new_infections", "new_cases", "hospital_admissions", "ic_admissions")) %>%
  group_by(outcome, scenario) %>%
  mutate(rolling_avg = zoo::rollmean(value, k = 7, fill = NA),
         outcome = case_when(
           outcome == "hospital_admissions" ~ "Hospital Admissions",
           outcome == "ic_admissions" ~ "IC Admissions",
           outcome == "new_cases" ~ "New Cases",
           outcome == "new_infections" ~ "New Infections"
         ),
         outcome = factor(outcome, levels = c("New Infections", "New Cases", "Hospital Admissions", "IC Admissions")),
         scenario = case_when(
           scenario == "ccm_alternative" ~ "Alternative",
           scenario == "ccm_no_vac" ~ "No Vaccination",
           scenario == "ccm_old_to_young" ~ "Old to Young",
           scenario == "ccm_young_to_old" ~ "Young to Old"
         ),
         scenario = factor(scenario, levels = c("Old to Young", "Young to Old", "Alternative", "No Vaccination")))
  
# cases threshhold
g_sum <- ggplot(for_plot, aes(x = date, y = rolling_avg, color = scenario)) +
  geom_line(position=position_dodge(width=2.5)) +
  labs(y = "Value", x = "Time (days)", color = "Vaccination Scenario") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  #geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~outcome, scales = "free")
g_sum

ggsave("inst/extdata/results/vac_alloc_results_cases_thresh.jpg",
       plot = g_sum,
       height = 6,
       width = 12,
       dpi = 300)

# ic threshhold
for_plot2 <- results_ic_thresh %>%
  filter(outcome %in% c("new_infections", "new_cases", "hospital_admissions", "ic_admissions", "new_deaths")) %>%
  group_by(outcome, scenario) %>%
  mutate(rolling_avg = zoo::rollmean(value, k = 7, fill = NA),
         outcome = case_when(
           outcome == "hospital_admissions" ~ "Hospital Admissions",
           outcome == "ic_admissions" ~ "IC Admissions",
           outcome == "new_cases" ~ "New Cases",
           outcome == "new_infections" ~ "New Infections",
           outcome == "new_deaths" ~ "New Deaths"
         ),
         outcome = factor(outcome, levels = c("New Infections", "New Cases", "Hospital Admissions", "IC Admissions", "New Deaths")),
         scenario = case_when(
           scenario == "ccm_alternative_w_deaths_ic" ~ "Alternative",
           scenario == "ccm_no_vac_w_deaths_ic" ~ "No Vaccination",
           scenario == "ccm_old_to_young_w_deaths_ic" ~ "Old to Young",
           scenario == "ccm_young_to_old_w_deaths_ic" ~ "Young to Old"
         ),
         scenario = factor(scenario, levels = c("Old to Young", "Young to Old", "Alternative", "No Vaccination")))

g_sum2 <- ggplot(for_plot2, aes(x = date, y = value, color = scenario)) +
  geom_line(position=position_dodge(width=2.5)) +
  labs(y = "Value", x = "Time (days)", color = "Vaccination Scenario") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  #geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~outcome, scales = "free")
g_sum2

ggsave("inst/extdata/results/vac_alloc_results_ic_w_deaths.jpg",
       plot = g_sum2,
       height = 6,
       width = 12,
       dpi = 300)

# defer second dose
for_plot3 <- results_dat %>%
  filter(outcome %in% c("new_infections", "new_cases", "hospital_admissions", "ic_admissions", "new_deaths")) %>%
  group_by(outcome, scenario) %>%
  mutate(rolling_avg = zoo::rollmean(value, k = 7, fill = NA),
         outcome = case_when(
           outcome == "hospital_admissions" ~ "Hospital Admissions",
           outcome == "ic_admissions" ~ "IC Admissions",
           outcome == "new_cases" ~ "New Cases",
           outcome == "new_infections" ~ "New Infections",
           outcome == "new_deaths" ~ "New Deaths"
         ),
         outcome = factor(outcome, levels = c("New Infections", "New Cases", "Hospital Admissions", "IC Admissions", "New Deaths")),
         scenario = case_when(
           scenario == "original" ~ "Original",
           scenario == "second dose deferred" ~ "Second Dose Deferred"
         ),
         scenario = factor(scenario, levels = c("Original", "Second Dose Deferred")))

g_sum3 <- ggplot(for_plot3, aes(x = date, y = value, color = scenario)) +
  geom_line(position=position_dodge(width=2.5)) +
  labs(y = "Value", x = "Time (days)", color = "Vaccination Scenario") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  #geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~outcome, scales = "free")
g_sum3

ggsave("inst/extdata/results/second_dose_deferred.jpg",
       plot = g_sum3,
       height = 6,
       width = 12,
       dpi = 300)