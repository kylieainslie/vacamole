# Plot results of different AstraZeneca vaccine allocation strategies--------------
# 1) old to young
# 2) young to old
# 3) alternative
# 4) no vaccination (control)

# load results --------------------------------------------------------------------
# constant contact matrix
o2y <- readRDS("inst/extdata/results/res_old_to_young.rds") 
y2o <- readRDS("inst/extdata/results/res_young_to_old.rds")
alt <- readRDS("inst/extdata/results/res_alternative.rds")
nov <- readRDS("inst/extdata/results/res_no_vac.rds")

# changing contact matrix
o2y_ccm <- readRDS("inst/extdata/results/res_ccm_old_to_young.rds") 
y2o_ccm <- readRDS("inst/extdata/results/res_ccm_young_to_old.rds")
alt_ccm <- readRDS("inst/extdata/results/res_ccm_alternative.rds")
nov_ccm <- readRDS("inst/extdata/results/res_ccm_no_vac.rds")

# changing contact matrix
o2y_1cp <- readRDS("inst/extdata/results/res_1cp_old_to_young.rds") 
y2o_1cp <- readRDS("inst/extdata/results/res_1cp_young_to_old.rds")
alt_1cp <- readRDS("inst/extdata/results/res_1cp_alternative.rds")
nov_1cp <- readRDS("inst/extdata/results/res_1cp_no_vac.rds")

# data wrangle for plotting/summary table -----------------------------------------
results <- bind_rows(o2y, y2o, alt, nov)
results2 <- bind_rows(o2y_ccm, y2o_ccm, alt_ccm, nov_ccm)

# summary table -------------------------------------------------------------------
summary_tab <- results %>%
  group_by(scenario, outcome) %>%
  summarise_at(.vars = "value", .funs = sum)

inc_sum <- summary_tab %>%
  filter(outcome == "incidence") %>%
  mutate(perc_diff = ((value*100)/2211084) - 100)

hosp_sum <- summary_tab %>%
  filter(outcome == "hosp_admissions") %>%
  mutate(perc_diff = ((value*100)/18898) - 100)

summary_tab2 <- results2 %>%
  group_by(scenario, outcome) %>%
  summarise_at(.vars = "value", .funs = sum)

inc_sum <- summary_tab2 %>%
  filter(outcome == "incidence") %>%
  mutate(perc_diff = ((value*100)/2756636) - 100)

hosp_sum <- summary_tab2 %>%
  filter(outcome == "hosp_admissions") %>%
  mutate(perc_diff = ((value*100)/21511) - 100)


# plot ----------------------------------------------------------------------------
# constant contact matrix
g_sum <- ggplot(results %>% filter(time > 3), aes(x = date, y = value, color = scenario)) +
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

ggsave("inst/extdata/results/vac_alloc_results_lockdown_cm_only.jpg",
       plot = g_sum,
       height = 6,
       width = 12,
       dpi = 300)

# changing contact matrix
g_sum2 <- ggplot(results2 %>% filter(time > 3), aes(x = date, y = value, color = scenario)) +
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

ggsave("inst/extdata/results/vac_alloc_results_lockdown_cm_only.jpg",
       plot = g_sum2,
       height = 6,
       width = 12,
       dpi = 300)
