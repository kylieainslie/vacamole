# ----------------------------------------------------------------
# Plot Vaccination coverage to check vac schedules
# ----------------------------------------------------------------

# read in vac schedules --------------------------------------------
vac_scheduleAC <- readRDS("inst/extdata/inputs/vaccination_schedules/vac_schedule_scenario_hub_round2_AC.rds") 
vac_scheduleBD <- readRDS("inst/extdata/inputs/vaccination_schedules/vac_schedule_scenario_hub_round2_BD.rds") 

# data wrangling -------------------------------------------------
age_dist_10 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
                 0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
n <- 17407585 # Dutch population size
n_vec_10 <- n * age_dist_10

vac_sched_comb <- bind_rows(vac_scheduleAC, vac_scheduleBD, .id = "Scenario") %>%
  mutate(Scenario = case_when(
    Scenario == 1 ~ "Booster campaign 60+",
    Scenario == 2 ~ "Booster campaign 18+"))

vac_sched_long <- vac_sched_comb %>%
  group_by(Scenario) %>%
  mutate(pf_d1_9 = (.data$pf_d1_9 * n_vec_10[9] + .data$pf_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         pf_d2_9 = (.data$pf_d2_9 * n_vec_10[9] + .data$pf_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         pf_d3_9 = (.data$pf_d3_9 * n_vec_10[9] + .data$pf_d3_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         pf_d4_9 = (.data$pf_d4_9 * n_vec_10[9] + .data$pf_d4_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         pf_d5_9 = (.data$pf_d5_9 * n_vec_10[9] + .data$pf_d5_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         mo_d1_9 = (.data$mo_d1_9 * n_vec_10[9] + .data$mo_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         mo_d2_9 = (.data$mo_d2_9 * n_vec_10[9] + .data$mo_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         mo_d3_9 = (.data$mo_d3_9 * n_vec_10[9] + .data$mo_d3_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         mo_d4_9 = (.data$mo_d4_9 * n_vec_10[9] + .data$mo_d4_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         mo_d5_9 = (.data$mo_d5_9 * n_vec_10[9] + .data$mo_d5_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         az_d1_9 = (.data$az_d1_9 * n_vec_10[9] + .data$az_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         az_d2_9 = (.data$az_d2_9 * n_vec_10[9] + .data$az_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         ja_d1_9 = (.data$ja_d1_9 * n_vec_10[9] + .data$ja_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
         ja_d2_9 = (.data$ja_d2_9 * n_vec_10[9] + .data$ja_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10])
  ) %>%
  select(.data$Scenario, .data$date, .data$pf_d1_1:.data$ja_d2_9, 
         -.data$pf_d1_10, -.data$pf_d2_10, -.data$pf_d3_10, -.data$pf_d4_10, -.data$pf_d5_10,
         -.data$mo_d1_10, -.data$mo_d2_10, -.data$mo_d3_10, -.data$mo_d4_10, -.data$mo_d5_10,
         -.data$az_d1_10, -.data$az_d2_10, -.data$ja_d1_10, -.data$ja_d2_10) %>%
  pivot_longer(pf_d1_1:ja_d2_9, 
               names_to = c("vaccine", "dose", "age_group"),
               names_sep = "_",
               values_to = "coverage") %>%
  mutate(dose = as.factor(case_when(
    dose == "d1" ~ "Dose 1",
    dose == "d2" ~ "Dose 2",
    dose == "d3" ~ "Dose 3",
    dose == "d4" ~ "Dose 4",
    dose == "d5" ~ "Dose 5"
  )),
  coverage = ifelse(coverage == 0, NA, coverage),
  Vaccine = case_when(
    vaccine == "az" ~ "AstraZeneca",
    vaccine == "pf" ~ "Pfizer/BioNTech",
    vaccine == "mo" ~ "Moderna",
    vaccine == "ja" ~ "Janssen"
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
  )) %>%
  filter(date > as.Date("2020-12-31"))

# plot ------------------------------------------------------------
# 18+ -------------------------------------------------------------
fig_AC <- ggplot(data = vac_sched_long %>%
                       filter(Scenario == "Booster campaign 60+" ,
                              age_group %in% c("10-19", "20-29", "30-39", "40-49", "50-59",
                                               "60-69", "70-79", "80+"),
                              dose %in% c("Dose 3", "Dose 4", "Dose 5")
                       ), 
                     aes(x = date, y = coverage, fill = Vaccine, color = Vaccine)) +
  geom_bar(stat = "identity") +
  labs(y = "Vaccine Coverage", x = "Date") +
  ylim(0,1) +
  scale_x_date(date_breaks = "3 months", date_labels = "%d %b %Y") +
  facet_grid(dose~age_group) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_blank(), #element_text(angle = 45, hjust = 1, size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) #+
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70")
fig_AC

# all vac scenarios (only 0-9 and 10-19) ---------------------------------------
fig_BD <- ggplot(data = vac_sched_long %>%
                      filter(Scenario == "Booster campaign 18+",
                             age_group %in% c("10-19", "20-29", "30-39", "40-49", "50-59",
                                              "60-69", "70-79", "80+"),
                             dose %in% c("Dose 3", "Dose 4", "Dose 5")
                             ), 
                    aes(x = date, y = coverage, fill = Vaccine, color = Vaccine)) +
  geom_bar(stat = "identity") +
  labs(y = "Vaccine Coverage", x = "Date") +
  ylim(0,1) +
  scale_x_date(date_breaks = "3 months", date_labels = "%d %b %Y") +
  facet_grid(dose~age_group) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        #axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) #+
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70")
fig_BD


# combine the plots ------------------------------------------------------------
# first for 0-9 and 10-19 age groups

fig_ext_no_legend <- plot_grid(fig_AC + theme(legend.position = "none"), 
                               fig_BD + theme(legend.position = "none"),
                               rel_heights = c(0.8, 1),
                               labels = "AUTO", nrow = 2)

# fig_all_no_legend <- plot_grid(fig_ext_no_legend, fig_AC + theme(legend.position = "none"),
#                                rel_widths = c(0.6, 1),
#                                labels = c("", "B"), nrow = 1)
legend <- get_legend(
  fig_AC + theme(legend.box.margin = margin(12, 0, 0, 0))
)

fig_all <- plot_grid(fig_ext_no_legend, legend, rel_heights = c(3, .4), nrow = 2)
fig_all

ggsave(filename = "/rivm/s/ainsliek/results/scenario_hub/round2/vac_coverage_plot.jpg", plot = fig_all,
       units = "in", height = 8, width = 13, dpi = 300)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# plot comp_ve
vac_ratesCD <- bind_rows(vac_ratesC, vac_ratesD, .id = "Scenario") %>%
  mutate(Scenario = case_when(
    Scenario == 1 ~ "C",
    Scenario == 2 ~ "D"
  )) %>%
  filter(param == "comp_ve")

fig_cv <- ggplot(data = vac_ratesCD %>%
                   filter(dose %in% c("d3", "d4", "d5"),
                          date > as.Date("2022-07-08"),
                          outcome == "hospitalisation"
                   ), 
                 aes(x = date, y = value, color = Scenario)) +
  geom_line() +
  labs(y = "VE", x = "Date") +
  ylim(0,1) +
  scale_x_date(date_breaks = "3 months", date_labels = "%d %b %Y") +
  facet_grid(dose~age_group) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        #axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) #+
#geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70")
fig_cv
