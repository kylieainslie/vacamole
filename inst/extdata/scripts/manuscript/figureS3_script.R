# ----------------------------------------------------------------
# Vaccination coverage figure (for supplemental material)
# ----------------------------------------------------------------

# read in vac schedules --------------------------------------------
basis_5plus <- read_csv("inst/extdata/inputs/vac_schedule_5plus.csv") %>%
  select(-starts_with("X")) %>%
  mutate(Scenario = "Vaccination of 5+")

basis_12plus <- read_csv("inst/extdata/inputs/vac_schedule_12plus.csv") %>%
  select(-starts_with("X")) %>%
  mutate(Scenario = "Vaccination of 12+",
         date = as.Date(date, format = "%m/%d/%Y"))

basis_18plus <- read_csv("inst/extdata/inputs/vac_schedule_18plus.csv") %>%
  select(-starts_with("X")) %>%
  mutate(Scenario = "Vaccination of 18+",
         date = as.Date(date, format = "%m/%d/%Y"))

# data wrangling -------------------------------------------------
age_dist_10 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
                 0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
n <- 17407585 # Dutch population size
n_vec_10 <- n * age_dist_10

vac_sched_comb <- bind_rows(basis_5plus, basis_12plus, basis_18plus) %>%
  mutate(Scenario = factor(Scenario, levels = c("Vaccination of 5+", "Vaccination of 12+", "Vaccination of 18+")))

vac_sched_long <- vac_sched_comb %>%
  group_by(Scenario) %>%
  mutate(pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d1_9 = (ja_d1_9 * n_vec_10[9] + ja_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d2_9 = (ja_d2_9 * n_vec_10[9] + ja_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
  ) %>%
  select(Scenario, date, pf_d1_1:ja_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10, 
         -az_d2_10, -ja_d1_10, -ja_d2_10) %>%
  pivot_longer(pf_d1_1:ja_d2_9, 
               names_to = c("vaccine", "dose", "age_group"),
               names_sep = "_",
               values_to = "coverage") %>%
  mutate(dose = as.factor(case_when(
    dose == "d1" ~ "Dose 1",
    dose == "d2" ~ "Dose 2"
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
fig_18plus <- ggplot(data = vac_sched_long %>%
                       filter(Scenario == "Vaccination of 18+",
                              age_group %in% c("20-29", "30-39", "40-49", "50-59",
                                               "60-69", "70-79", "80+")
                              ), 
                     aes(x = date, y = coverage, fill = Vaccine, color = Vaccine)) +
  geom_bar(stat = "identity") +
  labs(y = "Vaccine Coverage", x = "Date") +
  ylim(0,1) +
  scale_x_date(date_breaks = "3 months", date_labels = "%d %b %Y") +
  facet_grid(dose~age_group) +
  theme(legend.position = "bottom",
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    #axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.title=element_text(size=14,face="bold"))
fig_18plus

# all vac scenarios (only 0-9 and 10-19) ---------------------------
fig_5plus <- ggplot(data = vac_sched_long %>%
                       filter(Scenario == "Vaccination of 5+",
                              age_group %in% c("0-9", "10-19")), 
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
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 14),
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        axis.title=element_text(size=12,face="bold"))
fig_5plus

fig_12plus <- ggplot(data = vac_sched_long %>%
                      filter(Scenario == "Vaccination of 12+",
                             age_group %in% c("0-9", "10-19")), 
                    aes(x = date, y = coverage, fill = Vaccine, color = Vaccine)) +
  geom_bar(stat = "identity") +
  labs(y = "Vaccine Coverage", x = "Date") +
  ylim(0,1) +
  scale_x_date(date_breaks = "3 months", date_labels = "%d %b %Y") +
  facet_grid(dose~age_group) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1, size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 14),
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        axis.title=element_text(size=12,face="bold"))
fig_12plus

fig_18plus1 <- ggplot(data = vac_sched_long %>%
                                 filter(Scenario == "Vaccination of 18+",
                                        age_group %in% c("0-9", "10-19")), 
                               aes(x = date, y = coverage, fill = Vaccine, color = Vaccine)) +
  geom_bar(stat = "identity") +
  labs(y = "Vaccine Coverage", x = "Date") +
  ylim(0,1) +
  scale_x_date(date_breaks = "3 months", date_labels = "%d %b %Y") +
  facet_grid(dose~age_group) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        #axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 14),
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        axis.title=element_text(size=12,face="bold"))
fig_18plus1

# combine the plots -----------------------------------------------
# first for 0-9 and 10-19 age groups

fig_ext_no_legend <- plot_grid(fig_5plus + theme(legend.position = "none"), 
                               fig_12plus + theme(legend.position = "none"),
                               fig_18plus1 + theme(legend.position = "none"),
                               rel_heights = c(0.6,0.6, 1),
                               labels = "AUTO", nrow = 3)

fig_all_no_legend <- plot_grid(fig_ext_no_legend, fig_18plus + theme(legend.position = "none"),
                               rel_widths = c(0.6, 1),
                               labels = c("", "D"), nrow = 1)
legend <- get_legend(
  fig_18plus + theme(legend.box.margin = margin(12, 0, 0, 0))
)

fig_vc <- plot_grid(fig_all_no_legend, legend, rel_heights = c(3, .4), nrow = 2)
#fig_vc

ggsave(filename = "/rivm/s/ainsliek/results/impact_vac/vac_coverage_bar_plot.jpg", plot = fig_vc,
       units = "in", height = 8, width = 13, dpi = 300)
