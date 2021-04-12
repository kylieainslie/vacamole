# Plot results of different AstraZeneca vaccine allocation strategies--------------
# 1) old to young
# 2) young to old
# 3) alternative
# 4) no vaccination (control)

# load results --------------------------------------------------------------------
# main analysis
o2y <- readRDS("inst/extdata/results/res_o2y_15March.rds") 
y2o <- readRDS("inst/extdata/results/res_y2o_15March.rds") 
alt <- readRDS("inst/extdata/results/res_alt_15March.rds") 
nov_healthy <- readRDS("inst/extdata/results/res_no_vac_healthy_15March.rds") 
nov <- readRDS("inst/extdata/results/res_no_vac_15March.rds")

# sensitivity analysis
o2y_sa <- readRDS("inst/extdata/results/res_o2y_re-lockdown_15March.rds") 
y2o_sa <- readRDS("inst/extdata/results/res_y2o_re-lockdown_15March.rds") 
alt_sa <- readRDS("inst/extdata/results/res_alt_re-lockdown_15March.rds") 
nov_healthy_sa <- readRDS("inst/extdata/results/res_no_vac_healthy_re-lockdown_15March.rds") 
nov_sa <- readRDS("inst/extdata/results/res_no_vac_re-lockdown_15March.rds")

# 20%
o2y_sa20 <- readRDS("inst/extdata/results/res_o2y_sa20_14March.rds") 
y2o_sa20 <- readRDS("inst/extdata/results/res_y2o_sa20_14March.rds") 
alt_sa20 <- readRDS("inst/extdata/results/res_alt_sa20_14March.rds") 
nov_healthy_sa20 <- readRDS("inst/extdata/results/res_no_vac_healthy_sa20_14March.rds") 
nov_sa20 <- readRDS("inst/extdata/results/res_no_vac_sa20_14March.rds")

# 50%
o2y_sa50 <- readRDS("inst/extdata/results/res_o2y_sa50_14March.rds") 
y2o_sa50 <- readRDS("inst/extdata/results/res_y2o_sa50_14March.rds") 
alt_sa50 <- readRDS("inst/extdata/results/res_alt_sa50_14March.rds") 
nov_healthy_sa50 <- readRDS("inst/extdata/results/res_no_vac_healthy_sa50_14March.rds") 
nov_sa50 <- readRDS("inst/extdata/results/res_no_vac_sa50_14March.rds")

# defer second dose
basis_res <- readRDS("inst/extdata/results/res_basis_no_interventions_02April.rds") 
defer_res <- readRDS("inst/extdata/results/res_defer_no_interventions_02April.rds") 
# with Vasileiou et al. 2021 estimates
basis_vasileiou_res <- readRDS("inst/extdata/results/res_basis_vasileiou_no_interventions_02April.rds") 
defer_vasileiou_res <- readRDS("inst/extdata/results/res_defer_vasileiou_no_interventions_02April.rds") 

# data wrangle for plotting/summary table -----------------------------------------
results_main <- bind_rows(o2y$df_summary,
                     y2o$df_summary,
                     alt$df_summary,
                     nov_healthy$df_summary,
                     nov$df_summary,
                     .id = "scenario") %>%
  mutate(scenario = case_when(
    scenario == 1 ~ "old to young",
    scenario == 2 ~ "young to old",
    scenario == 3 ~ "alternative",
    scenario == 4 ~ "no vaccination (healthy adults)",
    scenario == 5 ~ "no vaccination (at all)"
  ),
  analysis = "main")

results_sa <- bind_rows(o2y_sa$df_summary,
                     y2o_sa$df_summary,
                     alt_sa$df_summary,
                     nov_healthy_sa$df_summary,
                     nov_sa$df_summary,
                     .id = "scenario") %>%
  mutate(scenario = case_when(
    scenario == 1 ~ "old to young",
    scenario == 2 ~ "young to old",
    scenario == 3 ~ "alternative",
    scenario == 4 ~ "no vaccination (healthy adults)",
    scenario == 5 ~ "no vaccination (at all)"
  ),
  analysis = "sa20")

results_sa50 <- bind_rows(o2y_sa50$df_summary,
                          y2o_sa50$df_summary,
                          alt_sa50$df_summary,
                          nov_healthy_sa50$df_summary,
                          .id = "scenario") %>%
  mutate(scenario = case_when(
    scenario == 1 ~ "old to young",
    scenario == 2 ~ "young to old",
    scenario == 3 ~ "alternative",
    scenario == 4 ~ "no vaccination (healthy adults)"
  ),
  analysis = "sa50")

results_defer <- bind_rows(basis_res$df_summary,
                          defer_res$df_summary,
                          basis_vasileiou_res$df_summary,
                          defer_vasileiou_res$df_summary,
                          .id = "scenario") %>%
  mutate(scenario = case_when(
    scenario == 1 ~ "Basis",
    scenario == 2 ~ "Defer 2nd dose",
    scenario == 3 ~ "Basis (Vasileiou VE)",
    scenario == 4 ~ "Defer 2nd dose (Vasileiou VE)"
  ),
  analysis = "dose deferral")

results_all <- bind_rows(results_main, results_sa20, results_sa50)
# summary table -------------------------------------------------------------------
results_dat <- results_defer

summary_tab <- results_dat %>%
  filter(date >= as.Date("2021-04-01"),
         date < as.Date("2021-09-01"),
         outcome %in% c("new_cases","hospital_admissions")) %>%
  group_by(analysis, scenario, outcome) %>%
  summarise_at(.vars = "value", .funs = sum)

ref_no_vac <- summary_tab %>%
  filter(scenario == "no vaccination (healthy adults)") %>%
  group_by(analysis, outcome) %>%
  summarise_at(.vars = "value", .funs = sum)

infections_sum <- summary_tab %>%
  filter(outcome == "new_infections") %>%
  group_by(analysis) %>%
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
for_plot <- results_dat %>%
  filter(outcome %in% c(#"new_infections", 
                        #"new_cases", 
                        "hospital_admissions"#, 
                        #"ic_admissions", "new_deaths"
                        )) %>%
  group_by(outcome, scenario) %>%
  mutate(rolling_avg = zoo::rollmean(value, k = 7, fill = NA),
         outcome = case_when(
           outcome == "hospital_admissions" ~ "Hospital Admissions",
           # outcome == "ic_admissions" ~ "IC Admissions",
           outcome == "new_cases" ~ "New Cases",
           # outcome == "new_infections" ~ "New Infections",
           # outcome == "new_deaths" ~ "New Deaths"
         )#,
         # outcome = factor(outcome, levels = c("New Infections", "New Cases", "Hospital Admissions", 
         #                                      "IC Admissions", "New Deaths")),
         # scenario = case_when(
         #   scenario == "alternative" ~ "Alternative",
         #   scenario == "no vaccination (healthy adults)" ~ "No Vaccination (Healthy Adults)",
         #   scenario == "old to young" ~ "Old to Young",
         #   scenario == "young to old" ~ "Young to Old",
         #   scenario == "no vaccination (at all)" ~ "No Vaccination (At All)"
         # ),
         # scenario = factor(scenario, levels = c("Old to Young", "Young to Old", "Alternative", 
         #                                        "No Vaccination (Healthy Adults)", "No Vaccination (At All)"))
         )
# plot
for_plot1 <- for_plot %>% 
  filter(date >= as.Date("2021-04-01"),
         date < as.Date("2021-09-01"),
         scenario %in% c("Basis", "Defer 2nd dose"))
custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
                "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
  #filter(scenario != "No Vaccination (At All)")
g_sum <- ggplot(for_plot1, 
                aes(x = date, y = rolling_avg, color = scenario)) +
  geom_line(#aes(linetype = "Main"),
            position=position_dodge(width=1)
            ) +
  # geom_line(data = for_plot1 %>% filter(analysis == "sa20"), 
  #           aes(x = date, y = value, color = scenario, linetype = "Reduction of 20%"), 
  #           #linetype = "dashed", 
  #           position=position_dodge(width=2.5)) +
  # geom_line(data = for_plot1 %>% filter(analysis == "sa50"), 
  #           aes(x = date, y = value, color = scenario, linetype = "Reduction of 50%"), 
  #           #linetype = "dotted", 
  #           position=position_dodge(width=2.5)) +
  labs(y = "Value", x = "Time (days)", color = "Vaccination Scenario") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  #scale_color_manual(values = custom.col[c(5,6)]) + #6, 8 
  #geom_vline(xintercept = as.Date("2021-04-25"),linetype="dashed", color = "grey70") +
  #geom_hline(yintercept = 6214.508, linetype = "dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~outcome, scales = "free")
g_sum

# save plot to file
ggsave("inst/extdata/results/defer_2nd_dose_02April.jpg",
       plot = g_sum,
       height = 6,
       width = 12,
       dpi = 300)

# calculate life years lost ---------------------------------------------------------
# data wrangling --------------------------------------------------------------------
results_by_ag <- bind_rows(o2y_sa$df,y2o_sa$df,alt_sa$df,nov_healthy_sa$df, nov_sa$df,.id = "scenario") %>%
  mutate(scenario = case_when(
    scenario == 1 ~ "old to young",
    scenario == 2 ~ "young to old",
    scenario == 3 ~ "alternative",
    scenario == 4 ~ "no vaccination (healthy adults)",
    scenario == 5 ~ "no vaccination (at all)"
  ))

# calculation -----------------------------------------------------------------------
life_expectancy <- c(77.89, 67.93, 58.08, 48.28, 38.6, 29.22, 20.52, 12.76, 4.35)
deaths_by_ag <- results_by_ag %>%
  group_by(scenario, outcome, age_group) %>%
  summarise_at(.vars = "value", .funs = sum) %>%
  filter(outcome == "new_deaths")

no_vac_deaths_by_ag <- deaths_by_ag %>%
  filter(scenario == "no vaccination (at all)")

deaths_prevented <- deaths_by_ag %>%
  mutate(deaths_prevented = no_vac_deaths_by_ag$value - value,
         life_years_lost_prevented = deaths_prevented * life_expectancy,
         life_years_lost = value * life_expectancy)

total_life_years_lost <- deaths_prevented %>%
  group_by(scenario) %>%
  summarise_at(.vars = c("value","life_years_lost", "life_years_lost_prevented"), .funs = sum) %>%
  mutate(perc_diff = ((life_years_lost*100)/life_years_lost[2]) - 100)

