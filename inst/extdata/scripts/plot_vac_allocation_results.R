# Plot results of different AstraZeneca vaccine allocation strategies--------------
# 1) old to young
# 2) young to old
# 3) alternative
# 4) no vaccination (control)

# load results --------------------------------------------------------------------
# main analysis
basis_res <- readRDS("inst/extdata/results/basis_w_waning_no_child_vac_9june.rds") 
july_res <- readRDS("inst/extdata/results/basis_w_waning_child_vac_july_start_cov07_9june.rds") 
aug_res <- readRDS("inst/extdata/results/basis_w_waning_child_vac_aug_start_cov07_9june.rds") 
sept_res <- readRDS("inst/extdata/results/basis_w_waning_child_vac_sept_start_cov07_9june.rds") 

# data wrangle for plotting/summary table -----------------------------------------
results_all <- bind_rows(basis_res$df_summary,
                         july_res$df_summary,
                         aug_res$df_summary,
                         sept_res$df_summary,
                          .id = "scenario") %>%
  mutate(scenario = case_when(
    scenario == 1 ~ "No vaccination of 12-17 year olds",
    scenario == 2 ~ "Vaccination of 12-17 year olds (mid-July start)",
    scenario == 3 ~ "Vaccination of 12-17 year olds (August start)",
    scenario == 4 ~ "Vaccination of 12-17 year olds (September start)",
  ),
   scenario = factor(scenario, levels = c("No vaccination of 12-17 year olds",
                                          "Vaccination of 12-17 year olds (mid-July start)",
                                          "Vaccination of 12-17 year olds (August start)",
                                          "Vaccination of 12-17 year olds (September start)")))

# by age group
results_ag <- bind_rows(basis_res$df,
                         az_res$df,
                         mRNA_res$df,
                         janssen60_res$df,
                         janssen50_res$df,
                         janssen40_res$df,
                         .id = "scenario") %>%
  mutate(scenario = case_when(
    scenario == 1 ~ "Basis",
    scenario == 2 ~ "AZ",
    scenario == 3 ~ "mRNA",
    scenario == 4 ~ "Janssen60",
    scenario == 5 ~ "Janssen50",
    scenario == 6 ~ "Janssen40"
  ),
    scenario = factor(scenario, levels = c("Basis", "AZ", "mRNA",
                                           "Janssen60", "Janssen50", "Janssen40")),
    age_group = case_when(
      age_group == 1 ~ "0-9",
      age_group == 2 ~ "10-19",
      age_group == 3 ~ "20-29",
      age_group == 4 ~ "30-39",
      age_group == 5 ~ "40-49",
      age_group == 6 ~ "50-59",
      age_group == 7 ~ "60-69",
      age_group == 8 ~ "70-79",
      age_group == 9 ~ "80+",
    )) 

# summary table -------------------------------------------------------------------
summary_tab <- results_all %>%
  filter(# date >= as.Date("2021-04-01"),
         # date < as.Date("2021-09-01"),
         outcome %in% c("new_cases","hospital_admissions", "new_deaths")) %>%
  group_by(scenario, outcome) %>%
  summarise_at(.vars = "value", .funs = sum)

summary_tab_ag <- results_ag %>%
  filter(outcome %in% c("new_cases","hospital_admissions", "new_deaths")) %>%
  group_by(scenario, 
           #age_group, 
           outcome) %>%
  summarise_at(.vars = "value", .funs = sum) 

# plot ----------------------------------------------------------------------------
# subset for plotting 
for_plot <- results_all %>% #results_ag
  filter(outcome %in% c("new_cases"#, 
                        #"hospital_admissions", 
                        #"ic_admissions", 
                        #"new_deaths"
                        ),
         scenario %in% c("No vaccination of 12-17 year olds",
                         "Vaccination of 12-17 year olds (mid-July start)")) %>%
  group_by(outcome, scenario) %>%
  mutate(rolling_avg = zoo::rollmean(value, k = 7, fill = NA),
         outcome = case_when(
           #outcome == "hospital_admissions" ~ "Hospital Admissions",
           # outcome == "ic_admissions" ~ "IC Admissions",
           outcome == "new_cases" ~ "New Cases"#,
           #outcome == "new_infections" ~ "New Infections",
           #outcome == "new_deaths" ~ "New Deaths"
         ),
         outcome = factor(outcome, levels = c("New Cases", "Hospital Admissions", 
                                               #"IC Admissions", 
                                              "New Deaths"))
         )

# overall
g_sum <- ggplot(for_plot, 
                aes(x = date, y = value, color = scenario)) +
  geom_line(#aes(linetype = "Main"),
    position=position_dodge(width=2)
  ) +
  labs(y = "Value", x = "Time (days)", color = "Vaccination Scenario") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  #scale_color_manual(values = custom.col[c(5,6)]) + #6, 8 
  #geom_vline(xintercept = as.Date("2021-04-25"),linetype="dashed", color = "grey70") +
  #geom_hline(yintercept = 6214.508, linetype = "dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
  ) +
  facet_wrap(~ outcome, 
             scales = "free", 
             ncol = 1
  )
g_sum

# save plot to file
ggsave("inst/extdata/results/res_plot_figure1_9june.jpg",
       plot = g_sum,
       height = 8,
       width = 12,
       dpi = 300)

# combine with model fit
model_fit$time <- model_fit$time - 1

model_forward_cases <- results_all %>%
  ungroup() %>%
  filter(outcome == "new_cases") %>%
  select(scenario, time, value) %>%
  filter(time != tail(model_fit$time,1)) %>% # remove first time point because it's a repeat of the last time point of the data fit
  rename(cases = value)

no_vac <- model_forward_cases %>%
  filter(scenario == "No vaccination of 12-17 year olds") %>%
  select(-scenario) %>%
  bind_rows(model_fit, .) %>%
  mutate(date = time + as.Date("2021-01-31"),
         scenario = "No vaccination of 12-17 year olds") %>%
  select(date, cases)
  
july_start <- model_forward_cases %>%
  filter(scenario == "Vaccination of 12-17 year olds (mid-July start)") %>%
  select(-scenario) %>%
  bind_rows(model_fit, .) %>%
  mutate(date = time + as.Date("2021-01-31"),
         scenario = "Vaccination of 12-17 year olds (mid-July start)") %>%
  select(date, cases)

p_fit <- ggplot(cases_fit_and_model, aes(x = date, y = cases)) +
  geom_line() +
  geom_point(data = osiris1 %>% filter(date < as.Date("2021-05-26")), aes(x = date, y = inc, color = "Osiris data")) +
  labs(y = "Daily Cases", x = "Time (days)") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  #geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_fit

# save plot to file
ggsave(paste0("inst/extdata/results/plot_",tag,".jpg"),
       plot = p_fit,
       height = 6,
       width = 12,
       dpi = 300)



# facet by age group
g_age <- ggplot(for_plot %>% filter(outcome == "Hospital Admissions"), 
                aes(x = date, y = value, color = scenario)) +
  geom_line(#aes(linetype = "Main"),
            position=position_dodge(width=2)
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
  scale_x_date(date_breaks = "month", date_labels = "%d %b") +
  #scale_color_manual(values = custom.col[c(5,6)]) + #6, 8 
  #geom_vline(xintercept = as.Date("2021-04-25"),linetype="dashed", color = "grey70") +
  #geom_hline(yintercept = 6214.508, linetype = "dashed", color = "grey70") +
  theme(legend.position = "bottom",
        #panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
        ) +
  facet_wrap(~ age_group, 
             #scales = "free", 
             ncol = 5
             )
g_age

# save plot to file
ggsave("inst/extdata/results/res_plot_hospital_admissions_10May.jpg",
       plot = g_age,
       height = 8,
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

