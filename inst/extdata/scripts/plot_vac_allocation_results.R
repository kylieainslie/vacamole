# Plot results of different AstraZeneca vaccine allocation strategies--------------
# 1) old to young
# 2) young to old
# 3) alternative
# 4) no vaccination (control)

# load results --------------------------------------------------------------------
# main analysis
# basis_res <- readRDS("inst/extdata/results/basis_no_child_vac_14june.rds") 
# july_res <- readRDS("inst/extdata/results/basis_w_waning_child_vac_july_start_cov07_9june.rds") 
# aug_res <- readRDS("inst/extdata/results/basis_w_waning_child_vac_aug_start_cov07_9june.rds") 
# sept_res <- readRDS("inst/extdata/results/basis_w_waning_child_vac_sept_start_cov07_9june.rds") 

# interplay between waning and child vax
basis_no_child_vac <- readRDS("inst/extdata/results/basis_no_child_vac_14june.rds")
basis_no_wane <- readRDS("inst/extdata/results/basis_no_wane_16june.rds")
basis_child_vac <- readRDS("inst/extdata/results/basis_child_vac_14june.rds")
basis_no_wane_child_vac <- readRDS("inst/extdata/results/basis_no_wane_child_vac_16june.rds")

# vac coverage sensitivity analysis
# july_cov07_res <- readRDS("inst/extdata/results/basis_w_waning_child_vac_july_start_cov07_9june.rds") 
# july_cov06_res <- readRDS("inst/extdata/results/basis_child_vac_july_start_cov06_9june.rds") 
# july_cov05_res <- readRDS("inst/extdata/results/basis_child_vac_july_start_cov06_9june.rds") 
# july_cov08_res <- readRDS("inst/extdata/results/basis_child_vac_july_start_cov08_9june.rds") 
# data wrangle for plotting/summary table -----------------------------------------
results_all <- bind_rows(#basis_no_child_vac$df_summary,
                         basis_no_wane$df_summary,
                         #basis_child_vac$df_summary,
                         basis_no_wane_child_vac$df_summary,
                          .id = "scenario") %>%
  mutate(scenario = case_when(
    #scenario == 1 ~ "No vaccination, waning",
    scenario == 1 ~ "No vaccination, no waning",
    #scenario == 3 ~ "Vaccination, waning",
    scenario == 2 ~ "Vaccination, no waning",
  ),
   scenario = factor(scenario, levels = c(#"No vaccination, waning",
                                          "No vaccination, no waning",
                                          #"Vaccination, waning",
                                          "Vaccination, no waning"
                                          ))
  )

results_sa <- bind_rows(basis_res$df_summary,
                        july_res$df_summary,
                        basis_res_no_wane$df_summary,
                        july_res_no_wane$df_summary,
                        .id = "scenario") %>%
  mutate(scenario = case_when(
    scenario == 1 ~ "No vaccination, waning",
    scenario == 2 ~ "Vacination, waning",
    scenario == 3 ~ "No vaccination, no waning",
    scenario == 4 ~ "Vaccination, no waning",
  ),
  scenario = factor(scenario, levels = c("No vaccination, waning",
                                         "Vaccination, waning",
                                         "No vaccination, no waning",
                                         "Vaccination, no waning")))

results_sa2 <- bind_rows(july_cov07_res$df_summary,
                         july_cov06_res$df_summary,
                         july_cov05_res$df_summary,
                         july_cov08_res$df_summary,
                        .id = "scenario") %>%
  mutate(scenario = case_when(
    scenario == 1 ~ "Coverage = 70%",
    scenario == 2 ~ "Coverage = 60%",
    scenario == 3 ~ "Coverage = 50%",
    scenario == 4 ~ "Coverage = 80%",
  ),
  scenario = factor(scenario, levels = c("Coverage = 50%",
                                         "Coverage = 60%",
                                         "Coverage = 70%",
                                         "Coverage = 80%")))

# by age group
results_ag <- bind_rows(#basis_no_child_vac$df,
                        basis_no_wane$df,
                        #basis_child_vac$df,
                        basis_no_wane_child_vac$df,
                        .id = "scenario") %>%
  mutate(scenario = case_when(
    #scenario == 1 ~ "No vaccination, waning",
    scenario == 1 ~ "No vaccination, no waning",
    #scenario == 3 ~ "Vaccination, waning",
    scenario == 2 ~ "Vaccination, no waning",
  ),
  scenario = factor(scenario, levels = c(#"No vaccination, waning",
                                         "No vaccination, no waning",
                                         #"Vaccination, waning",
                                         "Vaccination, no waning"
  )),
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
         outcome %in% c("new_cases","hospital_admissions", "ic_admissions","new_deaths")) %>%
  group_by(scenario, outcome) %>%
  summarise_at(.vars = "value", .funs = sum)

summary_tab_ag <- results_ag %>%
  filter(outcome %in% c("new_cases","hospital_admissions", "ic_admissions","new_deaths"),
         age_group != "10-19"
         ) %>%
  group_by(scenario,
           #date,
           #age_group, 
           outcome) %>%
  summarise_at(.vars = "value", .funs = sum) 

percent_diff <- function(x, ref){
  (x * 100)/ref - 100
}
# plot ----------------------------------------------------------------------------
# subset for plotting 
for_plot <- summary_tab_ag %>% #results_ag
  filter(outcome %in% c("new_cases", 
                        "hospital_admissions", 
                        "ic_admissions", 
                        "new_deaths"
                        )#,
         # date > as.Date("2021-10-28")#,
         # scenario %in% c("No vaccination, waning",
         #                 "Vaccination, waning"
         #                 ),
         #age_group == "10-19"
         ) %>%
  group_by(outcome, scenario) %>%
  mutate(#rolling_avg = zoo::rollmean(value, k = 7, fill = NA),
         outcome = case_when(
           outcome == "hospital_admissions" ~ "Hospital Admissions",
           outcome == "ic_admissions" ~ "IC Admissions",
           outcome == "new_cases" ~ "New Cases",
           # outcome == "new_infections" ~ "New Infections",
           outcome == "new_deaths" ~ "New Deaths"
         ),
         outcome = factor(outcome, levels = c("New Cases", 
                                              "Hospital Admissions", 
                                              "IC Admissions", 
                                              "New Deaths"
                                              )),
         scenario = case_when(
           scenario == "No vaccination, no waning" ~ "No vaccination of 12-17 year olds",
           scenario == "Vaccination, no waning" ~ "Vaccination of 12-17 year olds (mid-July start)"
         )
         )

# overall
library(scales)
custom.col <- hue_pal()(4)

g_sum <- ggplot(for_plot, 
                aes(x = date, y = value, color = scenario)) +
  geom_line(#aes(linetype = "Main"),
    #position=position_dodge(width=2)
  ) +
  labs(y = "Value", x = "Time (days)", color = "Vaccination Scenario") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  scale_color_manual(values = custom.col[c(1,2)]) + #6, 8 
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
#g_sum <- g_sum + guides(color=guide_legend(nrow=2,byrow=TRUE))
g_sum
# save plot to file
ggsave("inst/extdata/results/indirect_effects_no_waning_beta1_0_15_plot_16june.jpg",
       plot = g_sum,
       height = 12,
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
  select(date, cases, scenario)
  
july_start <- model_forward_cases %>%
  filter(scenario == "Vaccination of 12-17 year olds (mid-July start)") %>%
  select(-scenario) %>%
  bind_rows(model_fit, .) %>%
  mutate(date = time + as.Date("2021-01-31"),
         scenario = "Vaccination of 12-17 year olds (mid-July start)") %>%
  select(date, cases, scenario)

for_plotting <- bind_rows(no_vac, july_start)

p_fit <- ggplot(for_plotting, aes(x = date, y = cases, color = scenario)) +
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
tag <- "9june"
ggsave(paste0("inst/extdata/results/child_vac_figure2_",tag,".jpg"),
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

