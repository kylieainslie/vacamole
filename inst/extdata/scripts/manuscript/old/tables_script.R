# -----------------------------------------------------------
# Manuscript tables
# -----------------------------------------------------------

source("/inst/extdata/scripts/manuscript/data_wrangling_for_figures.R")
#save_path <- "inst/extdata/results"

# percent reduction ------------------------------------------------------------
# total population
peak_inc <- df_all %>%
  group_by(scenario_id, target_variable, date, sample) %>%
  summarise(sum = sum(value)) %>%
  ungroup() %>%
  group_by(scenario_id, target_variable, sample) %>%
  summarise(max = max(sum)) %>%
  ungroup() %>%
  group_by(scenario_id, target_variable) %>%
  summarise(max_mean = mean(max),
            max_q025 = quantile(max, probs = 0.025),
            max_q975 = quantile(max, probs = 0.975))

peak_inc_12v18 <- peak_inc %>%
  filter(scenario_id %in% c("Vaccination in 18+", "Vaccination in 12+")) %>%
  group_by(target_variable) %>%
  mutate(perc_change_mean = scales::percent((max_mean - max_mean[scenario_id == "Vaccination in 18+"])/max_mean[scenario_id == "Vaccination in 18+"]),
         perc_change_q025 = scales::percent((max_q025 - max_q025[scenario_id == "Vaccination in 18+"])/max_q025[scenario_id == "Vaccination in 18+"]),
         perc_change_q975 = scales::percent((max_q975 - max_q975[scenario_id == "Vaccination in 18+"])/max_q975[scenario_id == "Vaccination in 18+"]))

peak_inc_5v18 <- peak_inc %>%
  filter(scenario_id %in% c("Vaccination in 18+", "Vaccination in 5+")) %>%
  group_by(target_variable) %>%
  mutate(perc_change_mean = scales::percent((max_mean - max_mean[scenario_id == "Vaccination in 18+"])/max_mean[scenario_id == "Vaccination in 18+"]),
         perc_change_q025 = scales::percent((max_q025 - max_q025[scenario_id == "Vaccination in 18+"])/max_q025[scenario_id == "Vaccination in 18+"]),
         perc_change_q975 = scales::percent((max_q975 - max_q975[scenario_id == "Vaccination in 18+"])/max_q975[scenario_id == "Vaccination in 18+"]))


# Table 1 ---------------------------------------------------
table1 <- all_res_for_plot %>%
  filter(outcome != "Daily Deaths",
         Immunity = "No Waning") %>%
  group_by(Scenario, age_group2, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

# calculate percent difference
table1a <- table1 %>%
  group_by(age_group2, outcome) %>%
  mutate(abs_diff = mle - mle[Scenario == "Vaccination of 18+"],
         abs_diff_lower = lower - lower[Scenario == "Vaccination of 18+"],
         abs_diff_upper = upper - upper[Scenario == "Vaccination of 18+"],
         perc_diff = (mle * 100)/mle[Scenario == "Vaccination of 18+"] - 100,
         perc_diff_lower = (lower * 100)/lower[Scenario == "Vaccination of 18+"] - 100,
         perc_diff_upper = (upper * 100)/upper[Scenario == "Vaccination of 18+"] - 100) %>%
  mutate_if(is.numeric, round, 1) %>%
  as.data.frame()

# write.csv(table1, file = paste0(save_path, "table1.csv"))

# overall 
table1 %>%
  group_by(Scenario, outcome) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  group_by(outcome) %>%
  mutate(abs_diff = mle - mle[Scenario == "18+"],
         abs_diff_lower = lower - lower[Scenario == "18+"],
         abs_diff_upper = upper - upper[Scenario == "18+"],
         perc_diff = (mle * 100)/mle[Scenario == "18+"] - 100,
         perc_diff_lower = (lower * 100)/lower[Scenario == "18+"] - 100,
         perc_diff_upper = (upper * 100)/upper[Scenario == "18+"] - 100) %>%
  mutate_if(is.numeric, round, 1) %>%
  as.data.frame()

# Table 2 ---------------------------------------------------
table2 <- all_res_for_plot %>%
  filter(outcome != "Daily Deaths",
         Immunity = "Waning") %>%
  group_by(Scenario, age_group2, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum") %>%
  ungroup()

# calculate percent difference
table2a <- table2 %>%
  group_by(age_group2, outcome) %>%
  mutate(abs_diff = mle - mle[Scenario == "18+"],
         abs_diff_lower = lower - lower[Scenario == "18+"],
         abs_diff_upper = upper - upper[Scenario == "18+"],
         perc_diff = (mle * 100)/mle[Scenario == "18+"] - 100,
         perc_diff_lower = (lower * 100)/lower[Scenario == "18+"] - 100,
         perc_diff_upper = (upper * 100)/upper[Scenario == "18+"] - 100) %>%
  mutate_if(is.numeric, round, 1) %>%
  as.data.frame()

#write.csv(table2a, file = paste0(file_path, "table2.csv"))

# all age groups together
table2 %>%
  group_by(Scenario, outcome) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  group_by(outcome) %>%
  mutate(abs_diff = mle - mle[Scenario == "18+"],
         abs_diff_lower = lower - lower[Scenario == "18+"],
         abs_diff_upper = upper - upper[Scenario == "18+"],
         perc_diff = (mle * 100)/mle[Scenario == "18+"] - 100,
         perc_diff_lower = (lower * 100)/lower[Scenario == "18+"] - 100,
         perc_diff_upper = (upper * 100)/upper[Scenario == "18+"] - 100) %>%
  mutate_if(is.numeric, round, 1) %>%
  as.data.frame()