july1_orig <- cum_vac_schedule_orig %>% 
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(date == as.Date("2021-07-01"))

july1_delay3mo <- cum_vac_schedule_delay3mo %>% 
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(date == as.Date("2021-07-01"))

july1_no2dose <- cum_vac_schedule_no2dose %>% 
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(date == as.Date("2021-07-01"))

total_vac_num <- bind_rows(sweep(july1_orig[,-1], 2, n_vec_10, "*"),
                           sweep(july1_delay3mo[,-1], 2, n_vec_10, "*"),
                           sweep(july1_no2dose[,-1], 2, n_vec_10, "*"),
                           .id = "scenario") %>%
  mutate(scenario = case_when(
    scenario == 1 ~ "original",
    scenario == 2 ~ "delay 3 months",
    scenario == 3 ~ "no 2nd doses"
  )) %>%
  pivot_longer(cols = pf_d1_1:az_d2_10,names_to = c("vaccine", "dose", "age_group"),
               names_pattern = "(.*)_(.*)_(.*)", values_to = "n") %>%
  mutate(vaccine = case_when(
    vaccine == "pf" ~ "pfizer",
    vaccine == "mo" ~ "moderna",
    vaccine == "az" ~ "astrazeneca"),
    dose = case_when(
      dose == "d1" ~ 1,
      dose == "d2" ~ 2),
    dose = as.factor(dose),
    age_group = case_when(
      age_group == 1 ~ "0-9",
      age_group == 2 ~ "10-19",
      age_group == 3 ~ "20-29",
      age_group == 4 ~ "30-39",
      age_group == 5 ~ "40-49",
      age_group == 6 ~ "50-59",
      age_group == 7 ~ "60-69",
      age_group == 8 ~ "70-79",
      age_group == 9 ~ "80-89",
      age_group == 10 ~ "90+"
    )
  )

# calculate total pfizer doses
total_vac_num %>%
  filter(vaccine == "pfizer") %>%
  group_by(scenario, dose) %>%
  summarise_at(.vars = "n", sum)

# look at doses by week
delay3mo <- cum_vac_schedule_delay3mo %>% 
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(date <= as.Date("2021-07-01"))
vac_num_delay3mo <- sweep(july1_delay3mo[,-1], 2, n_vec_10, "*") %>%
  mutate(date = as.Date(cum_vac_schedule_delay3mo$date[1:179], "%m/%d/%Y"))
