age_dist_10 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
                 0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
n <- 17407585 # Dutch population size
n_vec_10 <- n * age_dist_10

basis_19April <- basis %>% 
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(date == as.Date("2021-04-19")) %>%
  select(date:ja_d2_10)

basis_30Sept <- basis %>% 
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(date == as.Date("2021-09-30")) %>%
  select(date:ja_d2_10)

az_30Sept <- az %>% 
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(date == as.Date("2021-09-30")) %>%
  select(date:ja_d2_10)

mRNA_30Sept <- mRNA %>% 
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(date == as.Date("2021-09-30")) %>%
  select(date:ja_d2_10)

janssen60_30Sept <- janssen60 %>% 
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(date == as.Date("2021-09-30")) %>%
  select(date:ja_d2_10)

janssen50_30Sept <- janssen50 %>% 
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(date == as.Date("2021-09-30")) %>%
  select(date:ja_d2_10)

janssen40_30Sept <- janssen40 %>% 
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  filter(date == as.Date("2021-09-30")) %>%
  select(date:ja_d2_10)


total_vac_num <- bind_rows(sweep(basis_30Sept[,-1], 2, n_vec_10, "*"),
                           sweep(az_30Sept[,-1], 2, n_vec_10, "*"),
                           sweep(mRNA_30Sept[,-1], 2, n_vec_10, "*"),
                           sweep(janssen60_30Sept[,-1], 2, n_vec_10, "*"),
                           sweep(janssen50_30Sept[,-1], 2, n_vec_10, "*"),
                           sweep(janssen40_30Sept[,-1], 2, n_vec_10, "*"),
                           .id = "scenario") %>%
  mutate( date = c(rep(basis_30Sept$date,6)),
    scenario = case_when(
    scenario == 1 ~ "basis",
    scenario == 2 ~ "az" ,
    scenario == 3 ~ "mRNA",
    scenario == 4 ~ "janssen60",
    scenario == 5 ~ "janssen50",
    scenario == 6 ~ "janssen40"
  )) %>%
  pivot_longer(cols = pf_d1_1:ja_d2_10,names_to = c("vaccine", "dose", "age_group"),
               names_pattern = "(.*)_(.*)_(.*)", values_to = "n") %>%
  mutate(vaccine = case_when(
    vaccine == "pf" ~ "Pfizer",
    vaccine == "mo" ~ "Moderna",
    vaccine == "az" ~ "AstraZeneca",
    vaccine == "ja" ~ "Janssen"),
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
    ),
    n_red = n/1000
  )

# calculate total doses
total_doses <- total_vac_num %>%
  #filter(vaccine == "pfizer") %>%
  group_by(date, scenario, dose, vaccine) %>%
  summarise_at(.vars = "n", sum) %>%
  mutate(n_red = n/1000)

p_total <- ggplot(data = total_doses, aes(x=dose, y=n_red, fill=vaccine)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=round(n_red,0)), position=position_dodge(width=0.9), vjust=-0.25) +
  ylab("Number of doses (x 1000)") +
  theme(legend.position = "bottom",
        panel.background = element_blank()) +
  facet_wrap(~scenario)
p_total

# save plot to file
ggsave("inst/extdata/results/total_vac_allocation_plot_07May.jpg",
       plot = p_total,
       height = 8,
       width = 12,
       dpi = 300)


p_total_age_dose1 <- ggplot(data = total_vac_num %>%
                              filter(dose == 1), aes(x=age_group, y=n_red, fill=vaccine)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ylab("Number of first doses (x 1000)") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~scenario)
p_total_age_dose1

# save plot to file
ggsave("inst/extdata/results/total_vac_allocation_plot_by_age_dose1_07May.jpg",
       plot = p_total_age_dose1,
       height = 8,
       width = 12,
       dpi = 300)


p_total_age_dose2 <- ggplot(data = total_vac_num %>%
                              filter(dose == 2), aes(x=age_group, y=n_red, fill=vaccine)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ylab("Number of second doses (x 1000)") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~scenario)
p_total_age_dose2

# save plot to file
ggsave("inst/extdata/results/total_vac_allocation_plot_by_age_dose2_07May.jpg",
       plot = p_total_age_dose2,
       height = 8,
       width = 12,
       dpi = 300)

# --------------------------------------------
# Calculate adverse events for AZ and Janssen
# --------------------------------------------
# vectors of tts_rates by 10 year age group
# rates for AZ are per 100,000
tts_rates_az <- data.frame(age_group = c("0-9","10-19","20-29","30-39",
                                         "40-49","50-59","60-69","70-79",
                                         "80+"),
                           raw = c(0, 0, 1.9, 1.8, 2.1, 1.1, 1.0, 0.5, 0.4),
                           corrected = c(0, 0, 2.7, 2.1, 2.7, 1.4, 1.3, 0.9, 0.4),
                           fatalities = c(0, 0, 0.1, 0.5, 0.9, 0.2, 0.3, 0.2, 0.2)
                           )
az_dose1 <- total_vac_num %>% 
  filter(vaccine == "AstraZeneca",
         dose == 1, 
         age_group != "90+") %>% # no AZ is given to 90+, so won't matter if we exclude it
  select(scenario, age_group, n) %>%
  group_by(scenario) %>%
  mutate(ae_raw = (n/100000) * tts_rates_az$raw,
         ae_corrected = (n/100000) * tts_rates_az$corrected,
         ae_fatal = (n/100000) * tts_rates_az$fatalities
         ) 

tab <- az_dose1 %>%
  summarise_at(.vars = c("n", "ae_raw", "ae_corrected", "ae_fatal"),
               .funs = "sum")

# janssen tts events
# rates are per 1 million
ja_tts_rates <- c(0, 0.51, 2.55, 6.94, 2.53, 0.59, 0.29, 0, 0)
ja_dose1 <- total_vac_num %>% 
  filter(vaccine == "Janssen",
         dose == 1, 
         age_group != "90+") %>% # no Janssen is given to 90+, so won't matter if we exclude it
  select(scenario, age_group, n) %>%
  group_by(scenario) %>%
  mutate(ae = (n/1000000) * ja_tts_rates) 

tab2 <- ja_dose1 %>%
  summarise_at(.vars = c("n", "ae"),
               .funs = "sum")




