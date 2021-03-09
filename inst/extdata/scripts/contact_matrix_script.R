# Some wrangling of contact matrices
library(dplyr)
library(tidyr)
# read in all contact matrices
contact_matrices_all <- read.delim("inst/extdata/data/S2_contact_matrices_withPico3_10y.tsv", header=TRUE, allowEscapes=FALSE, sep="\t")

# demographic info to make contact matrices symmetrical
# population distributions
demo <- data.frame(age_group = c("[0,10)", "[10,20)", "[20,30)", "[30,40)","[40,50)", "[50,60)",
                                 "[60,70)", "[70,80)", "[80,Inf]"), 
                   frac_2017 = c(0.105, 0.118, 0.126, 0.120, 0.138, 0.145, 0.122, 0.0808, 0.0447), 
                   frac_2019 = c(0.103, 0.116, 0.127, 0.122, 0.131, 0.145, 0.121, 0.0881, 0.0462),
                   frac_2020 = c(0.102, 0.115, 0.128, 0.123, 0.127, 0.145, 0.121, 0.0904, 0.0472))

# baseline (2017)
baseline <- contact_matrices_all %>%
  filter(survey == "baseline") %>%
  filter(contact_type == "all") %>%
  mutate(f_pop = rep(demo$frac_2017,9),
         c_est = round(m_est/f_pop,1)) %>%
  select(-survey, -contact_type, -m_est, -f_pop) %>%
  pivot_wider(., names_from = cont_age, values_from = c_est)
  
# April 2020 (lockdown)
april2020 <- contact_matrices_all %>%
  filter(survey == "April 2020") %>%
  filter(contact_type == "all") %>%
  mutate(f_pop = rep(demo$frac_2019,9),
         c_est = round(m_est/f_pop,1)) %>%
  select(-survey, -contact_type, -m_est, -f_pop) %>%
  pivot_wider(., names_from = cont_age, values_from = c_est)

# June 2020 
june2020 <- contact_matrices_all %>%
  filter(survey == "June 2020") %>%
  filter(contact_type == "all") %>%
  mutate(f_pop = rep(demo$frac_2019,9),
         c_est = round(m_est/f_pop,1)) %>%
  select(-survey, -contact_type, -m_est, -f_pop) %>%
  pivot_wider(., names_from = cont_age, values_from = c_est)

# September 2020 
september2020 <- contact_matrices_all %>%
  filter(survey == "September 2020") %>%
  filter(contact_type == "all") %>%
  mutate(f_pop = rep(demo$frac_2019,9),
         c_est = round(m_est/f_pop,1)) %>%
  select(-survey, -contact_type, -m_est, -f_pop) %>%
  pivot_wider(., names_from = cont_age, values_from = c_est) 

# February 2021
february2021 <- readRDS("inst/extdata/data/Contactpatterns_PICO4_10y.rds") %>%
  select(part_age, cnt_age, c_smt) %>%
  mutate(contact_type = c(rep("all", 81),
                          rep("community", 81),
                          rep("household", 81))) %>%
  filter(contact_type == "all") %>%
  select(-contact_type) %>%
  pivot_wider(., names_from = cnt_age, values_from = c_smt) 

# put in a list and write to rds
cm <- list(baseline = baseline,
           april2020 = april2020,
           june2020 = june2020,
           september2020 = september2020,
           february2021 = february2021)
saveRDS(cm,"inst/extdata/data/contact_matrices_for_model_input.rds")
