# Some wrangling of contact matrices
library(dplyr)
library(tidyr)
# read in all contact matrices
contact_matrices_all <- read.delim("inst/extdata/data/S2_contact_matrices_withPico3_10y.tsv", header=TRUE, allowEscapes=FALSE, sep="\t")

# baseline (2017)
baseline <- contact_matrices_all %>%
  filter(survey == "baseline") %>%
  filter(contact_type == "all") %>%
  select(-survey, -contact_type) %>%
  pivot_wider(., names_from = cont_age, values_from = m_est)
  
# April 2020 (lockdown)
april2020 <- contact_matrices_all %>%
  filter(survey == "April 2020") %>%
  filter(contact_type == "all") %>%
  select(-survey, -contact_type) %>%
  pivot_wider(., names_from = cont_age, values_from = m_est)

new_row <- april2020[1,-1] + april2020[2,-1]
april2020 <- april2020 %>%
  add_row(part_age = "[0,10)", `[0,10)` = new_row$`[0,10)`, `[10,20)` = new_row$`[10,20)`,
          `[20,30)` = new_row$`[20,30)`, `[30,40)` = new_row$`[30,40)`, 
          `[40,50)` = new_row$`[40,50)`, `[50,60)` = new_row$`[50,60)`, 
          `[60,70)` = new_row$`[60,70)`, `[70,80)` = new_row$`[70,80)`, 
          `[80,Inf]` = new_row$`[80,Inf]`, .before = 1) %>%
  slice(-(2:3))

# June 2020 
june2020 <- contact_matrices_all %>%
  filter(survey == "June 2020") %>%
  filter(contact_type == "all") %>%
  select(-survey, -contact_type) %>%
  pivot_wider(., names_from = cont_age, values_from = m_est)

# September 2020 
september2020 <- contact_matrices_all %>%
  filter(survey == "September 2020") %>%
  filter(contact_type == "all") %>%
  select(-survey, -contact_type) %>%
  pivot_wider(., names_from = cont_age, values_from = m_est) 

# put in a list and write to rds
cm <- list(baseline = baseline,
           april2020 = april2020,
           june2020 = june2020,
           september2020 = september2020)
saveRDS(cm,"inst/extdata/data/contact_matrices_for_model_input.rds")
