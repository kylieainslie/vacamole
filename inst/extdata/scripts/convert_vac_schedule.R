# convert cumulative vaccination schedule to non-cumulative ---------------------------------------
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
# read in vac schedule file -----------------------------------------------------------------------
#cum_vac_sched <- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210329 Basis.csv")

convert_vac_schedule <- function(vac_schedule){
# to combine age groups 9 and 10
age_dist_10 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
n <- 17407585 # Dutch population size
n_vec_10 <- n * age_dist_10

# take the difference for each row
# original
vac_schedule_orig <- data.frame(diff(as.matrix(vac_schedule[-1,-1]))) %>%
  add_row(vac_schedule[1,-1],.before = 1) %>%
  mutate(date = seq.Date(from = as.Date("2021-01-04"), to = as.Date("2021-12-30"), by = 1),
         pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         Ja_d1_9 = (Ja_d1_9 * n_vec_10[9] + Ja_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         Ja_d2_9 = (Ja_d2_9 * n_vec_10[9] + Ja_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
         ) %>%
  select(date, pf_d1_1:Ja_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10, 
         -az_d2_10, -Ja_d1_10, -Ja_d2_10) 

# filter for dates before 1 February 2021
before_feb <- vac_schedule_orig %>%
  filter(date < as.Date("2021-02-01")) %>%
  select(-date) %>%
  summarise_all(sum) %>%
  mutate(date = as.Date("2021-01-31")) %>%
  select(date, pf_d1_1:Ja_d2_9)

# filter so start date is 1 Feb 2021 with row for Jan 31 to reflect people who have already been vaccinated
vac_schedule_orig_new <- vac_schedule_orig %>%
  filter(date > as.Date("2021-01-31")) %>%
  add_row(before_feb, .before = 1)

vac_schedule_new_cs <- cumsum(vac_schedule_orig_new[,-1])

# pfizer
pf_dose1 <- vac_schedule_orig_new %>% select(date, pf_d1_1:pf_d1_9)
pf_dose1_cs <- vac_schedule_new_cs %>% select(pf_d1_1:pf_d1_9)
pf_dose2 <- vac_schedule_orig_new %>% select(date, pf_d2_1:pf_d2_9)
pf_dose2_cs <- vac_schedule_new_cs %>% select(pf_d2_1:pf_d2_9)

# moderna
mo_dose1 <- vac_schedule_orig_new %>% select(date, mo_d1_1:mo_d1_9)
mo_dose1_cs <- vac_schedule_new_cs %>% select(mo_d1_1:mo_d1_9)
mo_dose2 <- vac_schedule_orig_new %>% select(date, mo_d2_1:mo_d2_9)
mo_dose2_cs <- vac_schedule_new_cs %>% select(mo_d2_1:mo_d2_9)

# astrazeneca
az_dose1 <- vac_schedule_orig_new %>% select(date, az_d1_1:az_d1_9)
az_dose1_cs <- vac_schedule_new_cs %>% select(az_d1_1:az_d1_9)
az_dose2 <- vac_schedule_orig_new %>% select(date, az_d2_1:az_d2_9)
az_dose2_cs <- vac_schedule_new_cs %>% select(az_d2_1:az_d2_9)

# jansen
ja_dose1 <- vac_schedule_orig_new %>% select(date, Ja_d1_1:Ja_d1_9)
ja_dose1_cs <- vac_schedule_new_cs %>% select(Ja_d1_1:Ja_d1_9)
ja_dose2 <- vac_schedule_orig_new %>% select(date, Ja_d2_1:Ja_d2_9)
ja_dose2_cs <- vac_schedule_new_cs %>% select(Ja_d2_1:Ja_d2_9)


# calculate composite VE
name_suffix_d1 <- c(substr(names(pf_dose1[,-1]), 3, 7))
name_suffix_d2 <- c(substr(names(pf_dose2[,-1]), 3, 7))

alpha_dose1 <- pf_dose1[,-1] + mo_dose1[,-1] + az_dose1[,-1] + ja_dose1[,-1]
names(alpha_dose1) <- paste0("alpha", name_suffix_d1)
alpha_dose2 <- pf_dose2[,-1] + mo_dose2[,-1] + az_dose2[,-1] + ja_dose2[,-1]
names(alpha_dose2) <- paste0("alpha", name_suffix_d2)

total_dose1 <- pf_dose1_cs + mo_dose1_cs + az_dose1_cs + ja_dose1_cs
names(total_dose1) <- paste0("tot", name_suffix_d1)
total_dose2 <- pf_dose2_cs + mo_dose2_cs + az_dose2_cs + ja_dose2_cs
names(total_dose2) <- paste0("tot", name_suffix_d2)

frac_pf_dose1 <- pf_dose1_cs/total_dose1
frac_pf_dose1 <- ifelse(is.nan(frac_pf_dose1), 0, frac_pf_dose1)
frac_pf_dose2 <- pf_dose2_cs/total_dose2
#frac_pf_dose2 <- ifelse(is.nan(frac_pf_dose2), 0, frac_pf_dose2)

frac_mo_dose1 <- mo_dose1_cs/total_dose1
#frac_mo_dose1 <- ifelse(is.nan(frac_mo_dose1), 0, frac_mo_dose1)
frac_mo_dose2 <- mo_dose2_cs/total_dose2
#frac_mo_dose2 <- ifelse(is.nan(frac_mo_dose2), 0, frac_mo_dose2)

frac_az_dose1 <- az_dose1_cs/total_dose1
#frac_az_dose1 <- ifelse(is.nan(frac_az_dose1), 0, frac_az_dose1)
frac_az_dose2 <- az_dose2_cs/total_dose2
#frac_az_dose2 <- ifelse(is.nan(frac_az_dose2), 0, frac_az_dose2)

frac_ja_dose1 <- ja_dose1_cs/total_dose1
#frac_ja_dose1 <- ifelse(is.nan(frac_ja_dose1), 0, frac_az_dose1)
frac_ja_dose2 <- ja_dose2_cs/total_dose2
#frac_ja_dose2 <- ifelse(is.nan(frac_ja_dose2), 0, frac_ja_dose2)

comp_ve_dose1 <- frac_pf_dose1 * ve$pfizer[1] + 
  frac_mo_dose1 * ve$moderna[1] + 
  frac_az_dose1 * ve$astrazeneca[1] +
  frac_ja_dose1 * ve$jansen
comp_ve_dose2 <- frac_pf_dose2 * ve$pfizer[2] + 
  frac_mo_dose2 * ve$moderna[2] + 
  frac_az_dose2 * ve$astrazeneca[2] +
  frac_ja_dose2 * ve$jansen

# rate of hospitalisations multiplier
hosp_mult_dose1 <- frac_pf_dose1 * hosp_multiplier$pfizer[1] +
  frac_mo_dose1 * hosp_multiplier$moderna[1] +
  frac_az_dose1 * hosp_multiplier$astrazeneca[1] +
  frac_ja_dose1 * hosp_multiplier$jansen
hosp_mult_dose2 <- frac_pf_dose2 * hosp_multiplier$pfizer[2] +
  frac_mo_dose2 * hosp_multiplier$moderna[2] +
  frac_az_dose2 * hosp_multiplier$astrazeneca[2] +
  frac_ja_dose2 * hosp_multiplier$jansen

# composite delay to protection
delay_dose1 <- frac_pf_dose1 * delay$pfizer[1] +
  frac_mo_dose1 * delay$moderna[1] +
  frac_az_dose1 * delay$astrazeneca[1] +
  frac_ja_dose1 * delay$jansen
delay_dose1 <- ifelse(delay_dose1 == 0, 1, delay_dose1) # this prevents from dividing by 0 in the ODEs
delay_dose2 <- frac_pf_dose2 * delay$pfizer[2] +
  frac_mo_dose2 * delay$moderna[1] +
  frac_az_dose2 * delay$astrazeneca[2] +
  frac_ja_dose2 * delay$jansen
delay_dose2 <- ifelse(delay_dose2 == 0, 1, delay_dose2) # this prevents from dividing by 0 in the ODEs

eta_vec <- 1 - ifelse(is.nan(comp_ve_dose1), 0, comp_ve_dose1)
names(eta_vec) <- paste0("eta_",1:9)
eta2_vec <- 1- ifelse(is.nan(comp_ve_dose2), 0, comp_ve_dose2)
names(eta2_vec) <- paste0("eta2_",1:9)

comp_h_mult_dose1 <- 1 - ifelse(is.nan(hosp_mult_dose1), 0, hosp_mult_dose1)
names(comp_h_mult_dose1) <- paste0("hosp_mult_",1:9)
comp_h_mult_dose2 <- 1- ifelse(is.nan(hosp_mult_dose2), 0, hosp_mult_dose2)
names(comp_h_mult_dose2) <- paste0("hosp_mult2_",1:9)


return(vac_schedule_orig_new)
}

# convert to number of vaccines
# age_dist_9 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
#               0.14514332, 0.12092904, 0.08807406, 0.04622194)
# n_vec_9 <- n * age_dist_9
# 
# num_vacs <- sweep(vac_schedule_orig_new[,-1], 2, n_vec_9, "*") %>%
#   mutate(date = vac_schedule_orig_new$date) %>%
#   select(date, pf_d1_1:Ja_d2_9)
# 
# # read in num vacs and then convert back to proportion
# num_o2y <- read_csv("inst/extdata/data/num_vac_basis_scenario_20210311_o2y.csv")
# num_y2o <- read_csv("inst/extdata/data/num_vac_basis_scenario_20210311_y2o.csv")
# num_alt <- read_csv("inst/extdata/data/num_vac_basis_scenario_20210311_alt.csv")
# num_no_vac_healthy <- read_csv("inst/extdata/data/num_vac_basis_scenario_20210311_no_vac_healthy.csv")
# 
# prop_o2y <- sweep(num_o2y[,-1], 2, n_vec_9, "/") %>%
#   mutate(date = vac_schedule_orig_new$date) %>%
#   select(date, pf_d1_1:Ja_d2_9)
# 
# prop_y2o <- sweep(num_y2o[,-1], 2, n_vec_9, "/") %>%
#   mutate(date = vac_schedule_orig_new$date) %>%
#   select(date, pf_d1_1:Ja_d2_9)
# 
# prop_alt <- sweep(num_alt[,-1], 2, n_vec_9, "/") %>%
#   mutate(date = vac_schedule_orig_new$date) %>%
#   select(date, pf_d1_1:Ja_d2_9)
# 
# prop_no_vac_healthy <- sweep(num_no_vac_healthy[,-1], 2, n_vec_9, "/") %>%
#   mutate(date = vac_schedule_orig_new$date) %>%
#   select(date, pf_d1_1:Ja_d2_9)

