# convert cumulative vaccination schedule to non-cumulative ---------------------------------------
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
# read in vac schedule file -----------------------------------------------------------------------
cum_vac_sched <- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210311 Basis scenario.csv")

# to combine age groups 9 and 10
age_dist_10 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
n <- 17407585 # Dutch population size
n_vec_10 <- n * age_dist_10

# take the difference for each row
# original
vac_schedule_orig <- data.frame(diff(as.matrix(cum_vac_sched[-1,-1]))) %>%
  add_row(cum_vac_sched[1,-1],.before = 1) %>%
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

# convert to number of vaccines
age_dist_9 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.04622194)
n_vec_9 <- n * age_dist_9

num_vacs <- sweep(vac_schedule_orig_new[,-1], 2, n_vec_9, "*") %>%
  mutate(date = vac_schedule_orig_new$date) %>%
  select(date, pf_d1_1:Ja_d2_9)

# read in num vacs and then convert back to proportion
num_o2y <- read_csv("inst/extdata/data/num_vac_basis_scenario_20210311_o2y.csv")
num_y2o <- read_csv("inst/extdata/data/num_vac_basis_scenario_20210311_y2o.csv")
num_alt <- read_csv("inst/extdata/data/num_vac_basis_scenario_20210311_alt.csv")
num_no_vac_healthy <- read_csv("inst/extdata/data/num_vac_basis_scenario_20210311_no_vac_healthy.csv")

prop_o2y <- sweep(num_o2y[,-1], 2, n_vec_9, "/") %>%
  mutate(date = vac_schedule_orig_new$date) %>%
  select(date, pf_d1_1:Ja_d2_9)

prop_y2o <- sweep(num_y2o[,-1], 2, n_vec_9, "/") %>%
  mutate(date = vac_schedule_orig_new$date) %>%
  select(date, pf_d1_1:Ja_d2_9)

prop_alt <- sweep(num_alt[,-1], 2, n_vec_9, "/") %>%
  mutate(date = vac_schedule_orig_new$date) %>%
  select(date, pf_d1_1:Ja_d2_9)

prop_no_vac_healthy <- sweep(num_no_vac_healthy[,-1], 2, n_vec_9, "/") %>%
  mutate(date = vac_schedule_orig_new$date) %>%
  select(date, pf_d1_1:Ja_d2_9)

