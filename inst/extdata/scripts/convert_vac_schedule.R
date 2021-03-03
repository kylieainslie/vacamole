# convert cumulative vaccination schedule to non-cumulative ---------------------------------------
library(readr)
library(tidyr)
library(dplyr)

# read in vac schedule file -----------------------------------------------------------------------
cum_vac_schedule_orig <- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210225.csv")
cum_vac_schedule_no2dose<- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210225 scno2nddose.csv")
cum_vac_schedule_delay3mo <- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210225 sc2nddose3mo.csv")

#cum_vac_schedule_orig_w_jansen <- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210225 + Jansen.csv")
# to combine age groups 9 and 10
age_dist_10 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
n <- 17407585 # Dutch population size
n_vec_10 <- n * age_dist_10

# calculate vac coverage (up to 1 July 2021)

# take the difference for each row
# original
vac_schedule_orig <- data.frame(diff(as.matrix(cum_vac_schedule_orig[-1,-1]))) %>%
  add_row(cum_vac_schedule_orig[1,-1],.before = 1) %>%
  mutate(date = seq.Date(from = as.Date("2021-01-04"), to = as.Date("2021-12-30"), by = 1),
         pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         # ja_d1_9 = (Ja_d1_9 * n_vec_10[9] + Ja_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         # ja_d2_9 = (Ja_d2_9 * n_vec_10[9] + Ja_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
         ) %>%
  select(date, pf_d1_1:az_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10 #, 
         #-Ja_d1_10, -Ja_d2_10
         ) 

# no second dose
vac_schedule_no2dose <- data.frame(diff(as.matrix(cum_vac_schedule_no2dose[-1,-1]))) %>%
  add_row(cum_vac_schedule_no2dose[1,-1],.before = 1) %>%
  mutate(date = seq.Date(from = as.Date("2021-01-04"), to = as.Date("2021-12-30"), by = 1),
         pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
  ) %>%
  select(date, pf_d1_1:az_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10) %>%
  filter(date > as.Date("2021-03-08"))

# delay second dose for 3 months
vac_schedule_delay3mo <- data.frame(diff(as.matrix(cum_vac_schedule_delay3mo[-1,-1]))) %>%
  add_row(cum_vac_schedule_delay3mo[1,-1],.before = 1) %>%
  mutate(date = seq.Date(from = as.Date("2021-01-04"), to = as.Date("2021-12-30"), by = 1),
         pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
  ) %>%
  select(date, pf_d1_1:az_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10) %>%
  filter(date > as.Date("2021-03-08"))


# put the allocation up to 8 March from original schedule into new delay schedules
vac_schedule_orig_8March <- vac_schedule_orig %>%
  filter(date <= as.Date("2021-03-08"))

# filter for dates before 1 February 2021
before_feb <- vac_schedule_orig %>%
  filter(date < as.Date("2021-02-01")) %>%
  select(-date) %>%
  summarise_all(sum) %>%
  mutate(date = as.Date("2021-01-31")) %>%
  select(date, pf_d1_1:az_d2_9)

# filter so start date is 1 Feb 2021 with row for Jan 31 to reflect people who have already been vaccinated
vac_schedule_orig_new <- vac_schedule_orig %>%
  filter(date > as.Date("2021-01-31")) %>%
  add_row(before_feb, .before = 1)
vac_schedule_no2dose_new <- bind_rows(vac_schedule_orig_8March, vac_schedule_no2dose) %>%
  filter(date > as.Date("2021-01-31")) %>%
  add_row(before_feb, .before = 1)
vac_schedule_delay3mo_new <- bind_rows(vac_schedule_orig_8March, vac_schedule_delay3mo) %>%
  filter(date > as.Date("2021-01-31")) %>%
  add_row(before_feb, .before = 1)













