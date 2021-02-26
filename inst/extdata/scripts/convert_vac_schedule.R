# convert cumulative vaccination schedule to non-cumulative ---------------------------------------
library(readr)
library(tidyr)
library(dplyr)

# read in vac schedule file -----------------------------------------------------------------------
cum_vac_schedule_orig <- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210225.csv")
cum_vac_schedule_no2dose<- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210225 scno2nddose.csv")
cum_vac_schedule_delay3mo <- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210225 sc2nddose3mo.csv")

# to combine age groups 9 and 10
age_dist_10 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
n <- 17407585 # Dutch population size
n_vec_10 <- n * age_dist_10

# take the difference for each row
# original
vac_schedule_orig <- data.frame(diff(as.matrix(cum_vac_schedule_orig[-1,-1]))) %>%
  add_row(cum_vac_schedule[1,-1],.before = 1) %>%
  mutate(date = seq.Date(from = as.Date("2021-01-04"), to = as.Date("2021-12-30"), by = 1),
         pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
         ) %>%
  select(date, pf_d1_1:az_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10) 

# no second dose
vac_schedule_no2dose <- data.frame(diff(as.matrix(cum_vac_schedule_no2dose[-1,-1]))) %>%
  add_row(cum_vac_schedule[1,-1],.before = 1) %>%
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
  add_row(cum_vac_schedule[1,-1],.before = 1) %>%
  mutate(date = seq.Date(from = as.Date("2021-01-04"), to = as.Date("2021-12-30"), by = 1),
         pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
  ) %>%
  select(date, pf_d1_1:az_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10) %>%
  filter(date < as.Date("2021-03-08"))


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
  filter(date > as.Date("2021-01-31"))  %>%
  add_row(before_feb, .before = 1)
vac_schedule_delay3mo_new <- bind_rows(vac_schedule_orig_8March, vac_schedule_delay3mo) %>%
  filter(date > as.Date("2021-01-31"))  %>%
  add_row(before_feb, .before = 1)

vac_schedule <- vac_schedule_orig %>%
  filter(date > as.Date("2021-02-01"))%>%
  add_row(before_feb, .before = 1) 


#%>%
  # mutate(pf_d1_9 = ifelse(pf_d1_9 == 0.0001, 0, pf_d1_9),
  #        pf_d2_9 = ifelse(pf_d2_9 == 0.0001, 0, pf_d2_9),
  #        mo_d1_9 = ifelse(mo_d1_9 == 0.0001, 0, mo_d1_9),
  #        mo_d2_9 = ifelse(mo_d2_9 == 0.0001, 0, mo_d2_9),
  #        az_d1_9 = ifelse(az_d1_9 == 0.0001, 0, az_d1_9),
  #        az_d2_9 = ifelse(az_d2_9 == 0.0001, 0, az_d2_9))

# re-allocate all second doses beyond 2 March 2021 to first doses
# this results in some age groups with >100%
# vac_schedule_no_second_doses <- vac_schedule %>%
#   mutate(# pfizer dose 1
#          pf_d1_1 = ifelse(date > as.Date("2021-03-08"), pf_d1_1 + pf_d2_1, pf_d1_1),
#          pf_d1_2 = ifelse(date > as.Date("2021-03-08"), pf_d1_2 + pf_d2_2, pf_d1_2),
#          pf_d1_3 = ifelse(date > as.Date("2021-03-08"), pf_d1_3 + pf_d2_3, pf_d1_3),
#          pf_d1_4 = ifelse(date > as.Date("2020-03-08"), pf_d1_4 + pf_d2_4, pf_d1_4),
#          pf_d1_5 = ifelse(date > as.Date("2020-03-08"), pf_d1_5 + pf_d2_5, pf_d1_5),
#          pf_d1_6 = ifelse(date > as.Date("2020-03-08"), pf_d1_6 + pf_d2_6, pf_d1_6),
#          pf_d1_7 = ifelse(date > as.Date("2020-03-08"), pf_d1_7 + pf_d2_7, pf_d1_7),
#          pf_d1_8 = ifelse(date > as.Date("2020-03-08"), pf_d1_8 + pf_d2_8, pf_d1_8),
#          pf_d1_9 = ifelse(date > as.Date("2020-03-08"), pf_d1_9 + pf_d2_9, pf_d1_9),
#          # pfizer dose 2 (make zero)
#          pf_d2_1 = ifelse(date > as.Date("2021-03-08"), 0, pf_d2_1),
#          pf_d2_2 = ifelse(date > as.Date("2021-03-08"), 0, pf_d2_2),
#          pf_d2_3 = ifelse(date > as.Date("2021-03-08"), 0, pf_d2_3),
#          pf_d2_4 = ifelse(date > as.Date("2020-03-08"), 0, pf_d2_4),
#          pf_d2_5 = ifelse(date > as.Date("2020-03-08"), 0, pf_d2_5),
#          pf_d2_6 = ifelse(date > as.Date("2020-03-08"), 0, pf_d2_6),
#          pf_d2_7 = ifelse(date > as.Date("2020-03-08"), 0, pf_d2_7),
#          pf_d2_8 = ifelse(date > as.Date("2020-03-08"), 0, pf_d2_8),
#          pf_d2_9 = ifelse(date > as.Date("2020-03-08"), 0, pf_d2_9),
#          # moderna dose 1
#          mo_d1_1 = ifelse(date > as.Date("2021-03-08"), mo_d1_1 + mo_d2_1, mo_d1_1),
#          mo_d1_2 = ifelse(date > as.Date("2021-03-08"), mo_d1_2 + mo_d2_2, mo_d1_2),
#          mo_d1_3 = ifelse(date > as.Date("2021-03-08"), mo_d1_3 + mo_d2_3, mo_d1_3),
#          mo_d1_4 = ifelse(date > as.Date("2020-03-08"), mo_d1_4 + mo_d2_4, mo_d1_4),
#          mo_d1_5 = ifelse(date > as.Date("2020-03-08"), mo_d1_5 + mo_d2_5, mo_d1_5),
#          mo_d1_6 = ifelse(date > as.Date("2020-03-08"), mo_d1_6 + mo_d2_6, mo_d1_6),
#          mo_d1_7 = ifelse(date > as.Date("2020-03-08"), mo_d1_7 + mo_d2_7, mo_d1_7),
#          mo_d1_8 = ifelse(date > as.Date("2020-03-08"), mo_d1_8 + mo_d2_8, mo_d1_8),
#          mo_d1_9 = ifelse(date > as.Date("2020-03-08"), mo_d1_9 + mo_d2_9, mo_d1_9),
#          # moderna dose 2 (make zero)
#          mo_d2_1 = ifelse(date > as.Date("2021-03-08"), 0, mo_d2_1),
#          mo_d2_2 = ifelse(date > as.Date("2021-03-08"), 0, mo_d2_2),
#          mo_d2_3 = ifelse(date > as.Date("2021-03-08"), 0, mo_d2_3),
#          mo_d2_4 = ifelse(date > as.Date("2020-03-08"), 0, mo_d2_4),
#          mo_d2_5 = ifelse(date > as.Date("2020-03-08"), 0, mo_d2_5),
#          mo_d2_6 = ifelse(date > as.Date("2020-03-08"), 0, mo_d2_6),
#          mo_d2_7 = ifelse(date > as.Date("2020-03-08"), 0, mo_d2_7),
#          mo_d2_8 = ifelse(date > as.Date("2020-03-08"), 0, mo_d2_8),
#          mo_d2_9 = ifelse(date > as.Date("2020-03-08"), 0, mo_d2_9),
#          # astrazeneca dose 1
#          az_d1_1 = ifelse(date > as.Date("2021-03-08"), az_d1_1 + az_d2_1, az_d1_1),
#          az_d1_2 = ifelse(date > as.Date("2021-03-08"), az_d1_2 + az_d2_2, az_d1_2),
#          az_d1_3 = ifelse(date > as.Date("2021-03-08"), az_d1_3 + az_d2_3, az_d1_3),
#          az_d1_4 = ifelse(date > as.Date("2020-03-08"), az_d1_4 + az_d2_4, az_d1_4),
#          az_d1_5 = ifelse(date > as.Date("2020-03-08"), az_d1_5 + az_d2_5, az_d1_5),
#          az_d1_6 = ifelse(date > as.Date("2020-03-08"), az_d1_6 + az_d2_6, az_d1_6),
#          az_d1_7 = ifelse(date > as.Date("2020-03-08"), az_d1_7 + az_d2_7, az_d1_7),
#          az_d1_8 = ifelse(date > as.Date("2020-03-08"), az_d1_8 + az_d2_8, az_d1_8),
#          az_d1_9 = ifelse(date > as.Date("2020-03-08"), az_d1_9 + az_d2_9, az_d1_9),
#          # astrazeneca dose 2 (make zero)
#          az_d2_1 = ifelse(date > as.Date("2021-03-08"), 0, az_d2_1),
#          az_d2_2 = ifelse(date > as.Date("2021-03-08"), 0, az_d2_2),
#          az_d2_3 = ifelse(date > as.Date("2021-03-08"), 0, az_d2_3),
#          az_d2_4 = ifelse(date > as.Date("2020-03-08"), 0, az_d2_4),
#          az_d2_5 = ifelse(date > as.Date("2020-03-08"), 0, az_d2_5),
#          az_d2_6 = ifelse(date > as.Date("2020-03-08"), 0, az_d2_6),
#          az_d2_7 = ifelse(date > as.Date("2020-03-08"), 0, az_d2_7),
#          az_d2_8 = ifelse(date > as.Date("2020-03-08"), 0, az_d2_8),
#          az_d2_9 = ifelse(date > as.Date("2020-03-08"), 0, az_d2_9)
#          )
# 
# 
# # for loop approach
# dose2_stop_date <- as.Date("2021-03-08")
# second_dose_scenario <- "none" # can also be "double" for doubling time in between
# dates <- vac_schedule$date
# vac_schedule2 <- vac_schedule[,-1] # we'll write over vac_schedule2
# for (i in 1:dim(vac_schedule)[1]){
#   if(dates[i] < dose2_stop_date){next}
#   # check cumulative coverage so far
#   cumulative_sum <- vac_schedule %>%
#     filter(date <= dates[i]) %>%
#     select(-date) %>%
#     summarise_all(sum)
#   
#   pf_dose1 <- cumulative_sum %>% select(pf_d1_1:pf_d1_9)
#   pf_dose2 <- cumulative_sum %>% select(pf_d2_1:pf_d2_9)
#   mo_dose1 <- cumulative_sum %>% select(mo_d1_1:mo_d1_9)
#   mo_dose2 <- cumulative_sum %>% select(mo_d2_1:mo_d2_9)
#   az_dose1 <- cumulative_sum %>% select(az_d1_1:az_d1_9)
#   az_dose2 <- cumulative_sum %>% select(az_d2_1:az_d2_9)
#   
#   dose1_all <- pf_dose1 + mo_dose1 + az_dose1
#   # loop over age groups
#   for (j in 1:9){
#     # check if vac converage is < 1
#     if(pf_dose1[j] < 1){
#       # pfizer
#       vac_schedule2[i,j] <- vac_schedule2[i,j] + vac_schedule2[i,j + 9]
#       vac_schedule2[i,j + 9] <- 0
#       # moderna
#       if(pf_dose1[j] + mo_dose1[j])
#       vac_schedule2[i,j + 18] <- vac_schedule2[i,j + 18] + vac_schedule2[i,j + 27]
#       vac_schedule2[i,j + 27] <- 0
#       #astrazeneca
#       vac_schedule2[i,j + 36] <- vac_schedule2[i,j + 36] + vac_schedule2[i,j + 45]
#       vac_schedule2[i,j + 45] <- 0
#       
#     } else if (dose1_all[j] >= 1){
#       # pfizer
#       new_perc_pf <- (vac_schedule2[i,j + 9] * n_vec[j])/n_vec[j-1]
#       vac_schedule2[i,j-1] <- vac_schedule2[i,j-1] + new_perc_pf
#       vac_schedule2[i,j + 9] <- 0
#       # moderna
#       new_perc_mo <- (vac_schedule2[i,j + 27] * n_vec[j])/n_vec[j-1]
#       vac_schedule2[i,(j + 18) - 1] <- vac_schedule2[i,(j + 18) - 1] + new_perc_mo
#       vac_schedule2[i,j + 27] <- 0
#       #astrazeneca
#       new_perc_az <- (vac_schedule2[i,j + 45] * n_vec[j])/n_vec[j-1]
#       vac_schedule2[i,(j + 36) - 1] <- vac_schedule2[i,(j + 36) - 1] + new_perc_az
#       vac_schedule2[i,j + 45] <- 0
#     }
#     
#     
#     
#   } # end loop over age groups
#   
# } # end loop over dates
# 
# # check cumulative sums of vac_schedule2
# cumulative_sum2 <- vac_schedule2 %>%
#   summarise_all(sum)













