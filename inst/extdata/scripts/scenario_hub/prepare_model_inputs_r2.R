# Script for preparing and specifying model inputs for running scenarios for the 
# European Scenario Hub
# URL: https://github.com/covid19-forecast-hub-europe/covid19-scenario-hub-europe#readme

# preamble ---------------------------------------------------------------------
# This script will specify all of the model inputs and (if necessary)
# write them as RDS files, which will be called in main_script.R
# ------------------------------------------------------------------------------

library(fst)
library(tidyverse)
# get file path ----------------------------------------------------------------
path <- "/rivm/r/COVID-19/Surveillance/Data/OSIRIS/Geschoond/"
file <- list.files(path, pattern = ".fst")
if (identical(file, character(0))) {
  path <- paste0(path,"Previous/")
  file <- list.files(path, pattern = ".fst") %>%
    max()
}

# read in data -----------------------------------------------------------------
osiris <- read_fst(paste0(path,file))

# data wrangle -----------------------------------------------------------------
osiris_tally <- osiris %>%           # aggregate for number of cases per day
  # this removes any identifiable data
  select(OSIRISNR, INFECTIEZIEKTE, ZIE1eZiekteDt, Land) %>%
  filter(Land == "Nederland",
         INFECTIEZIEKTE %in% c("NCOV", "Weak Positive", 
                               "Antitgen Pos. + Symptoms", 
                               "PCR Positief", "Antigen Positief")) %>%
  select(-Land) %>%
  rename(date = ZIE1eZiekteDt) %>%
  group_by(date) %>%
  summarise(inc = n()) %>%
  filter(!is.na(date)) %>%
  complete(date = seq.Date(min(date), max(date), by="day"), fill = list(inc = 0))

cutoff_date <- as.Date("2022-07-24")

osiris1 <- osiris_tally %>%
  filter(date <= cutoff_date)

# or read from directory
osiris1 <- readRDS("inst/extdata/data/case_data_upto_2022-07-24.rds")

# Vaccination schedule ---------------------------------------------------------
# Use the following code on the file direct from Pieter
# Read in file and change column names for booster doses
# vac_path <- "C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/manuscripts/impact_vac/data/vaccination_scenarios/"
vac_path <- "/rivm/s/ainsliek/data/"

vac_sched <- read_csv(paste0(vac_path,"Cum_upt20220714.csv")) %>%
  rename_with(~ gsub("B1", "d3", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("B2", "d4", .x, fixed = TRUE)) %>%
  select(-starts_with("X")) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

# add columns for 5th dose of pfizer and moderna vaccines
new_columns <- c(paste0("pf_d5_", 1:10), paste0("mo_d5_", 1:10))
vac_sched[,new_columns] <- 0

# rearrange columns to preserve correct order (also exclude Novovax doses)
vac_sched1 <- vac_sched %>%
  select(date, pf_d1_1:pf_d4_10, pf_d5_1:pf_d5_10,
         mo_d1_1:mo_d4_10, mo_d5_1:mo_d5_10,
         az_d1_1:az_d2_10, ja_d1_1:ja_d2_10)

# create "empty" values from 1/1/2020 until start of vac sched (1/4/2021)
n_cols <- dim(vac_sched1)[2]-1 # exclude date column
n_rows <- vac_sched1$date[1] - osiris1$date[1]
empty_mat <- matrix(rep(0, n_cols * n_rows), nrow = n_rows)
dates <- seq.Date(osiris1$date[1], vac_sched1$date[1]-1, by = "day")
my_df <- data.frame(date = dates, empty_mat)
names(my_df) <- names(vac_sched1)
vac_schedule <- bind_rows(my_df, vac_sched1) %>%
  filter(date < as.Date("2022-07-10"))

# write out to directory
write.csv(vac_schedule,"inst/extdata/inputs/vaccination_schedules/vac_schedule_real_20220709.csv")
saveRDS(vac_schedule,"inst/extdata/inputs/vaccination_schedules/vac_schedule_real_20220709.rds")

# ------------------------------------------------------------------------------

# No booster campaign ----------------------------------------------------------
vac_schedule <- readRDS("inst/extdata/inputs/vaccination_schedules/vac_schedule_real_20220709.rds")
extra_start_date <- tail(vac_schedule$date,1) + 1
extra_end_date <- as.Date("2023-07-06")
extra_dates <- seq.Date(from = as.Date(extra_start_date), 
                        to = as.Date(extra_end_date), by = 1)
vac_schedule_no_boost <- data.frame(date = extra_dates) %>%
  full_join(vac_schedule, ., by = "date") %>%
  fill(-.data$date)

saveRDS(vac_schedule_no_boost, 
        "inst/extdata/inputs/vaccination_schedules/vac_schedule_scenario_hub_round2_no_boost.rds")

# Update vac schedule for Round 2 ----------------------------------------------
# Booster campaign in 60+ (only) -----------------------------------------------
# 1) increase vaccination of 4th dose to achieve 50% coverage by 15 September
# 2) begin 5th dose administration on 15 September and reach 50% coverage by 
# 15 December

# 1) increase coverage of 4th dose
vac_schedule <- readRDS("inst/extdata/inputs/vaccination_schedules/vac_schedule_real_20220709.rds")
current_4d_prop <- vac_schedule %>% 
  tail(.,1) %>%
  select(date, pf_d4_7:pf_d4_9, mo_d4_7:mo_d4_9) %>%
  mutate(d4_7 = pf_d4_7 + mo_d4_7,
         d4_8 = pf_d4_8 + mo_d4_8,
         d4_9 = pf_d4_9 + mo_d4_9) 
to_get_to_50 <- as.numeric(c(rep(0.5,3)) - current_4d_prop[c("d4_7", "d4_8", "d4_9")])
n_days <- as.Date("2022-09-15") - current_4d_prop$date
daily_prop <- ifelse(to_get_to_50 < 0, 0, to_get_to_50) / (as.numeric(n_days))
end_mo_prop <- 0.5 - (current_4d_prop$mo_d4_7 + current_4d_prop$pf_d4_7)
# sequence of increasing vac coverage
vac_cov_vec <- seq(from = tail(vac_schedule$mo_d4_7,1) + daily_prop[1], 
                   to = tail(vac_schedule$mo_d4_7,1) + end_mo_prop, 
                   length.out = n_days)

# add 4th dose vaccinations to get to 50% by 15 September 2022
# assume all 4th doses are moderna
extra_start_date <- tail(vac_schedule$date,1) + 1
extra_end_date <- as.Date("2022-09-15")
extra_dates <- seq.Date(from = as.Date(extra_start_date), 
                        to = as.Date(extra_end_date), by = 1)
vac_schedule_4d <- data.frame(date = extra_dates) %>%
  full_join(vac_schedule, ., by = "date") %>%
  fill(-.data$date)

vac_schedule_4d$mo_d4_7[which(vac_schedule_4d$date %in% extra_start_date:extra_end_date)] <- vac_cov_vec
# add more extra dates
extra_dates2 <- seq.Date(from = as.Date(extra_end_date)+1, 
                         to = as.Date("2023-07-29"), by = 1)
vac_schedule_4da <- data.frame(date = extra_dates2) %>%
  full_join(vac_schedule_4d, ., by = "date") %>%
  fill(-.data$date)

#saveRDS(vac_schedule_4da, "inst/extdata/inputs/vaccination_schedules/vac_schedule_scenario_hub_round1_AC.rds")

# 2) add 5th doses 
# assume 5th doses distribution is 25% pfizer and 75% moderna (up to 50% coverage)
# 5th doses start 15 September and end 15 December
n_days_5d <- as.numeric(as.Date("2022-12-15") - as.Date("2022-09-15"))
vac_cov_vec_pf_5d <- seq(from = 0, to = (0.5*0.25), length.out = n_days_5d)
vac_cov_vec_mo_5d <- seq(from = 0, to = (0.5*0.75), length.out = n_days_5d)

extra_start_date_5d <- tail(vac_schedule_4d$date,1) + 1
extra_end_date_5d <- as.Date("2022-12-15")
extra_dates_5d <- seq.Date(from = as.Date(extra_start_date_5d), 
                        to = as.Date(extra_end_date_5d), by = 1)
vac_schedule_5d <- data.frame(date = extra_dates_5d) %>%
  full_join(vac_schedule_4d, ., by = "date") %>%
  fill(-.data$date)
# moderna 5th doses
vac_schedule_5d$mo_d5_7[which(vac_schedule_5d$date %in% extra_start_date_5d:extra_end_date_5d)] <- vac_cov_vec_mo_5d
vac_schedule_5d$mo_d5_8[which(vac_schedule_5d$date %in% extra_start_date_5d:extra_end_date_5d)] <- vac_cov_vec_mo_5d
vac_schedule_5d$mo_d5_9[which(vac_schedule_5d$date %in% extra_start_date_5d:extra_end_date_5d)] <- vac_cov_vec_mo_5d
# pfizer 5th doses
vac_schedule_5d$pf_d5_7[which(vac_schedule_5d$date %in% extra_start_date_5d:extra_end_date_5d)] <- vac_cov_vec_pf_5d
vac_schedule_5d$pf_d5_8[which(vac_schedule_5d$date %in% extra_start_date_5d:extra_end_date_5d)] <- vac_cov_vec_pf_5d
vac_schedule_5d$pf_d5_9[which(vac_schedule_5d$date %in% extra_start_date_5d:extra_end_date_5d)] <- vac_cov_vec_pf_5d
# add more extra dates
extra_dates_5d2 <- seq.Date(from = as.Date(extra_end_date_5d)+1, 
                         to = as.Date("2023-07-29"), by = 1)
vac_schedule_5da <- data.frame(date = extra_dates_5d2) %>%
  full_join(vac_schedule_5d, ., by = "date") %>%
  fill(-.data$date)

saveRDS(vac_schedule_5da, "inst/extdata/inputs/vaccination_schedules/vac_schedule_scenario_hub_round2_AC.rds")
# ------------------------------------------------------------------------------
# Booster campaign in 18+ ------------------------------------------------------
# 1) increase vaccination of 3rd dose to achieve 50% coverage by 15 September
# 2) begin 4th dose administration on 15 September and reach 50% coverage by 
# 15 December

# 1) increase coverage of 3rd dose in 18-59 year olds
# use scenarios A/C input file as a start (so we only have to deal with 18-59 
# year olds)
vac_schedule <- readRDS("inst/extdata/inputs/vaccination_schedules/vac_schedule_scenario_hub_round2_AC.rds")
current_3d_prop <- vac_schedule %>% 
  filter(date <= as.Date("2022-07-09")) %>%
  tail(.,1) %>%
  select(date, pf_d3_2:pf_d3_6, mo_d3_2:mo_d3_6) %>%
  mutate(d3_2 = pf_d3_2 + mo_d3_2,
         d3_3 = pf_d3_3 + mo_d3_3,
         d3_4 = pf_d3_4 + mo_d3_4,
         d3_5 = pf_d3_5 + mo_d3_5,
         d3_6 = pf_d3_6 + mo_d3_6) 
to_get_to_50 <- as.numeric(c(0.5*0.2,rep(0.5,4)) - current_3d_prop[,12:16])
n_days <- as.Date("2022-09-15") - current_3d_prop$date
daily_prop <- ifelse(to_get_to_50 < 0, 0, to_get_to_50) / (as.numeric(n_days))

# sequence of increasing vac coverage
vac_cov_vec_2 <- seq(from = tail(vac_schedule$mo_d3_2,1) + daily_prop[1], 
                     to = tail(vac_schedule$mo_d3_2,1) + to_get_to_50[1], 
                     length.out = n_days)
vac_cov_vec_3 <- seq(from = tail(vac_schedule$mo_d3_3,1) + daily_prop[2], 
                     to = tail(vac_schedule$mo_d3_3,1) + to_get_to_50[2], 
                     length.out = n_days)
vac_cov_vec_4 <- seq(from = tail(vac_schedule$mo_d3_4,1) + daily_prop[3], 
                     to = tail(vac_schedule$mo_d3_4,1) + to_get_to_50[3], 
                     length.out = n_days)

# add 3rd dose vaccinations to get to 50% by 15 September 2022
# assume all remaining 3rd doses are moderna
extra_start_date <- as.Date("2022-07-09") + 1
extra_end_date <- as.Date("2022-09-15")
# extra_dates <- seq.Date(from = as.Date(extra_start_date), 
#                         to = as.Date(extra_end_date), by = 1)
# vac_schedule_4d <- data.frame(date = extra_dates) %>%
#   full_join(vac_schedule, ., by = "date") %>%
#   fill(-.data$date)

vac_schedule$mo_d3_2[which(vac_schedule$date %in% extra_start_date:extra_end_date)] <- vac_cov_vec_2
vac_schedule$mo_d3_3[which(vac_schedule$date %in% extra_start_date:extra_end_date)] <- vac_cov_vec_3
vac_schedule$mo_d3_4[which(vac_schedule$date %in% extra_start_date:extra_end_date)] <- vac_cov_vec_4

# fill remaining rows with coverage from 15 September 2022
vac_schedule_3d <- vac_schedule %>%
  mutate(mo_d3_2 = ifelse(date > extra_end_date, NA, mo_d3_2),
         mo_d3_3 = ifelse(date > extra_end_date, NA, mo_d3_3),
         mo_d3_4 = ifelse(date > extra_end_date, NA, mo_d3_4)
  ) %>%
  fill(.data$mo_d3_2:.data$mo_d3_4)

#saveRDS(vac_schedule_4da, "inst/extdata/inputs/vaccination_schedules/vac_schedule_scenario_hub_round1_AC.rds")

# 2) add 4th doses 
# assume 4th doses are all pfizer (moderna and pfizer have same characteristics) 
# 5th doses start 15 September and end 15 December
extra_start_date_4d <- as.Date("2022-09-15") + 1
extra_end_date_4d <- as.Date("2022-12-15")
n_days_4d <- as.numeric(as.Date("2022-12-15") - as.Date("2022-09-15"))

# create dataframe for vac coverage vectors
cov_dat <- data.frame(date = seq.Date(extra_start_date_4d, extra_end_date_4d, 
                                      by = 1),
                      vac_cov_d4_2 = seq(from = tail(vac_schedule_3d$pf_d4_2,1), 
                                         to = (0.5*0.2) - tail(vac_schedule_3d$mo_d4_2,1), 
                                         length.out = n_days_4d),
                      vac_cov_d4_3 = seq(from = tail(vac_schedule_3d$pf_d4_3,1), 
                                         to = 0.5 - tail(vac_schedule_3d$mo_d4_3,1), 
                                         length.out = n_days_4d),
                      vac_cov_d4_4 = seq(from = tail(vac_schedule_3d$pf_d4_4,1), 
                                         to = 0.5 - tail(vac_schedule_3d$mo_d4_4,1), 
                                         length.out = n_days_4d),
                      vac_cov_d4_5 = seq(from = tail(vac_schedule_3d$pf_d4_5,1), 
                                         to = 0.5 - tail(vac_schedule_3d$mo_d4_5,1), 
                                         length.out = n_days_4d),
                      vac_cov_d4_6 = seq(from = tail(vac_schedule_3d$pf_d4_6,1), 
                                         to = 0.5 - tail(vac_schedule_3d$mo_d4_6,1), 
                                         length.out = n_days_4d)
                        )


# pfizer 4th doses
vac_schedule_3d$pf_d4_2[which(vac_schedule_3d$date %in% extra_start_date_4d:extra_end_date_4d)] <- cov_dat$vac_cov_d4_2
vac_schedule_3d$pf_d4_3[which(vac_schedule_3d$date %in% extra_start_date_4d:extra_end_date_4d)] <- cov_dat$vac_cov_d4_3
vac_schedule_3d$pf_d4_4[which(vac_schedule_3d$date %in% extra_start_date_4d:extra_end_date_4d)] <- cov_dat$vac_cov_d4_4
vac_schedule_3d$pf_d4_5[which(vac_schedule_3d$date %in% extra_start_date_4d:extra_end_date_4d)] <- cov_dat$vac_cov_d4_5
vac_schedule_3d$pf_d4_6[which(vac_schedule_3d$date %in% extra_start_date_4d:extra_end_date_4d)] <- cov_dat$vac_cov_d4_6
# add more extra dates
# extra_dates_5d2 <- seq.Date(from = as.Date(extra_end_date_5d)+1, 
#                             to = as.Date("2023-05-20"), by = 1)
vac_schedule1 <- vac_schedule_3d %>%
  mutate(pf_d4_2 = ifelse(date > extra_end_date_4d, NA, pf_d4_2),
         pf_d4_3 = ifelse(date > extra_end_date_4d, NA, pf_d4_3),
         pf_d4_4 = ifelse(date > extra_end_date_4d, NA, pf_d4_4),
         pf_d4_5 = ifelse(date > extra_end_date_4d, NA, pf_d4_5),
         pf_d4_6 = ifelse(date > extra_end_date_4d, NA, pf_d4_6)
         ) %>%
  fill(.data$pf_d4_2:.data$pf_d4_6)

saveRDS(vac_schedule1, "inst/extdata/inputs/vaccination_schedules/vac_schedule_scenario_hub_round2_BD.rds")
