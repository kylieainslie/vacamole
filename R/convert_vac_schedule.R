#' convert cumulative vaccination schedule to non-cumulative ----------------------------------
#' @param vac_schedule 
#' @param ve 
#' @param hosp_multiplier
#' @param delay
#' @param ve_trans
#' @param add_child_vac
#' @return 
#' @keywords vacamole
#' @export
convert_vac_schedule <- function(vac_schedule, 
                                 ve, 
                                 hosp_multiplier, 
                                 delay, 
                                 ve_trans,
                                 add_child_vac = FALSE){
# to combine age groups 9 and 10 --------------------------------------------------------------
age_dist_10 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
n <- 17407585 # Dutch population size
n_vec_10 <- n * age_dist_10

# take the difference for each row ------------------------------------------------------------
vac_schedule_orig <- data.frame(diff(as.matrix(vac_schedule[-1,-1]))) %>%
  add_row(vac_schedule[1,-1],.before = 1) %>%
  mutate(date = seq.Date(from = as.Date("2021-01-04"), to = as.Date("2021-12-30"), by = 1),
         pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d1_9 = (ja_d1_9 * n_vec_10[9] + ja_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d2_9 = (ja_d2_9 * n_vec_10[9] + ja_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
         ) %>%
  select(date, pf_d1_1:ja_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10, 
         -az_d2_10, -ja_d1_10, -ja_d2_10) 

# filter for dates before 1 February 2021 -----------------------------------------------------
before_feb <- vac_schedule_orig %>%
  filter(date < as.Date("2021-02-01")) %>%
  select(-date) %>%
  summarise_all(sum) %>%
  mutate(date = as.Date("2021-01-31")) %>%
  select(date, pf_d1_1:ja_d2_9)

# filter so start date is 1 Feb 2021 with row for Jan 31 to reflect 
# people who have already been vaccinated -----------------------------------------------------
vac_schedule_orig_new <- vac_schedule_orig %>%
  filter(date > as.Date("2021-01-31")) %>%
  add_row(before_feb, .before = 1)

# add vaccination of children 12 - 17 starting September 1, 2021
# assume 50,000 vaccines per day for 21 days to get 85% coverage in this age group
# this equates to an increase in vaccination coverage in this age group of 
# 0.0243 per day for a total of 0.51 after 21 days.
# we assume the second dose is given 6 weeks after the first dose
if(add_child_vac){
  p_child_12_17_doses <- 0.6 * 0.85
  n_child_12_17_doses <- n_vec_10[2] * p_child_12_17_doses
  days_child_12_17 <- n_child_12_17_doses/50000
  p_child_12_17_doses_per_day <- p_child_12_17_doses/days_child_12_17

  vac_schedule_orig_new <- vac_schedule_orig_new %>%
    mutate(pf_d1_2 = ifelse(date >= as.Date("2021-09-01") & date <= as.Date("2021-09-21"), pf_d1_2 + p_child_12_17_doses_per_day, pf_d1_2),
          pf_d2_2 = ifelse(date >= as.Date("2021-10-13") & date <= as.Date("2021-11-02"), pf_d2_2 + p_child_12_17_doses_per_day, pf_d2_2))

  vac_schedule_new_cs <- cumsum(vac_schedule_orig_new[,-1])
} else {
  vac_schedule_new <- vac_schedule_orig_new
  vac_schedule_new_cs <- cumsum(vac_schedule_orig_new[,-1])
}
# create separate data frame for each vaccine -------------------------------------------------
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
ja_dose1 <- vac_schedule_orig_new %>% select(date, ja_d1_1:ja_d1_9)
ja_dose1_cs <- vac_schedule_new_cs %>% select(ja_d1_1:ja_d1_9)
ja_dose2 <- vac_schedule_orig_new %>% select(date, ja_d2_1:ja_d2_9)
ja_dose2_cs <- vac_schedule_new_cs %>% select(ja_d2_1:ja_d2_9)

# calculate composite VE ---------------------------------------------------------------------
name_suffix_d1 <- c(substr(names(pf_dose1[,-1]), 3, 7))
name_suffix_d2 <- c(substr(names(pf_dose2[,-1]), 3, 7))

# daily vaccination rate
alpha_dose1 <- pf_dose1[,-1] + mo_dose1[,-1] + az_dose1[,-1] + ja_dose1[,-1]
names(alpha_dose1) <- paste0("alpha", name_suffix_d1)
alpha_dose2 <- pf_dose2[,-1] + mo_dose2[,-1] + az_dose2[,-1] + ja_dose2[,-1]
names(alpha_dose2) <- paste0("alpha", name_suffix_d2)

# cumulative vaccination percentage
total_dose1 <- pf_dose1_cs + mo_dose1_cs + az_dose1_cs + ja_dose1_cs
names(total_dose1) <- paste0("tot", name_suffix_d1)
total_dose2 <- pf_dose2_cs + mo_dose2_cs + az_dose2_cs + ja_dose2_cs
names(total_dose2) <- paste0("tot", name_suffix_d2)

# fraction of each vaccine given up to current time point
frac_pf_dose1 <- data.matrix(pf_dose1_cs/total_dose1)
frac_pf_dose1 <- ifelse(is.nan(frac_pf_dose1), 0, frac_pf_dose1)
frac_pf_dose2 <- data.matrix(pf_dose2_cs/total_dose2)
frac_pf_dose2 <- ifelse(is.nan(frac_pf_dose2), 0, frac_pf_dose2)

frac_mo_dose1 <- data.matrix(mo_dose1_cs/total_dose1)
frac_mo_dose1 <- ifelse(is.nan(frac_mo_dose1), 0, frac_mo_dose1)
frac_mo_dose2 <- data.matrix(mo_dose2_cs/total_dose2)
frac_mo_dose2 <- ifelse(is.nan(frac_mo_dose2), 0, frac_mo_dose2)

frac_az_dose1 <- data.matrix(az_dose1_cs/total_dose1)
frac_az_dose1 <- ifelse(is.nan(frac_az_dose1), 0, frac_az_dose1)
frac_az_dose2 <- data.matrix(az_dose2_cs/total_dose2)
frac_az_dose2 <- ifelse(is.nan(frac_az_dose2), 0, frac_az_dose2)

frac_ja_dose1 <- data.matrix(ja_dose1_cs/total_dose1)
frac_ja_dose1 <- ifelse(is.nan(frac_ja_dose1), 0, frac_az_dose1)
frac_ja_dose2 <- data.matrix(ja_dose2_cs/total_dose2)
frac_ja_dose2 <- ifelse(is.nan(frac_ja_dose2), 0, frac_ja_dose2)

# composite VE (against infection)
comp_ve_dose1 <- frac_pf_dose1 * ve$pfizer[1] + 
  frac_mo_dose1 * ve$moderna[1] + 
  frac_az_dose1 * ve$astrazeneca[1] +
  frac_ja_dose1 * ve$jansen
colnames(comp_ve_dose1) <- paste0("ve", name_suffix_d1)
comp_ve_dose2 <- frac_pf_dose2 * ve$pfizer[2] + 
  frac_mo_dose2 * ve$moderna[2] + 
  frac_az_dose2 * ve$astrazeneca[2] +
  frac_ja_dose2 * ve$jansen
colnames(comp_ve_dose2) <- paste0("ve", name_suffix_d2)

# eta
eta_dose1 <- 1 - comp_ve_dose1
eta_dose2 <- 1- comp_ve_dose2

# rate of hospitalisations multiplier
hosp_mult_dose1 <- frac_pf_dose1 * hosp_multiplier$pfizer[1] +
  frac_mo_dose1 * hosp_multiplier$moderna[1] +
  frac_az_dose1 * hosp_multiplier$astrazeneca[1] +
  frac_ja_dose1 * hosp_multiplier$jansen
colnames(hosp_mult_dose1) <- paste0("hosp_mult", name_suffix_d1)
hosp_mult_dose2 <- frac_pf_dose2 * hosp_multiplier$pfizer[2] +
  frac_mo_dose2 * hosp_multiplier$moderna[2] +
  frac_az_dose2 * hosp_multiplier$astrazeneca[2] +
  frac_ja_dose2 * hosp_multiplier$jansen
colnames(hosp_mult_dose2) <- paste0("hosp_mult", name_suffix_d1)

eta_hosp_dose1 <- 1 - hosp_mult_dose1
eta_hosp_dose2 <- 1 - hosp_mult_dose2

# composite VE (against transmission)
comp_ve_trans_dose1 <- frac_pf_dose1 * ve_trans$pfizer[1] + 
  frac_mo_dose1 * ve_trans$moderna[1] + 
  frac_az_dose1 * ve_trans$astrazeneca[1] +
  frac_ja_dose1 * ve_trans$jansen
colnames(comp_ve_trans_dose1) <- paste0("ve_trans", name_suffix_d1)
comp_ve_trans_dose2 <- frac_pf_dose2 * ve_trans$pfizer[2] + 
  frac_mo_dose2 * ve_trans$moderna[2] + 
  frac_az_dose2 * ve_trans$astrazeneca[2] +
  frac_ja_dose2 * ve_trans$jansen
colnames(comp_ve_trans_dose2) <- paste0("ve_trans", name_suffix_d2)

# eta_trans
eta_trans_dose1 <- 1 - comp_ve_trans_dose1
eta_trans_dose2 <- 1- comp_ve_trans_dose2


# composite delay to protection
delay_dose1 <- frac_pf_dose1 * delay$pfizer[1] +
  frac_mo_dose1 * delay$moderna[1] +
  frac_az_dose1 * delay$astrazeneca[1] +
  frac_ja_dose1 * delay$jansen
delay_dose1 <- ifelse(delay_dose1 == 0, 1, delay_dose1) # this prevents from dividing by 0 in the ODEs
colnames(delay_dose1) <- paste0("delay", name_suffix_d1)
delay_dose2 <- frac_pf_dose2 * delay$pfizer[2] +
  frac_mo_dose2 * delay$moderna[1] +
  frac_az_dose2 * delay$astrazeneca[2] +
  frac_ja_dose2 * delay$jansen
delay_dose2 <- ifelse(delay_dose2 == 0, 1, delay_dose2) # this prevents from dividing by 0 in the ODEs
colnames(delay_dose2) <- paste0("delay", name_suffix_d2)

rtn <- list(alpha_dose1 = alpha_dose1, 
            alpha_dose2 = alpha_dose2, 
            eta_dose1 = eta_dose1, 
            eta_dose2 = eta_dose2,
            delay_dose1 = delay_dose1,
            delay_dose2 = delay_dose2,
            eta_hosp_dose1 = eta_hosp_dose1, 
            eta_hosp_dose2 = eta_hosp_dose2,
            eta_trans_dose1 = eta_trans_dose1, 
            eta_trans_dose2 = eta_trans_dose2)

return(rtn)
}
