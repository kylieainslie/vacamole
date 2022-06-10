#' convert cumulative vaccination schedule to non-cumulative -------------------
#' @param vac_schedule a data frame that the proportion of the population who 
#' receives vaccines
#' at each time point. Rows are time points, columns are 
#' <vaccine type>_<dose>_<age group>. For example,
#' the first column is the proportion of individuals in age group 1 who receive 
#' dose 1 of the first vaccine
#' type. The function assumes 9 or 10 age groups (1, .., 10) and four vaccine 
#' types (pf, mo, az, ja), each with a 2-dose regimen (d1, d2).
#' @param ve a named list of vaccine effectiveness against infection for each 
#' dose of each vaccine type.
#' @param hosp_multiplier a named list of the values for the vaccine 
#' effectiveness against hospitalization for each dose of each vaccine type. 
#' Vaccine effectiveness against hospitalization is incorporated as a multiplier
#' on the probability of being hospitalized after infection as 
#' (1 – VE_hospitalization) divided by (1-VE_infection) to account for the the 
#' inclusion of people who are never infected (and thus never hospitalized) 
#' included in the estimation of VE against hospitalization.
#' @param delay a named list of the time to protection for each dose of each 
#' vaccine type.
#' @param ve_trans a named list of vaccine effectiveness against transmission 
#' for each dose of each vaccine type.
#' @param wane logical, if TRUE vaccine effectiveness wanes by a logistic 
#' function parameterized by arguments
#' k and t0.
#' @param k logistic growth rate
#' @param t0 the time point at the midpoint of the logistic curve (where 50\%
#'  waning occurs)
#' @param add_extra_dates logical, if TRUE add extra rows to the vaccination 
#' schedule until extra_end_date
#' @param extra_end_date character string of the date (YYYY-MM-DD format) to 
#' end the vaccination schedule (which will
#' also be the last date of the simulation)
#' @return list of vaccination rate by day and dose and weighted VE and delay 
#' to protection by day and dose
#' @keywords vacamole
#' @import tidyr
#' @import dplyr
#' @export
convert_vac_schedule_debug <- function(vac_schedule,
                                  ve,
                                  hosp_multiplier,
                                  delay,
                                  ve_trans,
                                  wane = FALSE,
                                  add_extra_dates = FALSE,
                                  extra_start_date = "2022-01-01",
                                  extra_end_date = "2022-03-31"){
  
  # check if there are 9 age groups or 10 age groups ---------------------------
  # extract age group part of each relevant column name of vac_schedule
  names_check <- str_extract(names(vac_schedule), 
                             pattern = "(?<=[a-z]{2}_d[0-9]_)[0-9]+")
  # find the largest age group
  num_age_groups <- names_check %>%
    as.numeric() %>%
    max(na.rm = TRUE) %>%
    as.integer()
  
  # remove any redundant columns 
  vac_schedule <- vac_schedule %>% 
    select(date | matches("[a-z]{2}_d[0-9]_[0-9]+"))
  
  # size of each age group -----------------------------------------------------
  n <- 17407585 # Dutch population size
  n_vec_10 <- n * c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463,
                    0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
  
  # reformat date variable if necessary ----------------------------------------
  if (is.Date(vac_schedule$date)){ date_vec <- vac_schedule$date
  } else { date_vec <- as.Date(vac_schedule$date, format = "%m/%d/%Y")}
  
  # check that vac_schedule has 9 or 10 age groups -----------------------------
  `%!in%` <- Negate(`%in%`)
  if(num_age_groups %!in% c(9,10)){ 
    stop("number of age groups in vac_schedule is not 9 or 10, can't convert")
  }
  # if 10 age groups then: -----------------------------------------------------
  # combine age groups 9 and 10 into only one age group
  if (num_age_groups == 10) {
    vac_schedule <- vac_schedule %>%
      mutate(
        date = date_vec,
        pf_d1_9 = (.data$pf_d1_9 * n_vec_10[9] + .data$pf_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        pf_d2_9 = (.data$pf_d2_9 * n_vec_10[9] + .data$pf_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        pf_d3_9 = (.data$pf_d3_9 * n_vec_10[9] + .data$pf_d3_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        pf_d4_9 = (.data$pf_d4_9 * n_vec_10[9] + .data$pf_d4_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        pf_d5_9 = (.data$pf_d5_9 * n_vec_10[9] + .data$pf_d5_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        mo_d1_9 = (.data$mo_d1_9 * n_vec_10[9] + .data$mo_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        mo_d2_9 = (.data$mo_d2_9 * n_vec_10[9] + .data$mo_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        mo_d3_9 = (.data$mo_d3_9 * n_vec_10[9] + .data$mo_d3_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        mo_d4_9 = (.data$mo_d4_9 * n_vec_10[9] + .data$mo_d4_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        mo_d5_9 = (.data$mo_d5_9 * n_vec_10[9] + .data$mo_d5_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        az_d1_9 = (.data$az_d1_9 * n_vec_10[9] + .data$az_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        az_d2_9 = (.data$az_d2_9 * n_vec_10[9] + .data$az_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        ja_d1_9 = (.data$ja_d1_9 * n_vec_10[9] + .data$ja_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        ja_d2_9 = (.data$ja_d2_9 * n_vec_10[9] + .data$ja_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10])
      ) %>%
      select(
        .data$date, .data$pf_d1_1:.data$ja_d2_9, -.data$pf_d1_10, -.data$pf_d2_10, -.data$pf_d3_10,
        -.data$mo_d1_10, -.data$mo_d2_10, -.data$mo_d3_10, -.data$mo_d4_10, -.data$mo_d5_10,
        -.data$az_d1_10, -.data$az_d2_10, -.data$ja_d1_10, -.data$ja_d2_10
      )
  }
  
  # get daily proportion of people vaccinated ----------------------------------
  # take the difference for each row
  vac_schedule_daily <- data.frame(diff(as.matrix(vac_schedule %>% select(-date)))) %>%
    add_row(vac_schedule %>% select(-date) %>% slice(1), .before = 1) 
  
  # convert daily proportion to daily rate -------------------------------------
  # vac_schedule is the cumulative daily proportion 
  # cumulative sums need to be lagged by 1 day
  vac_schedule_rate <- vac_schedule_daily / (1 - lag(vac_schedule[,-1]))
  
  # add date variable back
  vac_schedule_daily <- vac_schedule_daily %>%
    mutate(date = date_vec) %>%
    select(date, everything()) # move date column to first position
  
  vac_schedule_rate <- vac_schedule_rate %>%
    mutate(date = date_vec) %>%
    select(date, everything()) # move date column to first position
  
  # add extra rows for dates further in the future -----------------------------
  if (add_extra_dates) {
    extra_dates <- seq.Date(from = as.Date(extra_start_date), 
                            to = as.Date(extra_end_date), by = 1)
    extra_dat <- data.frame(date = extra_dates) %>%
      full_join(vac_schedule_rate, extra_dates, by = "date") %>%
      mutate_at(vars(-.data$date), na_to_zero)
    vac_schedule_rate <- extra_dat
  } 
  
  # daily vaccination rate for each dose ---------------------------------------
  # first transform data frame to long format
  vac_rates_long <- vac_schedule_rate %>%
    pivot_longer(cols = !date, names_to = c("vac_product", "dose", "age_group"), 
                 names_sep = "_", values_to = "vac_rate") %>%
    mutate_at(vars(vac_rate), na_to_zero) 
  
  vac_rates_by_dose <- vac_rates_long %>%
    group_by(date, dose, age_group) %>%
    summarise(alpha = sum())
  
  # proportion vaccinated
  vac_prop_long <- vac_schedule_daily %>%
    pivot_longer(cols = !date, names_to = c("vac_product", "dose", "age_group"), 
                 names_sep = "_", values_to = "vac_prop") %>%
    mutate_at(vars(vac_prop), na_to_zero) 
  
  vac_prop_by_dose <- vac_prop_long %>%
    group_by(date, dose, age_group) %>%
    summarise(tot = sum(vac_prop))
  
  # get fraction of vaccines from each vaccine product administered on each day 
  vac_prop_long1 <- left_join(vac_prop_long, vac_prop_by_dose, by = c("date", "dose", "age_group")) %>%
    group_by(date, vac_product, dose, age_group) %>%
    mutate(frac = vac_prop/tot,
           frac = ifelse(is.nan(frac), 0, frac))

  # join with rate data frame
  vac_info_joined <- left_join(vac_rates_long, vac_prop_long1, by = c("date", "vac_product", "dose", "age_group")) 
  
  # get first day of vaccination with each dose
  # this will be used when calculating waning
  first_day_vac <- vac_info_joined %>%
    group_by(dose) %>%
    filter(vac_prop > 0, .preserve=TRUE) %>%
    summarise(first_day = first(date))
  
  # ----------------------------------------------------------------------------
  # Vaccine effectiveness
  # ----------------------------------------------------------------------------
  test <- vac_info_joined %>%
    filter(date >= as.Date("2021-01-04"),
           date <= as.Date("2021-03-04"),
           vac_product == "pf",
           dose == "d2",
           age_group %in% c(9)) 

  if (wane) {
    ve_dat <- left_join(test, first_day_vac, by = "dose") %>% # vac_info_joined %>%
      mutate(time_since_vac_start = ifelse(date >= first_day, date - first_day + 1, NA)) %>%
      group_by(vac_product, dose, age_group) %>%
      group_modify(~calc_waning(prop = .x$vac_prop, time_point = .x$time_since_vac_start))
  } else {
    waning <- c(rep(0, length(t_vec)))
  }
  

  
  # VE against infection -------------------------------------------------------
  # dose 1
  ve_p_dose1 <- calc_ve_w_waning(vac_rate = pf_dose1[, -1], ve_val = ve$pfizer[1], waning = waning)
  ve_m_dose1 <- calc_ve_w_waning(vac_rate = mo_dose1[, -1], ve_val = ve$moderna[1], waning = waning)
  ve_a_dose1 <- calc_ve_w_waning(vac_rate = az_dose1[, -1], ve_val = ve$astrazeneca[1], waning = waning)
  ve_j_dose1 <- calc_ve_w_waning(vac_rate = ja_dose1[, -1], ve_val = ve$jansen[1], waning = waning)
  
  # dose 2
  ve_p_dose2 <- calc_ve_w_waning(vac_rate = pf_dose2[, -1], ve_val = ve$pfizer[2], waning = waning)
  ve_m_dose2 <- calc_ve_w_waning(vac_rate = mo_dose2[, -1], ve_val = ve$moderna[2], waning = waning)
  ve_a_dose2 <- calc_ve_w_waning(vac_rate = az_dose2[, -1], ve_val = ve$astrazeneca[2], waning = waning)
  ve_j_dose2 <- calc_ve_w_waning(vac_rate = ja_dose2[, -1], ve_val = ve$jansen[1], waning = waning)
  
  # dose 3
  ve_p_dose3 <- calc_ve_w_waning(vac_rate = pf_dose3[, -1], ve_val = ve$pfizer[3], waning = waning)
  ve_m_dose3 <- calc_ve_w_waning(vac_rate = mo_dose3[, -1], ve_val = ve$moderna[3], waning = waning)
  
  # dose 4
  ve_p_dose4 <- calc_ve_w_waning(vac_rate = pf_dose4[, -1], ve_val = ve$pfizer[4], waning = waning)
  ve_m_dose4 <- calc_ve_w_waning(vac_rate = mo_dose4[, -1], ve_val = ve$moderna[4], waning = waning)
  
  # dose 3
  ve_p_dose5 <- calc_ve_w_waning(vac_rate = pf_dose5[, -1], ve_val = ve$pfizer[5], waning = waning)
  ve_m_dose5 <- calc_ve_w_waning(vac_rate = mo_dose5[, -1], ve_val = ve$moderna[5], waning = waning)
  
  # composite VE (against infection)
  comp_ve_dose1 <- frac_pf_dose1 * ve_p_dose1 +
    frac_mo_dose1 * ve_m_dose1 +
    frac_az_dose1 * ve_a_dose1 +
    frac_ja_dose1 * ve_j_dose1
  colnames(comp_ve_dose1) <- paste0("ve", name_suffix_d1)
  
  comp_ve_dose2 <- frac_pf_dose2 * ve_p_dose2 +
    frac_mo_dose2 * ve_m_dose2 +
    frac_az_dose2 * ve_a_dose2 +
    frac_ja_dose2 * ve_j_dose2
  colnames(comp_ve_dose2) <- paste0("ve", name_suffix_d2)
  
  comp_ve_dose3 <- frac_pf_dose3 * ve_p_dose3 + frac_mo_dose3 * ve_m_dose3 
  colnames(comp_ve_dose3) <- paste0("ve", name_suffix_d3)
  
  comp_ve_dose4 <- frac_pf_dose4 * ve_p_dose4 + frac_mo_dose4 * ve_m_dose4 
  colnames(comp_ve_dose4) <- paste0("ve", name_suffix_d4)
  
  comp_ve_dose5 <- frac_pf_dose5 * ve_p_dose5 + frac_mo_dose5 * ve_m_dose5 
  colnames(comp_ve_dose5) <- paste0("ve", name_suffix_d5)
  
  # eta
  eta_dose1 <- 1 - comp_ve_dose1
  eta_dose2 <- 1 - comp_ve_dose2
  eta_dose3 <- 1 - comp_ve_dose3
  eta_dose4 <- 1 - comp_ve_dose4
  eta_dose5 <- 1 - comp_ve_dose5
  
  # VE against hospitalisation
  # dose 1
  ve_hosp_p_dose1 <- calc_ve_w_waning(vac_rate = pf_dose1[, -1], ve_val = hosp_multiplier$pfizer[1], waning = waning)
  ve_hosp_m_dose1 <- calc_ve_w_waning(vac_rate = mo_dose1[, -1], ve_val = hosp_multiplier$moderna[1], waning = waning)
  ve_hosp_a_dose1 <- calc_ve_w_waning(vac_rate = az_dose1[, -1], ve_val = hosp_multiplier$astrazeneca[1], waning = waning)
  ve_hosp_j_dose1 <- calc_ve_w_waning(vac_rate = ja_dose1[, -1], ve_val = hosp_multiplier$jansen[1], waning = waning)
  
  # dose 2
  ve_hosp_p_dose2 <- calc_ve_w_waning(vac_rate = pf_dose2[, -1], ve_val = hosp_multiplier$pfizer[2], waning = waning)
  ve_hosp_m_dose2 <- calc_ve_w_waning(vac_rate = mo_dose2[, -1], ve_val = hosp_multiplier$moderna[2], waning = waning)
  ve_hosp_a_dose2 <- calc_ve_w_waning(vac_rate = az_dose2[, -1], ve_val = hosp_multiplier$astrazeneca[2], waning = waning)
  ve_hosp_j_dose2 <- calc_ve_w_waning(vac_rate = ja_dose2[, -1], ve_val = hosp_multiplier$jansen[2], waning = waning)
  
  # dose 3
  ve_hosp_p_dose3 <- calc_ve_w_waning(vac_rate = pf_dose3[, -1], ve_val = hosp_multiplier$pfizer[3], waning = waning)
  ve_hosp_m_dose3 <- calc_ve_w_waning(vac_rate = mo_dose3[, -1], ve_val = hosp_multiplier$moderna[3], waning = waning)
  
  # dose 4
  ve_hosp_p_dose4 <- calc_ve_w_waning(vac_rate = pf_dose4[, -1], ve_val = hosp_multiplier$pfizer[4], waning = waning)
  ve_hosp_m_dose4 <- calc_ve_w_waning(vac_rate = mo_dose4[, -1], ve_val = hosp_multiplier$moderna[4], waning = waning)
  
  # dose 5
  ve_hosp_p_dose5 <- calc_ve_w_waning(vac_rate = pf_dose5[, -1], ve_val = hosp_multiplier$pfizer[5], waning = waning)
  ve_hosp_m_dose5 <- calc_ve_w_waning(vac_rate = mo_dose5[, -1], ve_val = hosp_multiplier$moderna[5], waning = waning)
  
  # rate of hospitalisations multiplier
  hosp_mult_dose1 <- frac_pf_dose1 * ve_hosp_p_dose1 +
    frac_mo_dose1 * ve_hosp_m_dose1 +
    frac_az_dose1 * ve_hosp_a_dose1 +
    frac_ja_dose1 * ve_hosp_j_dose1
  colnames(hosp_mult_dose1) <- paste0("hosp_mult", name_suffix_d1)
  
  hosp_mult_dose2 <- frac_pf_dose2 * ve_hosp_p_dose2 +
    frac_mo_dose2 * ve_hosp_m_dose2 +
    frac_az_dose2 * ve_hosp_a_dose2 +
    frac_ja_dose2 * ve_hosp_j_dose2
  colnames(hosp_mult_dose2) <- paste0("hosp_mult", name_suffix_d2)
  
  hosp_mult_dose3 <- frac_pf_dose3 * ve_hosp_p_dose3 + frac_mo_dose3 * ve_hosp_m_dose3 
  colnames(hosp_mult_dose3) <- paste0("hosp_mult", name_suffix_d3)
  
  hosp_mult_dose4 <- frac_pf_dose4 * ve_hosp_p_dose4 + frac_mo_dose4 * ve_hosp_m_dose4 
  colnames(hosp_mult_dose4) <- paste0("hosp_mult", name_suffix_d4)
  
  hosp_mult_dose5 <- frac_pf_dose5 * ve_hosp_p_dose5 + frac_mo_dose5 * ve_hosp_m_dose5 
  colnames(hosp_mult_dose5) <- paste0("hosp_mult", name_suffix_d5)
  
  eta_hosp_dose1 <- hosp_mult_dose1
  eta_hosp_dose2 <- hosp_mult_dose2
  eta_hosp_dose3 <- hosp_mult_dose3
  eta_hosp_dose4 <- hosp_mult_dose4
  eta_hosp_dose5 <- hosp_mult_dose5
  
  # VE against hospitalisation
  # dose 1
  ve_trans_p_dose1 <- calc_ve_w_waning(vac_rate = pf_dose1[, -1], ve_val = ve_trans$pfizer[1], waning = waning)
  ve_trans_m_dose1 <- calc_ve_w_waning(vac_rate = mo_dose1[, -1], ve_val = ve_trans$moderna[1], waning = waning)
  ve_trans_a_dose1 <- calc_ve_w_waning(vac_rate = az_dose1[, -1], ve_val = ve_trans$astrazeneca[1], waning = waning)
  ve_trans_j_dose1 <- calc_ve_w_waning(vac_rate = ja_dose1[, -1], ve_val = ve_trans$jansen[1], waning = waning)
  
  # dose 2
  ve_trans_p_dose2 <- calc_ve_w_waning(vac_rate = pf_dose2[, -1], ve_val = ve_trans$pfizer[2], waning = waning)
  ve_trans_m_dose2 <- calc_ve_w_waning(vac_rate = mo_dose2[, -1], ve_val = ve_trans$moderna[2], waning = waning)
  ve_trans_a_dose2 <- calc_ve_w_waning(vac_rate = az_dose2[, -1], ve_val = ve_trans$astrazeneca[2], waning = waning)
  ve_trans_j_dose2 <- calc_ve_w_waning(vac_rate = ja_dose2[, -1], ve_val = ve_trans$jansen[1], waning = waning)
  
  # dose 3
  ve_trans_p_dose3 <- calc_ve_w_waning(vac_rate = pf_dose3[, -1], ve_val = ve_trans$pfizer[3], waning = waning)
  ve_trans_m_dose3 <- calc_ve_w_waning(vac_rate = mo_dose3[, -1], ve_val = ve_trans$moderna[3], waning = waning)
  
  # dose 4
  ve_trans_p_dose4 <- calc_ve_w_waning(vac_rate = pf_dose4[, -1], ve_val = ve_trans$pfizer[4], waning = waning)
  ve_trans_m_dose4 <- calc_ve_w_waning(vac_rate = mo_dose4[, -1], ve_val = ve_trans$moderna[4], waning = waning)
  
  # dose 5
  ve_trans_p_dose5 <- calc_ve_w_waning(vac_rate = pf_dose5[, -1], ve_val = ve_trans$pfizer[5], waning = waning)
  ve_trans_m_dose5 <- calc_ve_w_waning(vac_rate = mo_dose5[, -1], ve_val = ve_trans$moderna[5], waning = waning)
  
  # composite VE (against transmission)
  comp_ve_trans_dose1 <- frac_pf_dose1 * ve_trans_p_dose1 +
    frac_mo_dose1 * ve_trans_m_dose1 +
    frac_az_dose1 * ve_trans_a_dose1 +
    frac_ja_dose1 * ve_trans_j_dose1
  colnames(comp_ve_trans_dose1) <- paste0("ve_trans", name_suffix_d1)
  
  comp_ve_trans_dose2 <- frac_pf_dose2 * ve_trans_p_dose2 +
    frac_mo_dose2 * ve_trans_m_dose2 +
    frac_az_dose2 * ve_trans_a_dose2 +
    frac_ja_dose2 * ve_trans_j_dose2
  colnames(comp_ve_trans_dose2) <- paste0("ve_trans", name_suffix_d2)
  
  comp_ve_trans_dose3 <- frac_pf_dose3 * ve_trans_p_dose3 + frac_mo_dose3 * ve_trans_m_dose3
  colnames(comp_ve_trans_dose3) <- paste0("ve_trans", name_suffix_d3)
  
  comp_ve_trans_dose4 <- frac_pf_dose4 * ve_trans_p_dose4 + frac_mo_dose4 * ve_trans_m_dose4
  colnames(comp_ve_trans_dose4) <- paste0("ve_trans", name_suffix_d4)
  
  comp_ve_trans_dose5 <- frac_pf_dose5 * ve_trans_p_dose5 + frac_mo_dose5 * ve_trans_m_dose5
  colnames(comp_ve_trans_dose5) <- paste0("ve_trans", name_suffix_d5)
  
  # eta_trans
  eta_trans_dose1 <- 1 - comp_ve_trans_dose1
  eta_trans_dose2 <- 1 - comp_ve_trans_dose2
  eta_trans_dose3 <- 1 - comp_ve_trans_dose3
  eta_trans_dose4 <- 1 - comp_ve_trans_dose4
  eta_trans_dose5 <- 1 - comp_ve_trans_dose5
  
  # composite delay to protection
  delay_dose1 <- frac_pf_dose1 * delay$pfizer[1] +
    frac_mo_dose1 * delay$moderna[1] +
    frac_az_dose1 * delay$astrazeneca[1] +
    frac_ja_dose1 * delay$jansen[1]
  #delay_dose1 <- ifelse(delay_dose1 == 0, 1, delay_dose1) # this prevents from dividing by 0 in the ODEs
  colnames(delay_dose1) <- paste0("delay", name_suffix_d1)
  
  delay_dose2 <- frac_pf_dose2 * delay$pfizer[2] +
    frac_mo_dose2 * delay$moderna[2] +
    frac_az_dose2 * delay$astrazeneca[2] +
    frac_ja_dose2 * delay$jansen[2]
  #delay_dose2 <- ifelse(delay_dose2 == 0, 1, delay_dose2) # this prevents from dividing by 0 in the ODEs
  colnames(delay_dose2) <- paste0("delay", name_suffix_d2)
  
  delay_dose3 <- frac_pf_dose3 * delay$pfizer[3] + frac_mo_dose3 * delay$moderna[3]
  #delay_dose3 <- ifelse(delay_dose3 == 0, 1, delay_dose3) # this prevents from dividing by 0 in the ODEs
  colnames(delay_dose3) <- paste0("delay", name_suffix_d3)
  
  delay_dose4 <- frac_pf_dose4 * delay$pfizer[4] + frac_mo_dose4 * delay$moderna[4]
  #delay_dose4 <- ifelse(delay_dose4 == 0, 1, delay_dose4) # this prevents from dividing by 0 in the ODEs
  colnames(delay_dose4) <- paste0("delay", name_suffix_d4)
  
  delay_dose5 <- frac_pf_dose5 * delay$pfizer[5] + frac_mo_dose5 * delay$moderna[5]
  #delay_dose5 <- ifelse(delay_dose5 == 0, 1, delay_dose5) # this prevents from dividing by 0 in the ODEs
  colnames(delay_dose5) <- paste0("delay", name_suffix_d5)
  
  # output
  rtn <- list(
    alpha_dose1 = alpha_dose1,
    alpha_dose2 = alpha_dose2,
    alpha_dose3 = alpha_dose3,
    alpha_dose4 = alpha_dose4,
    alpha_dose5 = alpha_dose5,
    eta_dose1 = eta_dose1,
    eta_dose2 = eta_dose2,
    eta_dose3 = eta_dose3,
    eta_dose4 = eta_dose4,
    eta_dose5 = eta_dose5,
    delay_dose1 = delay_dose1,
    delay_dose2 = delay_dose2,
    delay_dose3 = delay_dose3,
    delay_dose4 = delay_dose4,
    delay_dose5 = delay_dose5,
    eta_hosp_dose1 = eta_hosp_dose1,
    eta_hosp_dose2 = eta_hosp_dose2,
    eta_hosp_dose3 = eta_hosp_dose3,
    eta_hosp_dose4 = eta_hosp_dose4,
    eta_hosp_dose5 = eta_hosp_dose5,
    eta_trans_dose1 = eta_trans_dose1,
    eta_trans_dose2 = eta_trans_dose2,
    eta_trans_dose3 = eta_trans_dose3,
    eta_trans_dose4 = eta_trans_dose4,
    eta_trans_dose5 = eta_trans_dose5
  )
  
  return(rtn)
}
