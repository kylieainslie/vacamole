#' convert cumulative vaccination schedule to non-cumulative ----------------------------------
#' @param vac_schedule a data frame that the proportion of the population who receives vaccines
#' at each time point. Rows are time points, columns are <vaccine type>_<dose>_<age group>. For example,
#' the first column is the proportion of individuals in age group 1 who receive dose 1 of the first vaccine
#' type. The function assumes 9 or 10 age groups (1, .., 10) and four vaccine types (pf, mo, az, ja), each with
#' a 2-dose regimen (d1, d2).
#' @param ve a named list of vaccine effectiveness against infection for each dose of each vaccine type.
#' @param hosp_multiplier a named list of the values for the vaccine effectiveness against hospitalization
#' for each dose of each vaccine type. Vaccine effectiveness against hospitalization is incorporated as a
#' multiplier on the probability of being hospitalized after infection as (1 – VE_hospitalization) divided by (1-VE_infection)
#' to account for the the inclusion of people who are never infected (and thus never hospitalized) included
#' in the estimation of VE against hospitalization.
#' @param delay a named list of the time to protection for each dose of each vaccine type.
#' @param ve_trans a named list of vaccine effectiveness against transmission for each dose of each vaccine type.
#' @param add_child_vac logical, if TRUE 5-11 year olds are vaccinated
#' @param child_vac_coverage total vaccine coverage in 5-11 to achieve (between 0 and 1) if add_child_vac = TRUE
#' @param child_doses_per_day the number of doses to be administered to 5-11 year olds if add_child_vac = TRUE
#' @param child_vac_start_date character string of the date (YYYY-MM-DD format) to start vaccinating 5-11 year
#' olds, if add_child_vac = TRUE
#' @param wane logical, if TRUE vaccine effectiveness wanes by a logistic function parameterized by arguments
#' k and t0.
#' @param k logistic growth rate
#' @param t0 the time point at the midpoint of the logistic curve (where 50\% waning occurs)
#' @param add_extra_dates logical, if TRUE add extra rows to the vaccination schedule until extra_end_date
#' @param extra_end_date character string of the date (YYYY-MM-DD format) to end the vaccination schedule (which will
#' also be the last date of the simulation)
#' @return list of vaccination rate by day and dose and weighted VE and delay to protection by day and dose
#' @keywords vacamole
#' @import tidyr
#' @import dplyr
#' @export
convert_vac_schedule <- function(vac_schedule,
                                 ve,
                                 hosp_multiplier,
                                 delay,
                                 ve_trans,
                                 add_child_vac = FALSE,
                                 child_vac_coverage = 0.75,
                                 child_doses_per_day = 50000,
                                 child_vac_start_date = "2021-09-01",
                                 wane = FALSE,
                                 k = 0.03,
                                 t0 = 180,
                                 add_extra_dates = FALSE,
                                 extra_start_date = "2022-01-01",
                                 extra_end_date = "2022-03-31"){
  
# check if there are 9 age groups or 10 age groups --------------------------------------------
  # extract age group part of each relevant column name of vac_schedule
  names_check <- str_extract(names(vac_schedule), pattern = "(?<=[a-z]{2}_d[0-9]_)[0-9]+")
  # find the largest age group
  num_age_groups <- names_check %>%
    as.numeric() %>%
    max(na.rm = TRUE) %>%
    as.integer()
  
# remove any redundant columns 
  vac_schedule <- vac_schedule %>% select(date | matches("[a-z]{2}_d[0-9]_[0-9]+"))
  
# some demographic info -----------------------------------------------------------------------
  age_dist_10 <- c(
    0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463,
    0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671
  )
  n <- 17407585 # Dutch population size
  n_vec_10 <- n * age_dist_10

  date_vec <- as.Date(vac_schedule$date, format = "%m/%d/%Y")

  # if 10 age groups then: ----------------------------------------------------------------------
  if (num_age_groups == 10) {

    # take the difference for each row ------------------------------------------------------------
    vac_schedule_orig <- data.frame(diff(as.matrix(vac_schedule %>% select(-date)))) %>%
      add_row(vac_schedule %>% select(-date) %>% slice(1), .before = 1) %>%
      mutate(
        date = date_vec,
        pf_d1_9 = (.data$pf_d1_9 * n_vec_10[9] + .data$pf_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        pf_d2_9 = (.data$pf_d2_9 * n_vec_10[9] + .data$pf_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        mo_d1_9 = (.data$mo_d1_9 * n_vec_10[9] + .data$mo_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        mo_d2_9 = (.data$mo_d2_9 * n_vec_10[9] + .data$mo_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        az_d1_9 = (.data$az_d1_9 * n_vec_10[9] + .data$az_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        az_d2_9 = (.data$az_d2_9 * n_vec_10[9] + .data$az_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        ja_d1_9 = (.data$ja_d1_9 * n_vec_10[9] + .data$ja_d1_10 * n_vec_10[10]) / sum(n_vec_10[9:10]),
        ja_d2_9 = (.data$ja_d2_9 * n_vec_10[9] + .data$ja_d2_10 * n_vec_10[10]) / sum(n_vec_10[9:10])
      ) %>%
      select(
        .data$date, .data$pf_d1_1:.data$ja_d2_9, -.data$pf_d1_10, -.data$pf_d2_10, -.data$mo_d1_10,
        -.data$mo_d2_10, -.data$az_d1_10, -.data$az_d2_10, -.data$ja_d1_10, -.data$ja_d2_10
      )
  } else if (num_age_groups == 9) { # if 9 age groups
    # keep only columns that are of the form <vaccine type>_<dose>_<age group> or the date column
    vac_schedule_orig <- data.frame(diff(as.matrix(vac_schedule %>% select(-date)))) %>%
      add_row(vac_schedule %>% select(-date) %>% slice(1), .before = 1) %>%
      mutate(date = date_vec) %>%
      select(date, everything()) # move date column to first position
  } else { # if not 9 or 10 age groups, can't convert vac_schedule
    stop("number of age groups in vac_schedule is not 9 or 10, can't convert")
  }

  vac_schedule_orig_new <- vac_schedule_orig

  # add extra rows for dates further in the future (so there's no error when running the model)
  if (add_extra_dates) {
    extra_dates <- seq.Date(from = as.Date(extra_start_date), to = as.Date(extra_end_date), by = 1)
    na_to_zero <- function(x) {
      ifelse(is.na(x), 0, x)
    }
    extra_dat <- data.frame(date = extra_dates) %>%
      full_join(vac_schedule_orig_new, extra_dates, by = "date") %>%
      mutate_at(vars(-.data$date), na_to_zero)
    vac_schedule_orig_new <- extra_dat
  }

  # add vaccination of children 5-11
  # we assume the second dose is given 3 weeks after the first dose
  if (add_child_vac) {

    # 5-9 year olds
    p_child_5_9_doses <- 0.5 * child_vac_coverage
    n_child_5_9_doses <- n_vec_10[1] * p_child_5_9_doses
    days_child_5_9 <- ceiling(n_child_5_9_doses / child_doses_per_day)
    p_child_5_9_doses_per_day <- p_child_5_9_doses / days_child_5_9

    # 10-11 year olds
    p_child_10_11_doses <- 0.2 * child_vac_coverage
    n_child_10_11_doses <- n_vec_10[2] * p_child_10_11_doses
    days_child_10_11 <- ceiling(n_child_10_11_doses / child_doses_per_day)
    p_child_10_11_doses_per_day <- p_child_10_11_doses / days_child_10_11

    # define start and end dates for vaccinating each group
    start_date_5_9 <- as.Date(child_vac_start_date) + days_child_10_11
    end_date_5_9 <- as.Date(child_vac_start_date) + days_child_10_11 + days_child_5_9
    start_date_10_11 <- as.Date(child_vac_start_date)
    end_date_10_11 <- as.Date(child_vac_start_date) + days_child_10_11

    vac_schedule_orig_new <- vac_schedule_orig_new %>%
      # vaccinate 10-11 first, then 5-9
      mutate(
        pf_d1_1 = ifelse(date >= start_date_5_9 & date <= end_date_5_9, .data$pf_d1_1 + p_child_5_9_doses_per_day, .data$pf_d1_1),
        pf_d2_1 = ifelse(date >= start_date_5_9 + 21 & date <= end_date_5_9 + 21, .data$pf_d2_1 + p_child_5_9_doses_per_day, .data$pf_d2_1),
        pf_d1_2 = ifelse(date >= start_date_10_11 & date <= end_date_10_11, .data$pf_d1_2 + p_child_10_11_doses_per_day, .data$pf_d1_2),
        pf_d2_2 = ifelse(date >= start_date_10_11 + 21 & date <= end_date_10_11 + 21, .data$pf_d2_2 + p_child_10_11_doses_per_day, .data$pf_d2_2)
      )

    vac_schedule_new_cs <- cumsum(vac_schedule_orig_new[, -1])
    # vac_out <- data.frame(date = vac_schedule_orig_new$date, vac_schedule_new_cs)
    # write.csv(vac_out, file = "inst/extdata/data/vaccination_scenarios/Cum_upt20210701 Basis 75% in 5+ KA.csv")
  } else {
    vac_schedule_new <- vac_schedule_orig_new
    vac_schedule_new_cs <- cumsum(vac_schedule_orig_new[, -1])
  }

  # create separate data frame for each vaccine -------------------------------------------------
  # pfizer
  pf_dose1 <- vac_schedule_orig_new %>% select(date, .data$pf_d1_1:.data$pf_d1_9)
  pf_dose1_cs <- vac_schedule_new_cs %>% select(.data$pf_d1_1:.data$pf_d1_9)
  pf_dose2 <- vac_schedule_orig_new %>% select(date, .data$pf_d2_1:.data$pf_d2_9)
  pf_dose2_cs <- vac_schedule_new_cs %>% select(.data$pf_d2_1:.data$pf_d2_9)

  # moderna
  mo_dose1 <- vac_schedule_orig_new %>% select(date, .data$mo_d1_1:.data$mo_d1_9)
  mo_dose1_cs <- vac_schedule_new_cs %>% select(.data$mo_d1_1:.data$mo_d1_9)
  mo_dose2 <- vac_schedule_orig_new %>% select(date, .data$mo_d2_1:.data$mo_d2_9)
  mo_dose2_cs <- vac_schedule_new_cs %>% select(.data$mo_d2_1:.data$mo_d2_9)

  # astrazeneca
  az_dose1 <- vac_schedule_orig_new %>% select(date, .data$az_d1_1:.data$az_d1_9)
  az_dose1_cs <- vac_schedule_new_cs %>% select(.data$az_d1_1:.data$az_d1_9)
  az_dose2 <- vac_schedule_orig_new %>% select(date, .data$az_d2_1:.data$az_d2_9)
  az_dose2_cs <- vac_schedule_new_cs %>% select(.data$az_d2_1:.data$az_d2_9)

  # jansen
  ja_dose1 <- vac_schedule_orig_new %>% select(date, .data$ja_d1_1:.data$ja_d1_9)
  ja_dose1_cs <- vac_schedule_new_cs %>% select(.data$ja_d1_1:.data$ja_d1_9)
  ja_dose2 <- vac_schedule_orig_new %>% select(date, .data$ja_d2_1:.data$ja_d2_9)
  ja_dose2_cs <- vac_schedule_new_cs %>% select(.data$ja_d2_1:.data$ja_d2_9)

  # calculate composite VE ---------------------------------------------------------------------
  name_suffix_d1 <- c(substr(names(pf_dose1[, -1]), 3, 7))
  name_suffix_d2 <- c(substr(names(pf_dose2[, -1]), 3, 7))

  # daily vaccination rate
  alpha_dose1 <- pf_dose1[, -1] + mo_dose1[, -1] + az_dose1[, -1] + ja_dose1[, -1]
  names(alpha_dose1) <- paste0("alpha", name_suffix_d1)
  alpha_dose2 <- pf_dose2[, -1] + mo_dose2[, -1] + az_dose2[, -1] + ja_dose2[, -1]
  names(alpha_dose2) <- paste0("alpha", name_suffix_d2)

  # cumulative vaccination percentage
  total_dose1 <- pf_dose1_cs + mo_dose1_cs + az_dose1_cs + ja_dose1_cs
  names(total_dose1) <- paste0("tot", name_suffix_d1)
  total_dose2 <- pf_dose2_cs + mo_dose2_cs + az_dose2_cs + ja_dose2_cs
  names(total_dose2) <- paste0("tot", name_suffix_d2)

  # fraction of each vaccine given up to current time point
  frac_pf_dose1 <- data.matrix(pf_dose1_cs / total_dose1)
  frac_pf_dose1 <- ifelse(is.nan(frac_pf_dose1), 0, frac_pf_dose1)
  frac_pf_dose2 <- data.matrix(pf_dose2_cs / total_dose2)
  frac_pf_dose2 <- ifelse(is.nan(frac_pf_dose2), 0, frac_pf_dose2)

  frac_mo_dose1 <- data.matrix(mo_dose1_cs / total_dose1)
  frac_mo_dose1 <- ifelse(is.nan(frac_mo_dose1), 0, frac_mo_dose1)
  frac_mo_dose2 <- data.matrix(mo_dose2_cs / total_dose2)
  frac_mo_dose2 <- ifelse(is.nan(frac_mo_dose2), 0, frac_mo_dose2)

  frac_az_dose1 <- data.matrix(az_dose1_cs / total_dose1)
  frac_az_dose1 <- ifelse(is.nan(frac_az_dose1), 0, frac_az_dose1)
  frac_az_dose2 <- data.matrix(az_dose2_cs / total_dose2)
  frac_az_dose2 <- ifelse(is.nan(frac_az_dose2), 0, frac_az_dose2)

  frac_ja_dose1 <- data.matrix(ja_dose1_cs / total_dose1)
  frac_ja_dose1 <- ifelse(is.nan(frac_ja_dose1), 0, frac_ja_dose1)
  frac_ja_dose2 <- data.matrix(ja_dose2_cs / total_dose2)
  frac_ja_dose2 <- ifelse(is.nan(frac_ja_dose2), 0, frac_ja_dose2)

  # calculate amount of waning
  # it is based on the weighted overage of the VE, daily vaccination rate, and
  # time since vaccination
  t_vec <- seq(1, dim(pf_dose1)[1], by = 1)

  if (wane) {
    waning <- (1 / (1 + exp(-k * (t_vec - t0))))
  } else {
    waning <- c(rep(0, length(t_vec)))
  }

  # VE against infection
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

  # eta
  eta_dose1 <- 1 - comp_ve_dose1
  eta_dose2 <- 1 - comp_ve_dose2

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
  ve_hosp_j_dose2 <- calc_ve_w_waning(vac_rate = ja_dose2[, -1], ve_val = hosp_multiplier$jansen[1], waning = waning)

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
  colnames(hosp_mult_dose2) <- paste0("hosp_mult", name_suffix_d1)

  eta_hosp_dose1 <- hosp_mult_dose1
  eta_hosp_dose2 <- hosp_mult_dose2

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

  # eta_trans
  eta_trans_dose1 <- 1 - comp_ve_trans_dose1
  eta_trans_dose2 <- 1 - comp_ve_trans_dose2


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

  rtn <- list(
    alpha_dose1 = alpha_dose1,
    alpha_dose2 = alpha_dose2,
    eta_dose1 = eta_dose1,
    eta_dose2 = eta_dose2,
    delay_dose1 = delay_dose1,
    delay_dose2 = delay_dose2,
    eta_hosp_dose1 = eta_hosp_dose1,
    eta_hosp_dose2 = eta_hosp_dose2,
    eta_trans_dose1 = eta_trans_dose1,
    eta_trans_dose2 = eta_trans_dose2
  )

  return(rtn)
}
