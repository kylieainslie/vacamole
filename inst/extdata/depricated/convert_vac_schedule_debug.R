#' convert cumulative vaccination schedule to non-cumulative -------------------
#' @param vac_schedule a data frame that the proportion of the population who 
#' receives vaccines
#' at each time point. Rows are time points, columns are 
#' <vaccine type>_<dose>_<age group>. For example,
#' the first column is the proportion of individuals in age group 1 who receive 
#' dose 1 of the first vaccine
#' type. The function assumes 9 or 10 age groups (1, .., 10) and four vaccine 
#' types (pf, mo, az, ja), each with a 2-dose regimen (d1, d2).
#' @param ve_pars a data frame with VE estimates by vaccine product, dose, outcome, 
#' and age group. The data frame should have the following columns: vac_product,	
#' dose,	age_group, outcome,	variant, delay,	ve
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
#' @return data frame of vaccination rate by day, dose, vaccine product, and age
#' group and weighted VE and delay to protection by day, dose, age group, and 
#' outcome
#' @keywords vacamole
#' @import tidyr
#' @import dplyr
#' @export
convert_vac_schedule_debug <- function(vac_schedule,
                                  ve_pars,
                                  wane = FALSE,
                                  add_extra_dates = FALSE,
                                  extra_start_date,
                                  extra_end_date){
  
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
        .data$date, .data$pf_d1_1:.data$ja_d2_9, 
        -.data$pf_d1_10, -.data$pf_d2_10, -.data$pf_d3_10, -.data$pf_d4_10, -.data$pf_d5_10,
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
  vac_info_joined <- left_join(vac_rates_long, vac_prop_long1, by = c("date", "vac_product", "dose", "age_group")) %>%
    mutate(age_group = as.numeric(age_group)) 
  
  # get first day of vaccination with each dose
  # this will be used when calculating waning
  first_day_vac <- vac_info_joined %>%
    group_by(dose) %>%
    filter(vac_prop > 0, .preserve=TRUE) %>%
    summarise(first_day = first(date))
  
  # ----------------------------------------------------------------------------
  # Vaccine effectiveness
  # ----------------------------------------------------------------------------
  if (wane) {
    ve_dat <- left_join(vac_info_joined, first_day_vac, by = "dose") %>% 
      mutate(time_since_vac_start = ifelse(date >= first_day, date - first_day + 1, NA)) %>%
      group_by(vac_product, dose, age_group) %>%
      # calculate waning
      group_modify(~calc_waning(prop = .x$vac_prop, time_point = .x$time_since_vac_start)) %>%
      # convert t to date
      mutate(date = row_number() + vac_info_joined$date[1] - 1) %>%
      ungroup() %>%
      # calculate VE with waning
      left_join(., ve_pars, by = c("vac_product", "dose", "age_group")) %>%
      left_join(., vac_info_joined, by = c("date","vac_product", "dose", "age_group")) %>%
      group_by(dose, age_group, outcome) %>%
      mutate(ve_wane = ve - (ve * w))
    
    # calculate composite VE ---------------------------------------------------
    ve_comp <- ve_dat %>%
      group_by(date, dose, age_group, outcome) %>%
      summarise(alpha = sum(vac_rate),
                comp_ve = sum(frac * ve_wane),
                comp_delay = sum(frac * delay)) %>%
      ungroup() %>%      
      # fill in zeros with previous value
      # even when no people are vaccinated on a given date, comp_ve should be > 0
      mutate(comp_ve = ifelse(comp_ve == 0, NA, comp_ve),
             comp_delay = ifelse(comp_ve == 0, NA, comp_delay)) %>% 
      tidyr::fill(comp_ve, .direction = c("down")) %>%
      tidyr::fill(comp_delay, .direction = c("down")) %>%
      mutate(eta = 1 - comp_ve,
             eta = ifelse(is.na(eta), 1, eta),
             comp_delay = ifelse(is.na(comp_delay), 1, comp_delay)) %>%
      pivot_longer(., alpha:eta, names_to = "param", values_to = "value")

  } else {
    ve_dat <- left_join(vac_info_joined, first_day_vac, by = "dose") %>% # vac_info_joined %>%
      mutate(time_since_vac_start = ifelse(date >= first_day, date - first_day + 1, NA)) %>%
      left_join(., ve_pars, by = c("vac_product", "dose", "age_group")) 
    
    # calculate composite VE ---------------------------------------------------
    ve_comp <- ve_dat %>%
      group_by(date, dose, age_group, outcome) %>%
      summarise(alpha = sum(vac_rate),
                comp_ve = sum(frac * ve),
                comp_delay = sum(frac * delay)) %>%
      ungroup() %>%
      # fill in zeros with previous value
      # even when no people are vaccinated on a given date, comp_ve should be > 0
      mutate(comp_ve = ifelse(comp_ve == 0, NA, comp_ve),
             comp_delay = ifelse(comp_ve == 0, NA, comp_delay)) %>% 
      tidyr::fill(comp_ve, .direction = c("down")) %>%
      tidyr::fill(comp_delay, .direction = c("down")) %>%
      mutate(eta = 1 - comp_ve,
             eta = ifelse(is.na(eta), 1, eta),
             comp_delay = ifelse(is.na(comp_delay), 1, comp_delay)) %>%
      pivot_longer(., alpha:eta, names_to = "param", values_to = "value")
    
  }

  # output ---------------------------------------------------------------------
  return(ve_comp)
}
