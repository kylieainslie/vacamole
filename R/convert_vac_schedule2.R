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
#' @param k_inf logistic growth rate for waning curve against infection
#' @param k_sev logistic growth rate for waning curve against severe disease
#' @param t0 the time point at the midpoint of the logistic curve (where 50\%
#'  waning occurs)
#' @return data frame of vaccination rate by day, dose, vaccine product, and age
#' group and weighted VE and delay to protection by day, dose, age group, and 
#' outcome
#' @keywords vacamole
#' @import tidyr
#' @import dplyr
#' @export
convert_vac_schedule2 <- function(vac_schedule,
                                  ve_pars,
                                  wane = FALSE,
                                  k_inf = 0.012,
                                  k_sev = 0.006,
                                  t0 = 365
                                  ){
  
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
  
  # daily vaccination rate for each dose ---------------------------------------
  # first transform data frame to long format
  vac_rates_long <- vac_schedule_rate %>%
    pivot_longer(cols = !date, names_to = c("vac_product", "dose", "age_group"), 
                 names_sep = "_", values_to = "vac_rate") %>%
    mutate_at(vars(vac_rate), na_to_zero) 
  
  vac_rates_by_dose <- vac_rates_long %>%
    group_by(date, dose, age_group) %>%
    summarise(alpha = sum())
  
  # cumulative vaccination proportion for each dose/product --------------------
  # convert vac_schedule to long format
  vac_schedule_cs <- vac_schedule %>%
    pivot_longer(cols = !date, names_to = c("vac_product", "dose", "age_group"), 
                 names_sep = "_", values_to = "vac_prop") %>%
    mutate_at(vars(vac_prop), na_to_zero)
    
  # proportion vaccinated up to current day by dose and age group (summed over
  # vac_products)
  vac_prop_total <- vac_schedule_cs %>%
    group_by(date, dose, age_group) %>%
    summarise(total = sum(vac_prop)) %>% 
    ungroup() 
  
  vac_prop <- left_join(vac_schedule_cs, vac_prop_total, by = c("date", "dose", "age_group")) %>%
    group_by(date, vac_product, dose, age_group) %>%
    mutate(frac = vac_prop/total,
           frac = ifelse(is.nan(frac), 0, frac)) 
  

  # join with rate data frame
  vac_info_joined <- left_join(vac_rates_long, vac_prop, by = c("date", "vac_product", "dose", "age_group")) %>%
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
    # calculate waning against infection
    ve_dat1 <- left_join(vac_info_joined, first_day_vac, by = "dose") %>% 
      mutate(time_since_vac_start = ifelse(date >= first_day, date - first_day + 1, NA)) %>%
      group_by(vac_product, dose, age_group) %>%
      group_modify(~calc_waning(prop = .x$vac_prop, time_point = .x$time_since_vac_start,
                   k = k_inf, t0 = t0)) %>%
      mutate(date = row_number() + vac_info_joined$date[1] - 1) %>% # convert t to date
      rename(w_inf = w) %>%                                         # rename waning column
      ungroup()
    
    # calculate waning against severe disease
    ve_dat2 <- left_join(vac_info_joined, first_day_vac, by = "dose") %>% 
      mutate(time_since_vac_start = ifelse(date >= first_day, date - first_day + 1, NA)) %>%
      group_by(vac_product, dose, age_group) %>%
      group_modify(~calc_waning(prop = .x$vac_prop, time_point = .x$time_since_vac_start,
                   k = k_sev, t0 = t0)) %>%
      mutate(date = row_number() + vac_info_joined$date[1] - 1) %>% # convert t to date
      rename(w_sev = w) %>%                                         # rename waning column
      ungroup()
      
      # calculate VE with waning
    ve_dat <- left_join(ve_dat1, ve_dat2, by = c("vac_product", "dose", "age_group", "t", "date")) %>%
      left_join(., ve_pars, by = c("vac_product", "dose", "age_group")) %>%
      left_join(., vac_info_joined, by = c("date","vac_product", "dose", "age_group")) %>%
      group_by(dose, age_group, outcome) %>%
      mutate(ve_wane_inf = ifelse(ve -  w_inf < 0, 0, ve - w_inf),
             ve_wane_sev = ifelse(ve -  w_sev < 0, 0, ve - w_sev))
    
    # calculate composite VE ---------------------------------------------------
    ve_comp <- ve_dat %>%
      group_by(date, dose, age_group, outcome) %>%
      summarise(alpha = sum(vac_rate, na.rm = TRUE),
                comp_ve_inf = sum(frac * ve_wane_inf, na.rm = TRUE),
                comp_ve_sev = sum(frac * ve_wane_sev, na.rm = TRUE),
                comp_delay = sum(frac * delay, na.rm = TRUE)) 
      
    # fill in NA values --------------------------------------------------------
    # even when no people are vaccinated on a given date, 
    # comp_ve should be > 0
    ve_comp_fill <- ve_comp %>%
      mutate(comp_ve_inf = ifelse(comp_ve_inf == 0, NA, comp_ve_inf),
             comp_ve_sev = ifelse(comp_ve_sev == 0, NA, comp_ve_sev),
             comp_delay = ifelse(comp_delay == 0, NA, comp_delay)) %>% 
      tidyr::fill(comp_ve_inf:comp_delay, .direction = c("down"))
      # tidyr::fill(comp_ve_sev, .direction = c("down")) %>%
      # tidyr::fill(comp_delay, .direction = c("down")) 
      
    # get etas -----------------------------------------------------------------
    ve_comp_eta <- ve_comp_fill %>%
      mutate(eta_inf = 1 - comp_ve_inf,
             eta_inf = ifelse(is.na(eta_inf), 1, eta_inf),
             eta_sev = 1 - comp_ve_sev,
             eta_sev = ifelse(is.na(eta_sev), 1, eta_sev),
             comp_delay = ifelse(is.na(comp_delay), 1, comp_delay)) %>%
      pivot_longer(., alpha:eta_sev, names_to = "param", values_to = "value")
    
    # output -------------------------------------------------------------------
    # subset VE comp so eta value matches outcome 
    # (ex: eta_inf only for outcome = infection)
    rtn <- ve_comp_eta %>%
      filter(!(outcome == "hospitalisation" & param == "eta_inf"),
             !(outcome == "hospitalisation" & param == "comp_ve_inf"),
             !(outcome == "infection" & param == "eta_sev"),
             !(outcome == "infection" & param == "comp_ve_sev"),
             !(outcome == "transmission" & param == "eta_sev"),
             !(outcome == "transmission" & param == "comp_ve_sev")
             ) %>%
      mutate(param = ifelse((param == "comp_ve_inf" | param == "comp_ve_sev"), "comp_ve", param),
             param = ifelse((param == "eta_inf" | param == "eta_sev"), "eta", param),
      ) %>%
      ungroup()
    

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
    rtn <- ve_comp
  }

  # output ---------------------------------------------------------------------
  return(rtn)
}
