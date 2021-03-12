#' Calculate force of infection 
#' @param times vector of times
#' @param params list of parameter values
#' @return list of percentage of each group to be vaccinated at time point
#' and the composite ve for each group (based on how much of each vaccine is
#' used and the ve of each vaccine)
#' @keywords vacamole
#' @export
get_vac_rate <- function(times,params){
  with(as.list(c(times,params)),{
    #print(times)
  if (no_vac){
    total_dose1 <- c(rep(0,9))
    total_dose2 <- c(rep(0,9))
    eta_vec <- c(rep(1,9))
    eta2_vec <- c(rep(1,9))
    delay_dose1 <- 1
    delay_dose2 <- 1
  } else {
  # pfizer
  pf_dose1 <- vac_schedule %>%
    select(date, pf_d1_1:pf_d1_9)
  pf_dose1_cs <- cumsum(pf_dose1[,-1])
  pf_dose2 <- vac_schedule %>%
    select(date, pf_d2_1:pf_d2_9)
  pf_dose2_cs <- cumsum(pf_dose2[,-1])
  
  # moderna
  mo_dose1 <- vac_schedule %>%
    select(date, mo_d1_1:mo_d1_9)
  mo_dose1_cs <- cumsum(mo_dose1[,-1])
  mo_dose2 <- vac_schedule %>%
    select(date, mo_d2_1:mo_d2_9)
  mo_dose2_cs <- cumsum(mo_dose2[,-1])
  
  # astrazeneca
  az_dose1 <- vac_schedule %>%
    select(date, az_d1_1:az_d1_9)
  az_dose1_cs <- cumsum(az_dose1[,-1])
  az_dose2 <- vac_schedule %>%
    select(date, az_d2_1:az_d2_9)
  az_dose2_cs <- cumsum(az_dose2[,-1])
  
  # jansen
  ja_dose1 <- vac_schedule %>%
    select(date, Ja_d1_1:Ja_d1_9)
  ja_dose1_cs <- cumsum(ja_dose1[,-1])
  ja_dose2 <- vac_schedule %>%
    select(date, Ja_d2_1:Ja_d2_9)
  ja_dose2_cs <- cumsum(ja_dose2[,-1])
  
  # calculate composite VE
  time_point <- floor(times) + 1
  
  total_dose1 <- unlist(pf_dose1_cs[time_point,] + mo_dose1_cs[time_point,] + 
                          az_dose1_cs[time_point,] + ja_dose1_cs[time_point,])
  names(total_dose1) <- paste0("tot_",c(substr(names(total_dose1), 4, 7)))
  total_dose2 <- unlist(pf_dose2_cs[time_point,] + mo_dose2_cs[time_point,] + 
                          az_dose2_cs[time_point,] + ja_dose2_cs[time_point,])
  names(total_dose2) <- paste0("tot_",c(substr(names(total_dose2), 4, 7)))
  
  frac_pf_dose1 <- unlist(pf_dose1_cs[time_point,]/total_dose1)
  frac_pf_dose1 <- ifelse(is.nan(frac_pf_dose1), 0, frac_pf_dose1)
  frac_pf_dose2 <- unlist(pf_dose2_cs[time_point,]/total_dose2)
  frac_pf_dose2 <- ifelse(is.nan(frac_pf_dose2), 0, frac_pf_dose2)
  
  frac_mo_dose1 <- unlist(mo_dose1_cs[time_point,]/total_dose1)
  frac_mo_dose1 <- ifelse(is.nan(frac_mo_dose1), 0, frac_mo_dose1)
  frac_mo_dose2 <- unlist(mo_dose2_cs[time_point,]/total_dose2)
  frac_mo_dose2 <- ifelse(is.nan(frac_mo_dose2), 0, frac_mo_dose2)
  
  frac_az_dose1 <- unlist(az_dose1_cs[time_point,]/total_dose1)
  frac_az_dose1 <- ifelse(is.nan(frac_az_dose1), 0, frac_az_dose1)
  frac_az_dose2 <- unlist(az_dose2_cs[time_point,]/total_dose2)
  frac_az_dose2 <- ifelse(is.nan(frac_az_dose2), 0, frac_az_dose2)
  
  frac_ja_dose1 <- unlist(ja_dose1_cs[time_point,]/total_dose1)
  frac_ja_dose1 <- ifelse(is.nan(frac_ja_dose1), 0, frac_az_dose1)
  frac_ja_dose2 <- unlist(ja_dose2_cs[time_point,]/total_dose2)
  frac_ja_dose2 <- ifelse(is.nan(frac_ja_dose2), 0, frac_ja_dose2)

  comp_ve_dose1 <- frac_pf_dose1 * ve$pfizer[1] + 
    frac_mo_dose1 * ve$moderna[1] + 
    frac_az_dose1 * ve$astrazeneca[1] +
    frac_ja_dose1 * ve$jansen
  comp_ve_dose2 <- frac_pf_dose2 * ve$pfizer[2] + 
    frac_mo_dose2 * ve$moderna[2] + 
    frac_az_dose2 * ve$astrazeneca[2] +
    frac_ja_dose2 * ve$jansen
  
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
  }
  list(c(total_dose1,total_dose2,eta_vec,eta2_vec,
         delay_dose1, delay_dose2))
  
  }) 
}
