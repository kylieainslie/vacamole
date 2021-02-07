#' Calculate force of infection 
#' @param vac_schedule data frame of vac schedule for all three vaccines and 
#' doses 1 and 2 by age group
#' @param ve list of vaccine efficacy estimates for doses 1 and 2
#' @param time time point
#' @return list of percentage of each group to be vaccinated at time point
#' and the composite ve for each group (based on how much of each vaccine is
#' used and the ve of each vaccine)
#' @keywords vacamole
#' @export
get_vac_rate <- function(times, params){
  with(as.list(c(times,params)),{
    #print(times)
  # pfizer
  pf_dose1 <- vac_schedule %>%
    select(date, pf_d1_1:pf_d1_9)
  pf_dose2 <- vac_schedule %>%
    select(date, pf_d2_1:pf_d2_9)
  
  # moderna
  mo_dose1 <- vac_schedule %>%
    select(date, mo_d1_1:mo_d1_9)
  mo_dose2 <- vac_schedule %>%
    select(date, mo_d2_1:mo_d2_9)
  
  # astrazeneca
  az_dose1 <- vac_schedule %>%
    select(date, az_d1_1:az_d1_9)
  az_dose2 <- vac_schedule %>%
    select(date, az_d2_1:az_d2_9)
  
  # calculate composite VE
  time_point <- ifelse(floor(times) == 0, floor(times) + 1, floor(times))
  
  total_dose1 <- unlist(pf_dose1[time_point,-1] + mo_dose1[time_point,-1] + az_dose1[time_point,-1])
  names(total_dose1) <- paste0("tot_",c(substr(names(total_dose1), 4, 7)))
  total_dose2 <- unlist(pf_dose2[time_point,-1] + mo_dose2[time_point,-1] + az_dose2[time_point,-1])
  names(total_dose2) <- paste0("tot_",c(substr(names(total_dose2), 4, 7)))
  
  comp_ve_dose1 <- unlist((pf_dose1[time_point,-1]/total_dose1) * ve$pfizer[1] + 
                   (mo_dose1[time_point,-1]/total_dose1) * ve$moderna[1] + 
                   (az_dose1[time_point,-1]/total_dose1) * ve$astrazeneca[1])
  
  comp_ve_dose2 <- unlist((pf_dose2[time_point,-1]/total_dose2) * ve$pfizer[2] + 
                   (mo_dose2[time_point,-1]/total_dose2) * ve$moderna[2] + 
                   (az_dose2[time_point,-1]/total_dose2) * ve$astrazeneca[2])
  
  eta_vec <- 1 - ifelse(is.nan(comp_ve_dose1), 0, comp_ve_dose1)
  names(eta_vec) <- paste0("eta_",1:9)
  eta2_vec <- 1- ifelse(is.nan(comp_ve_dose2), 0, comp_ve_dose2)
  names(eta2_vec) <- paste0("eta2_",1:9)
  
  list(c(total_dose1,total_dose2,eta_vec,eta2_vec))
  
  }) 
}
