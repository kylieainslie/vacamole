#' Calculate force of infection 
#' @param dat data frame of states
#' @param params list of parameter values
#' @param contact_matrix contact matrix
#' @return matrix of force of infection in each age group (columns) at each
#' time point (rows)
#' @keywords vacamole
#' @export
get_foi <- function(dat, params, vac_inputs){

    
    # sum over different I states for each time step and age group
    E_all <- as.matrix(dat$E + dat$Ev_1d + dat$Ev_2d)
    I_all <- as.matrix(dat$I + (vac_inputs$eta_trans_dose1 * dat$Iv_1d) + (vac_inputs$eta_trans_dose2 * dat$Iv_2d))
    H_all <- as.matrix(dat$H + dat$Hv_1d + dat$Hv_2d)
    
    
    # calculate force of infection for each time point
    time_vec <- 1:dim(E_all)[1]
    ic_admin <- rep(0, length(time_vec))
    cases    <- rep(0, length(time_vec))
    criteria <- rep(0, length(time_vec))
    cm_check <- rep(NA, length(time_vec))
    flag_relaxed <- 0
    flag_very_relaxed <- 0
    flag_normal <- 0
    
    for(t in time_vec){
      ic_admin[t] <- sum(params$i1 * H_all[t,])
      cases[t] <- sum(params$sigma * E_all[t,] * params$p_report)
      
      criteria[t] <- (params$use_cases) * cases[t] + (!params$use_cases) * ic_admin[t] 
      
      tmp2 <- choose_contact_matrix(params, t, criteria[t], flag_relaxed, 
                                    flag_very_relaxed, flag_normal)
      contact_mat <- tmp2$contact_matrix
      flag_relaxed <- tmp2$flag_relaxed 
      flag_very_relaxed <- tmp2$flag_very_relaxed 
      flag_normal <- tmp2$flag_normal
      
      # determine force of infection ----------------------------------
      calendar_day <- params$t_calendar_start + t
      
      # lambda <- params$beta * (contact_mat %*% (I_all[t,]))
      beta_t <- params$beta * (1 + params$beta1 * cos(2 * pi * calendar_day/365.24))
      lambda <- beta_t * (contact_mat %*% I_all[t,])
      
      # ---------------------------------------------------------------
      
      # check for which contact metrix was selected -------------------
      cm_check[t] <- ifelse(identical(contact_mat, params$c_lockdown), "c_lockdown", 
                            ifelse(identical(contact_mat, params$c_relaxed),"c_relaxed",
                                   ifelse(identical(contact_mat, params$c_very_relaxed),"c_very_relaxed", "c_normal")))
      # ---------------------------------------------------------------
  
      if (t == 1){ rtn <- data.frame(time = t-1, age_group = 1:9, foi = lambda)
      } else { 
        tmp <- data.frame(time = t-1, age_group = 1:9, foi = lambda)
        rtn <- bind_rows(rtn, tmp)}
    }
  
    rtn2 <- list(check = data.frame(time = time_vec - 1, 
                                criteria = criteria,
                                new_cases = cases,
                                ic_admin = ic_admin,
                                cm_check = cm_check),
                 lambda = rtn)
  
  return(rtn2)
}
