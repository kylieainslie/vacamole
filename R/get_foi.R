#' Calculate force of infection 
#' @param dat data frame of states
#' @param params list of parameter values
#' @param contact_matrix contact matrix
#' @return matrix of force of infection in each age group (columns) at each
#' time point (rows)
#' @keywords vacamole
#' @export
get_foi <- function(dat, params){

    
    # sum over different I states for each time step and age group
    # S <- as.matrix(dat$S)
    # Shold_1d <- as.matrix(dat$Shold_1d)
    # Sv_1d <- as.matrix(dat$Sv_1d)
    # Shold_2d <- as.matrix(dat$Shold_2d)
    # Sv_2d <- as.matrix(dat$Sv_2d)
    E_all <- as.matrix(dat$E + dat$Ev_1d + dat$Ev_2d)
    I_all <- as.matrix(dat$I + dat$Iv_1d + dat$Iv_2d)
    H_all <- as.matrix(dat$H + dat$Hv_1d + dat$Hv_2d)
    
    
    # calculate force of infection for each time point
    time_vec <- 1:dim(E_all)[1]
    ic_admin <- rep(0, length(time_vec))
    cases    <- rep(0, length(time_vec))
    criteria <- rep(0, length(time_vec))
    #slope    <- rep(0, length(time_vec))
    cm_check <- rep(NA, length(time_vec))
    flag_relaxed <- 0
    flag_very_relaxed <- 0
    flag_normal <- 0
    
    for(t in time_vec){
      ic_admin[t] <- sum(params$i1 * H_all[t,])
      cases[t] <- sum(params$sigma * E_all[t,] * params$p_report)
      
      criteria[t] <- (params$use_cases) * cases[t] + (!params$use_cases) * ic_admin[t] 
      
      #if(t == 0 ){cases <- sum(init_lambda * ((S[t,] + Shold_1d[t,]) + eta * (Sv_1d[t,] + Shold_2d[t,]) + eta2 * Sv_2d[t,]))}
      #criteria <- (use_cases) * cases + (!use_cases) * ic_admin 
      
      tmp2 <- choose_contact_matrix(params, t, criteria[t], flag_relaxed, 
                                    flag_very_relaxed, flag_normal)
      contact_mat <- tmp2$contact_matrix
      flag_relaxed <- tmp2$flag_relaxed 
      flag_very_relaxed <- tmp2$flag_very_relaxed 
      flag_normal <- tmp2$flag_normal
      
      # determine force of infection ----------------------------------
      lambda <- params$beta * params$delta * (contact_mat %*% (I_all[t,]))
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
