#' Calculate force of infection 
#' @param dat data frame of states
#' @param beta transmission parameter
#' @param contact_matrix contact matrix
#' @param N vector of total group sizes
#' @return matrix of force of infection in each age group (columns) at each
#' time point (rows)
#' @keywords vacamole
#' @export
get_foi <- function(dat, 
                    beta,
                    sigma,
                    i1,
                    p_report,
                    N,
                    c_lockdown = NULL,
                    c_relaxed = NULL,
                    c_very_relaxed = NULL,
                    c_normal = NULL,
                    thresh_l,
                    thresh_m,
                    thresh_u,
                    use_cases = TRUE,
                    force_relax = NULL
                    ){

    # sum over different I states for each time step and age group
    I_all <- as.matrix(dat$I + dat$Iv_1d + dat$Iv_2d)
    H_all <- as.matrix(dat$H + dat$Hv_1d + dat$Hv_2d)
    E_all <- as.matrix(dat$E + dat$Ev_1d + dat$Ev_2d)
    
    # calculate force of infection for each time point
    #foi <- t(apply(I_all, 1, function(x){beta * (c_main %*% (x/N))}))
    time_vec <- 1:dim(I_all)[1]
    ic_admin <- rep(0, length(time_vec))
    cases    <- rep(0, length(time_vec))
    criteria <- rep(0, length(time_vec))
    slope    <- rep(0, length(time_vec))
    cm_check <- rep(NA, length(time_vec))
    
    for(t in time_vec){
      ic_admin[t] <- sum(i1 * H_all[t,])
      cases[t] <- sum(sigma * E_all[t,] * p_report)
      
      criteria[t] <- (use_cases) * cases[t] + (!use_cases) * ic_admin[t] 
      slope[t] <- ifelse( t == 1, 0, criteria[t] - criteria[t-1])
      
      contact_mat <- (is.null(force_relax)) *
        ((criteria[t] < thresh_u & criteria[t] >= thresh_m & slope[t] < 0) * c_lockdown +
           (criteria[t] < thresh_m & criteria[t] >= thresh_l & slope[t] < 0) * c_relaxed +
           (criteria[t] < thresh_l) * c_very_relaxed +
           (criteria[t] >= thresh_l & criteria[t] <= thresh_m & slope[t] > 0) * c_very_relaxed +
           (criteria[t] > thresh_m & criteria[t] <= thresh_u & slope[t] > 0) * c_relaxed +
           (criteria[t] > thresh_u) * c_lockdown +
           (criteria[t] >= 0 & criteria[t] < thresh_l & slope[t] > 0) * c_normal) +
        (!is.null(force_relax) & (t < force_relax)) *
        ((criteria[t] < thresh_u & criteria[t] >= thresh_m & slope[t] < 0) * c_lockdown +
           (criteria[t] < thresh_m & criteria[t] >= thresh_l & slope[t] < 0) * c_relaxed +
           (criteria[t] < thresh_l) * c_very_relaxed +
           (criteria[t] >= thresh_l & criteria[t] <= thresh_m & slope[t] > 0) * c_very_relaxed +
           (criteria[t] > thresh_m & criteria[t] <= thresh_u & slope[t] > 0) * c_relaxed +
           (criteria[t] > thresh_u) * c_lockdown +
           (criteria[t] >= 0 & criteria[t] < thresh_l & slope[t] > 0) * c_normal) +
        (!is.null(force_relax) & (t >= force_relax)) *
        ((criteria[t] < thresh_u & criteria[t] >= thresh_m & slope[t] < 0) * c_relaxed +
           (criteria[t] < thresh_m & criteria[t] >= thresh_l & slope[t] < 0) * c_relaxed +
           (criteria[t] < thresh_l) * c_very_relaxed +
           (criteria[t] >= thresh_l & criteria[t] <= thresh_m & slope[t] > 0) * c_very_relaxed +
           (criteria[t] > thresh_m & criteria[t] <= thresh_u & slope[t] > 0) * c_relaxed +
           (criteria[t] > thresh_u) * c_relaxed +
           (criteria[t] >= 0 & criteria[t] < thresh_l & slope[t] > 0) * c_normal)
      
      # contact_mat <- (criteria[t] < thresh_u & criteria[t] >= thresh_m & slope[t] < 0) * c_lockdown +
      #   (criteria[t] < thresh_m & criteria[t] >= thresh_l & slope[t] < 0) * c_relaxed +
      #   (criteria[t] < thresh_l) * c_very_relaxed +
      #   (criteria[t] >= thresh_l & criteria[t] <= thresh_m & slope[t] > 0) * c_very_relaxed +
      #   (criteria[t] > thresh_m & criteria[t] <= thresh_u & slope[t] > 0) * c_relaxed +
      #   (criteria[t] > thresh_u) * c_lockdown +
      #   (criteria[t] >= 0 & criteria[t] < thresh_l & slope[t] > 0) * c_normal
      
      # contact_mat <- (t < 20) * c_lockdown +
      #   (criteria[t] < thresh_u & criteria[t] >= thresh_l & slope[t] < 0 & t >=20) * c_lockdown +
      #   (criteria[t] < thresh_l & t >=20) * c_relaxed +
      #   (criteria[t] >= thresh_l & criteria[t] <= thresh_u & slope[t] > 0 & t >=20) * c_relaxed +
      #   (criteria[t] > thresh_u & t >=20) * c_lockdown 
      
      cm_check[t] <- ifelse(identical(contact_mat, c_lockdown), "c_lockdown", 
                            ifelse(identical(contact_mat, c_relaxed),"c_relaxed",
                                   ifelse(identical(contact_mat, c_very_relaxed),"c_very_relaxed", "c_normal")))

      lambda <- beta * (contact_mat %*% I_all[t,])

      if (t == 1){ rtn <- data.frame(time = t-1, age_group = 1:9, foi = lambda)
      } else { 
        tmp <- data.frame(time = t-1, age_group = 1:9, foi = lambda)
        rtn <- bind_rows(rtn, tmp)}
    }
  
    rtn2 <- list(check = data.frame(time = time_vec - 1, 
                                criteria = criteria,
                                slope = slope,
                                cm_check = cm_check),
                 lambda = rtn)
  
  return(rtn2)
}
