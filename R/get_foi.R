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
                    i1,
                    N,
                    c_lockdown = NULL,
                    c_relaxed = NULL,
                    ic_thresh_l,
                    ic_thresh_u
                    ){

    # sum over different I states for each time step and age group
    I_all <- as.matrix(dat$I + dat$Iv_1d + dat$Iv_2d)
    H_all <- as.matrix(dat$H + dat$Hv_1d + dat$Hv_2d)
    
    # calculate force of infection for each time point
    #foi <- t(apply(I_all, 1, function(x){beta * (c_main %*% (x/N))}))
    time_vec <- 1:dim(I_all)[1]
    ic_admin <- rep(0, length(time_vec))
    slope <- rep(0, length(time_vec))
    cm_check <- rep(NA, length(time_vec))
    
    for(t in time_vec){
      ic_admin[t] <- sum(i1 * H_all[t,])
      slope[t] <- ifelse( t == 1, 0, ic_admin[t] - ic_admin[t-1])
      
      contact_mat <- (ic_admin[t] < ic_thresh_u & ic_admin[t] >= ic_thresh_l & slope[t] < 0) * c_lockdown +
                     (ic_admin[t] < ic_thresh_l) * c_relaxed +
                     (ic_admin[t] >= ic_thresh_l & ic_admin[t] < ic_thresh_u & slope[t] > 0) * c_relaxed +
                     (ic_admin[t] >= ic_thresh_u) * c_lockdown 
      cm_check[t] <- ifelse(identical(contact_mat, c_lockdown), "c_lockdown", "c_relaxed")

      lambda <- beta * (contact_mat %*% (I_all[t,]/N))

      if (t == 1){ rtn <- data.frame(time = t-1, age_group = 1:9, foi = lambda)
      } else { 
        tmp <- data.frame(time = t-1, age_group = 1:9, foi = lambda)
        rtn <- bind_rows(rtn, tmp)}
    }
  
    rtn2 <- list(cp = data.frame(time = time_vec - 1, 
                                ic_admin = ic_admin,
                                slope = slope,
                                cm_check = cm_check),
                 lambda = rtn)
  
  return(rtn2)
}
