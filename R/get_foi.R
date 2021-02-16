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
                    c_main, 
                    N 
                    # constant_contact_matrix = TRUE,
                    # c_lockdown = NULL,
                    # c_relaxed = NULL,
                    # times = NULL,
                    # eta = NULL,
                    # eta2 = NULL,
                    # one_cp = FALSE
                    ){

  # if(!constant_contact_matrix){
    # sum over different I states for each time step and age group
    I_all <- dat$I + dat$Iv_1d + dat$Iv_2d
    
    # calculate force of infection for each time point
    foi <- t(apply(I_all, 1, function(x){beta * (c_main %*% (x/N))}))
    
  # } else{
  #   rtn <- list()
  #   for(t in 1:length(times)){
  #     if(t == 1){ 
  #       C <- c_main 
  #       lambda <- beta * c_main %*% as.numeric((dat$I[t,] + dat$Iv_1d[t,] + dat$Iv_2d[t,])/N)
  #     } else {
  #     upper_thresh <- sum(N) * 21/100000
  #     lower_thresh <- sum(N) * 7/100000
  #     log_cm <- all.equal(C, c_lockdown)
  #   
  #     incidence <- (dat$S[t,] + dat$Shold_1d[t,] + (eta[t,] * (dat$Sv_1d[t,] + dat$Shold_2d[t,])) + 
  #                     (eta2[t,] * dat$Sv_2d[t,])) * lambda
  #     s_inc <- sum(incidence)
  #     if(one_cp){
  #       flag <- ifelse(sum(incidence) < lower_thresh & log_cm == TRUE, 0, 1)
  #     } else{
  #     flag <- ifelse((s_inc < upper_thresh & s_inc > lower_thresh & log_cm == TRUE) |
  #                    (sum(incidence) > upper_thresh & log_cm == FALSE) |
  #                    (sum(incidence) > upper_thresh & log_cm == TRUE), 1, 0)
  #     }    
  #     
  #     if(flag == 0){ C <- c_relaxed
  #     } else { C <- c_lockdown}
  #     lambda <- beta * C %*% as.numeric((dat$I[t,] + dat$Iv_1d[t,] + dat$Iv_2d[t,])/N)
  #     }
  #     tmp <- as.vector(lambda)
  #     names(tmp) <- names(dat$I)
  #     rtn[[t]] <- tmp
  #   }
  #   foi <- bind_rows(rtn)
  # }
  
  return(foi)
}
