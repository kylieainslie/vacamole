#' Calculate force of infection
#' @param x output from SEIR model with number of individuals in each compartment at each time point
#' @param y1 data frame of protection against transmission due to vaccination from dose 1 at each time point
#' @param y2 data frame of protection against transmission due to vaccination from dose 2 at each time point
#' @param y3 data frame of protection against transmission due to vaccination from dose 3 at each time point
#' @param y4 data frame of protection against transmission due to vaccination from dose 4 at each time point
#' @param y5 data frame of protection against transmission due to vaccination from dose 5 at each time point
#' @param beta vector of values of transmission rate for each time point
#' @param contact_mat contact matrix
#' @param times vector of time points that correspond to the rows of x
#' @return matrix of force of infection in each age group (columns) at each
#' time point (rows)
#' @keywords vacamole
#' @export
get_foi <- function(x, y1, y2, y3, y4, y5, beta, contact_mat, times){
  tmp <- list()
  for(t in 1:length(times)){
    tmp[[t]] <- t(beta[t] * (contact_mat %*% (unlist(x[t,c(paste0("I",1:9))]) + 
                                                (unlist(y1[t,]) * unlist(x[t, c(paste0("Iv_1d",1:9))])) + 
                                                (unlist(y2[t,]) * unlist(x[t, c(paste0("Iv_2d",1:9))])) + 
                                                (unlist(y3[t,]) * unlist(x[t, c(paste0("Iv_3d",1:9))])) + 
                                                (unlist(y4[t,]) * unlist(x[t, c(paste0("Iv_4d",1:9))])) + 
                                                (unlist(y5[t,]) * unlist(x[t, c(paste0("Iv_5d",1:9))]))))
    )
  }
  
  rtn <- do.call(rbind, tmp)
  return(rtn)
  
}