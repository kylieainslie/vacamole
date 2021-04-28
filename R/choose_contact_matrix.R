#' Determine contact matrix based on thresholds of cases or IC admissions
#' @param params list of parameter values
#' @param criteria 
#' @param slope calendar date of start of simulation
#' @return List of summary results
#' @keywords vacamole
#' @export

choose_contact_matrix <- function(params, times, criteria, flag_relaxed, 
                                  flag_very_relaxed, flag_normal){
  # define variables from params
  thresh_n <- params$thresh_n
  thresh_l <- params$thresh_l
  thresh_m <- params$thresh_m
  thresh_u <- params$thresh_u
  #c_start <- params$c_start
  c_lockdown <- params$c_lockdown
  c_relaxed <- params$c_relaxed
  c_very_relaxed <- params$c_very_relaxed
  c_normal <- params$c_normal
  #t_start_end <- params$t_start_end
  
  # if(!is.null(t_start_end) & times <= t_start_end){
  #   contact_matrix <- c_start
  # } else{
  #   # use simpler conditions where measures are only relaxed and not re-tightened
    # for flags
      if(criteria >= thresh_u){ 
        flag_relaxed <- 0
        flag_very_relaxed <- 0
        flag_normal <- 0
      }
      if(criteria <= thresh_m){flag_relaxed <- flag_relaxed + 1}
      if(criteria <= thresh_l){flag_very_relaxed <- flag_very_relaxed + 1}
      if(criteria <= thresh_n){flag_normal <- flag_normal + 1 }
    # for contact matrix
      if (flag_relaxed > 0 & flag_very_relaxed == 0) {contact_matrix <- c_relaxed
      } else if (flag_very_relaxed > 0 & flag_normal == 0) { contact_matrix <- c_very_relaxed
      } else if( flag_normal > 0 ){ contact_matrix <- c_normal
      } else {contact_matrix <- c_lockdown}
  # } 
  
  # cat("flag_relaxed: ", flag_relaxed, "\n")
  # cat("flag_very_relaxed: ", flag_very_relaxed, "\n")
  # cat("flag_normal: ", flag_normal, "\n")
  # if(identical(contact_matrix, c_lockdown)){ cat("criteria: ", criteria, "contact matrix: c_lockdown", "\n")
  # } else if(identical(contact_matrix, c_relaxed)){cat("criteria: ", criteria, "contact matrix: c_relaxed", "\n")
  # } else if (identical(contact_matrix, c_very_relaxed)){ cat("criteria: ", criteria, "contact matrix: c_very_relaxed", "\n")
  # } else if (identical(contact_matrix, c_normal)){cat("criteria: ", criteria, "contact matrix: c_normal", "\n")
  # } else {cat("criteria: ", criteria, "contact matrix: c_start", "\n")}

  rtn <- list(contact_matrix = contact_matrix,
              flag_relaxed = flag_relaxed,
              flag_very_relaxed = flag_very_relaxed,
              flag_normal = flag_normal)
  return(rtn)
}