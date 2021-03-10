#' Determine contact matrix based on thresholds of cases or IC admissions
#' @param params list of parameter values
#' @param criteria 
#' @param slope calendar date of start of simulation
#' @return List of summary results
#' @keywords vacamole
#' @export

choose_contact_matrix <- function(params, times, criteria, slope){
  # define variables from params
  force_relax <- params$force_relax
  thresh_l <- params$thresh_l
  thresh_m <- params$thresh_m
  thresh_u <- params$thresh_u
  c_start <- params$c_start
  c_lockdown <- params$c_lockdown
  c_relaxed <- params$c_relaxed
  c_very_relaxed <- params$c_very_relaxed
  c_normal <- params$c_normal
  t_start_end <- params$t_start_end
  
  if(times <= t_start_end){
    contact_mat <- c_start
  } else if(is.null(force_relax)){
    contact_mat <- (criteria >= thresh_u) * c_lockdown +
      (criteria < thresh_u & criteria >= thresh_m & slope < 0) * c_lockdown +
      (criteria < thresh_m & criteria >= thresh_l & slope < 0) * c_relaxed +
      (criteria < thresh_l & slope < 0) * c_very_relaxed +
      (criteria >= thresh_l & criteria <= thresh_m & slope > 0) * c_very_relaxed +
      (criteria > thresh_m & criteria <= thresh_u & slope > 0) * c_relaxed +
      (criteria >= 0 & criteria < thresh_l & slope > 0) * c_normal
    } else{
      if(t < force_relax){
        contact_mat <- (criteria > thresh_u) * c_lockdown +
          (criteria < thresh_u & criteria >= thresh_m & slope < 0) * c_lockdown +
          (criteria < thresh_m & criteria >= thresh_l & slope < 0) * c_relaxed +
          (criteria < thresh_l & slope < 0) * c_very_relaxed +
          (criteria >= thresh_l & criteria <= thresh_m & slope > 0) * c_very_relaxed +
          (criteria > thresh_m & criteria <= thresh_u & slope > 0) * c_relaxed +
          (criteria >= 0 & criteria < thresh_l & slope > 0) * c_normal
      } else {
        contact_mat <- (criteria > thresh_u) * c_relaxed +
          (criteria < thresh_u & criteria >= thresh_m & slope < 0) * c_relaxed +
          (criteria < thresh_m & criteria >= thresh_l & slope < 0) * c_relaxed +
          (criteria < thresh_l & slope < 0) * c_very_relaxed +
          (criteria >= thresh_l & criteria <= thresh_m & slope > 0) * c_very_relaxed +
          (criteria > thresh_m & criteria <= thresh_u & slope > 0) * c_relaxed +
          (criteria >= 0 & criteria < thresh_l & slope > 0) * c_normal
      }
    }
  
  if(identical(contact_mat, c_lockdown)){ cat("criteria: ", criteria, "contact matrix: c_lockdown", "\n")
  } else if(identical(contact_mat, c_relaxed)){cat("criteria: ", criteria, "contact matrix: c_relaxed", "\n")
  } else if (identical(contact_mat, c_very_relaxed)){ cat("criteria: ", criteria, "contact matrix: c_very_relaxed", "\n")
  } else if (identical(contact_mat, c_normal)){cat("criteria: ", criteria, "contact matrix: c_normal", "\n")
  } else {cat("criteria: ", criteria, "contact matrix: c_start", "\n")}


  return(contact_mat)
}