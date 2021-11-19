#' Determine contact matrix based on thresholds of cases or IC admissions
#' @param times vector of time points
#' @param params list of parameter values
#' @param criteria criteria by which to change contact matrix. User can choose between a number of
#' cases per day or a number of IC admissions.
#' @param flag_relaxed an integer value that becomes non-zero (positive) when the contact matrix consistent with
#' "relaxed" non-pharmaceutical interventions is being used.
#' @param flag_very_relaxed an integer value that becomes non-zero (positive) when the contact matrix consistent with
#' "very relaxed" non-pharmaceutical interventions is being used.
#' @param flag_normal an integer value that becomes non-zero (positive) when the contact matrix consistent with
#' normal (pre-pandemic) contact patterns is being used.
#' @param keep_fixed logical. if TRUE the contact matrix stays fixed over the
#' entire simulation period
#' @return List of summary results
#' @keywords vacamole
#' @export

choose_contact_matrix <- function(times,
                                  params,
                                  criteria,
                                  flag_relaxed,
                                  flag_very_relaxed,
                                  flag_normal,
                                  keep_fixed) {
  # define variables from params
  c_start <- params$c_start
  if (!is.null(c_start) & keep_fixed) {
    contact_matrix <- c_start
  } else {
    thresh_n <- params$thresh_n
    thresh_l <- params$thresh_l
    thresh_m <- params$thresh_m
    thresh_u <- params$thresh_u

    c_lockdown <- params$c_lockdown
    c_relaxed <- params$c_relaxed
    c_very_relaxed <- params$c_very_relaxed
    c_normal <- params$c_normal
    t_normal <- params$t_normal

    # use simpler conditions where measures are only relaxed and not re-tightened
    # for flags
    if (criteria >= thresh_u) {
      flag_relaxed <- 0
      flag_very_relaxed <- 0
      flag_normal <- 0
    }
    if (criteria <= thresh_m) {
      flag_relaxed <- flag_relaxed + 1
    }
    if (criteria <= thresh_l) {
      flag_very_relaxed <- flag_very_relaxed + 1
    }
    if (criteria <= thresh_n) {
      if (is.null(t_normal)) {
        flag_normal <- flag_normal + 1
      } else if (times < t_normal) {
        flag_normal <- 0
      } else {
        flag_normal <- flag_normal + 1
      }
    }
    # for contact matrix
    if (flag_relaxed > 0 & flag_very_relaxed == 0) {
      contact_matrix <- c_relaxed
    } else if (flag_very_relaxed > 0 & flag_normal == 0) {
      contact_matrix <- c_very_relaxed
    } else if (flag_normal > 0) {
      contact_matrix <- c_normal
    } else {
      contact_matrix <- c_start
    }
  }

  rtn <- list(
    contact_matrix = contact_matrix,
    flag_relaxed = flag_relaxed,
    flag_very_relaxed = flag_very_relaxed,
    flag_normal = flag_normal
  )
  return(rtn)
}
