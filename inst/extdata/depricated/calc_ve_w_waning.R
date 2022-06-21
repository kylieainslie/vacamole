#' Calculate amount of vaccine waning
#' @param vac_rate rate of vaccination at each timepoint
#' @param ve_val value of VE
#' @param waning vector of amount of waning over time
#' @return matrix of VE values with waning in each age group (columns) at each
#' time point (rows)
#' @keywords vacamole
#' @export
calc_ve_w_waning <- function(vac_rate, ve_val, waning) {
  vac_rate <- as.matrix(vac_rate)
  waning_tot <- matrix(, nrow = nrow(vac_rate), ncol = ncol(vac_rate))
  ve_tot <- matrix(, nrow = nrow(vac_rate), ncol = ncol(vac_rate))
  for (t_tot in 1:nrow(vac_rate)) { #
    waning_t <- matrix(, nrow = t_tot, ncol = ncol(vac_rate))
    for (t in 1:t_tot) {
      waning_t[t, ] <- ifelse(t_tot - t > 0, waning[t_tot - t], 0) * vac_rate[t, ]
    }
    waning_tot[t_tot, ] <- apply(waning_t, 2, sum)
    ve_tot[t_tot, ] <- ve_val - (ve_val * waning_tot[t_tot, ])
  }
  return(ve_tot)
}
