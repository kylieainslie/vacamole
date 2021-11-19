#' Function wrapper to run model for multiple time windows with set of
#' transmission rates
#' @param breakpoints list of breakpoints
#' @param beta_values vector of beta values
#' @param init_conditions vector of initial conditions
#' @param params vector of model input parameters
#' @param mle logical, if TRUE the current run is the maximum likelihood estimate
#' @keywords vacamole
#' @importFrom lubridate yday
#' @importFrom utils tail
#' @export
model_run_wrapper <- function(breakpoints,
                              beta_values,
                              init_conditions,
                              params,
                              mle = TRUE) {
  n_bp <- length(breakpoints$date)
  out_store <- list()
  daily_cases_store <- list()

  for (j in 1:n_bp) {
    print(j)
    out_store[[j]] <- list()
    daily_cases_store[[j]] <- list()
    # set time vector -------------------------------------
    if (j == 1) {
      # if first time window, start time at 0
      end_day <- yday(breakpoints$date[j]) - 1

      times <- seq(0, end_day, by = 1)
    } else {
      if (breakpoints$indicator_2021[j] == 1) {
        if (breakpoints$indicator_2021[j - 1] == 1) { # wait for two consecutive dates in 2021
          start_day <- yday(breakpoints$date[j - 1]) - 1 + 366 # shift days by 1 because we start time at 0 (not 1)
        } else {
          start_day <- yday(breakpoints$date[j - 1]) - 1
        }
        end_day <- yday(breakpoints$date[j]) - 1 + 366
      } else {
        start_day <- yday(breakpoints$date[j - 1]) - 1
        end_day <- yday(breakpoints$date[j]) - 1
      }
      times <- seq(start_day, end_day, by = 1)
    }

    # loop over beta values (only 1 if mle = TRUE)
    n_beta <- ifelse(mle, 1, dim(beta_values[[j]][1]))

    for (i in 1:n_beta) {

      # beta ------------------------------------------------
      params$beta <- ifelse(mle, beta_values[[j]][1], beta_values[[j]][i, 1])

      # contact matrices ------------------------------------
      if (mle) {
        if (j == n_bp) {
          params$c_start <- breakpoints$contact_matrix[[j - 1]]$mean
        } else {
          params$c_start <- breakpoints$contact_matrix[[j]]$mean
        }
      } else {
        if (j == n_bp) {
          params$c_start <- matrix(unlist(breakpoints$contact_matrix[[j - 1]][i]), nrow = 9)
        } else {
          params$c_start <- matrix(unlist(breakpoints$contact_matrix[[j]][i]), nrow = 9)
        }
      }

      if (j == 1) {
        # set initial conditions ----------------------------
        init_update <- init_conditions
      } else {
        # update initial conditions --------------------------
        init_update <- c(t = times[1], unlist(lapply(unname(out_store[[j - 1]][[i]]), utils::tail, 1)))
      }
      # --------------------------------------------------
      # run model
      # --------------------------------------------------
      seir_out <- lsoda(init_update, times, age_struct_seir_ode, params)
      seir_out <- as.data.frame(seir_out)
      out_store[[j]][[i]] <- postprocess_age_struct_model_output(seir_out)
      daily_cases_store[[j]][[i]] <- params$sigma * rowSums(out_store[[j]][[i]]$E + out_store[[j]][[i]]$Ev_1d + out_store[[j]][[i]]$Ev_2d) * params$p_report
      # --------------------------------------------------
    } # end loop over beta values
  } # end of for loop over breakpoints

  return(daily_cases_store)
}
