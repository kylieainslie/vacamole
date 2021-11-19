# ------------------------------------------------------------------
# Forward simulation function wrapper ------------------------------
# ------------------------------------------------------------------
#' Wrapper function to perform forward simulations
#' @param params List of model input parameters.
#' @param start_date Start date of simulations (character string in YYYY-MM-DD format)
#' @param end_date End date of simulations (character string in YYYY-MM-DD format)
#' @param init_cond Named vector of initial conditions for model compartments
#' @param beta_m Transmission rate mean.
#' @param vac_inputs Vaccination inputs for model. This should be the output of convert_vac_schedule().
#' @param beta_c Transmission rate when contact matrix changes to pre-pandemic (time > t_normal)
#' @param beta_d Transmission rate draws for determining uncertainty around the mean transmission rate
#' @param t_normal Time point to force model to switch to pre-pandemic contact matrix
#' @param contact_matrices List of contact matrices
#' @param tag Character string of name of output file
#' @return Data frame of model outputs for each simulation
#' @import tidyr
#' @import dplyr
#' @import lubridate
#' @import deSolve
#' @keywords vacamole
#' @export
forward_sim_func_wrap <- function(params,
                                  start_date,
                                  end_date,
                                  init_cond,
                                  beta_m,
                                  vac_inputs,
                                  beta_c,
                                  beta_d,
                                  t_normal,
                                  contact_matrices,
                                  tag) {

  # empty list for output
  out <- list()
  # specify time points ----------------------------------------------
  start_date <- lubridate::yday(as.Date(start_date) + 365) + 365
  end_date <- lubridate::yday(as.Date(end_date)) + (365 * 2)
  times <- seq(start_date, end_date, by = 1)

  initial_conditions <- c(t = times[1], init_cond)

  # Update parameter values for input into model solver ------
  params$beta <- beta_m
  #TODO if there's only one contact matrix then don't take the mean
  params$c_start <- contact_matrices$june_2021$mean
  params$c_lockdown <- contact_matrices$february_2021$mean
  params$c_relaxed <- contact_matrices$june_2020$mean
  params$c_very_relaxed <- contact_matrices$june_2021$mean
  params$c_normal <- contact_matrices$baseline_2017$mean
  params$vac_inputs <- vac_inputs
  params$beta_change <- beta_c[1]
  params$t_normal <- t_normal

  # if time doesn't start at 0 we need to initialise the contact
  # matrices flags ---------------------------------------------------
  assign("flag_relaxed", 0, envir = .GlobalEnv)
  assign("flag_very_relaxed", 0, envir = .GlobalEnv)
  assign("flag_normal", 0, envir = .GlobalEnv)

  #  Solve model ------------------------------------------------------
  # mle
  print("mle")
  seir_out <- lsoda(initial_conditions, times, age_struct_seir_ode, params)
  seir_out <- as.data.frame(seir_out)
  out_mle <- postprocess_age_struct_model_output(seir_out)
  tmp_mle <- lapply(out_mle, wide_to_long, times) %>%
    bind_rows() %>%
    mutate(sim = 0)

  # run for beta draws ------------------------------------------------
  n_beta <- length(beta_d)
  rtn_out <- list()

  for (i in 1:n_beta) {
    print(i)

    # reset flags
    # need to make sure they are saved to the global environment
    assign("flag_relaxed", 0, envir = .GlobalEnv)
    assign("flag_very_relaxed", 0, envir = .GlobalEnv)
    assign("flag_normal", 0, envir = .GlobalEnv)

    # change parameters
    params$beta <- beta_d[i]
    params$beta_c <- beta_c[i + 1]
    params$c_start <- contact_matrices$june_2021[[i]]
    params$normal <- contact_matrices$baseline_2017[[i]]

    # run model
    seir_out <- lsoda(initial_conditions, times, age_struct_seir_ode, params)
    seir_out <- as.data.frame(seir_out)
    out <- postprocess_age_struct_model_output(seir_out)

    # convert output from each simulation into long data frame
    tmp <- lapply(out, wide_to_long, times) %>%
      bind_rows() %>%
      mutate(sim = i)

    rtn_out[[i]] <- tmp
  }


  # return outouts
  rtn_out$mle <- tmp_mle
  rtn <- bind_rows(rtn_out)

  saveRDS(rtn, paste0("inst/extdata/results/", tag, ".rds"))

  return(rtn)
} # end of function
