#' convert cumulative vaccination schedule to non-cumulative ----------------------------------
#' @param breakpoints
#' @param params
#' @param init
#' @param fit_pars
#' @param case_data
#' @param contact_matrices
#' @param vac_info
#' @param save_output_to_file
#' @param path_out
#' @return list with output XX add more details XX
#' @keywords vacamole
#' @import tidyr
#' @import dplyr
#' @import lubridate
#' @import mvtnorm
#' @import rARPACK
#' @export
fit_to_data_func <- function(breakpoints, 
                             params, 
                             init,
                             fit_pars,
                             case_data,
                             contact_matrices,
                             vac_info,
                             save_output_to_file = TRUE,
                             path_out = NULL,
                             ...
                             ){
  

# create empty lists for storage ----------------------------
n_bp <- length(breakpoints$date)-1
mles <- matrix(rep(NA, 2*n_bp), nrow = n_bp)
colnames(mles) <- c("beta", "alpha")

out_mle <- list()
parameter_draws <- list()
beta_draws <- list()
daily_cases <- list()
susceptibles <- list()
exposed <- list()
infected <- list()
recovered <- list()
recovered1 <- list()
recovered2 <- list()
recovered3 <- list()

# begin loop over breakpoints --------------------------------
for (j in 1:n_bp) {
  
  print(paste("Fitting from", breakpoints$date[j], "to", breakpoints$date[j+1]))
  
  # set contact matrix for time window
  if (breakpoints$contact_matrix[j+1] == "april_2017"){contact_matrix <- contact_matrices$april_2017
  } else if (breakpoints$contact_matrix[j+1] == "april_2020"){contact_matrix <- contact_matrices$april_2020
  } else if (breakpoints$contact_matrix[j+1] == "june_2020"){contact_matrix <- contact_matrices$june_2020
  } else if (breakpoints$contact_matrix[j+1] == "september_2020"){contact_matrix <- contact_matrices$september_2020
  } else if (breakpoints$contact_matrix[j+1] == "february_2021"){contact_matrix <- contact_matrices$february_2021
  } else if (breakpoints$contact_matrix[j+1] == "june_2021"){contact_matrix <- contact_matrices$june_2021
  } else {contact_matrix <- contact_matrices$november_2021} 
  
  params$c_start <- contact_matrix
  
  if(!is.null(params$vac_inputs)){
    # set VE for time window depending on which variant was dominant
    if (breakpoints$variant[j+1] == "wildtype"){params$vac_inputs <- vac_info$wildtype
    } else if (breakpoints$variant[j+1] == "alpha"){params$vac_inputs <- vac_info$alpha
    } else if (breakpoints$variant[j+1] == "delta"){params$vac_inputs <- vac_info$delta
    } else {params$vac_inputs <- vac_info$omicron}
  }
  # set time sequence  
  times <- seq(breakpoints$time[j], breakpoints$time[j+1], by = 1)
  
  # update initial conditions based on last time window
  if (j == 1) {
    init_update <- init
    pars <- fit_pars$init_value
    S_diag <- diag(init_update[c(2:10)])
    rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  } else {
    init_update <- c(t = times[1], unlist(lapply(unname(out_mle[[j-1]]), tail,1)))
    beta_est <- (mles[j-1,1]/params$gamma)*rho
    pars <- c(beta_est, fit_pars$init_value[2]) # mles[j-1,-1]
    S_diag <- diag(init_update[c(2:10)])
    rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  }

  # subset data for time window
  case_data_sub <- case_data[times + 1, ]
  
  res <- optim(par = pars, 
               fn = likelihood_func2,
               method = "L-BFGS-B",
               lower = fit_pars$lower_bound,
               upper = fit_pars$upper_bound,
               t = times,
               data = case_data_sub,
               params = params,
               init = init_update,
               stochastic = FALSE,
               hessian = TRUE
  )
  
  # store MLE
  mles[j,1] <- (res$par[1] / rho) * params$gamma
  mles[j,2] <- res$par[2]
  
  print(mles[j,])
  # draw 200 parameter values
  parameter_draws[[j]] <- mvtnorm::rmvnorm(200, res$par, solve(res$hessian))
  beta_draws[[j]] <- data.frame(beta = (parameter_draws[[j]][,1] / rho) * params$gamma) %>%
    mutate(index = 1:200)
# --------------------------------------------------
  # run for mle to get initial conditions for next timepoint
  params$beta <- mles[j,1]

  seir_out <- lsoda(init_update, times, age_struct_seir_ode2, params)
  seir_out <- as.data.frame(seir_out)
  
  # store outputs
  out_mle[[j]] <- postprocess_age_struct_model_output2(seir_out)
  daily_cases[[j]] <- params$sigma * rowSums(out_mle[[j]]$E) * params$p_report #+ out_mle[[j]]$Ev_1d + out_mle[[j]]$Ev_2d + out_mle[[j]]$Ev_3d + out_mle[[j]]$Ev_4d + out_mle[[j]]$Ev_5d
  susceptibles[[j]] <- rowSums(out_mle[[j]]$S)
  exposed[[j]] <- rowSums(out_mle[[j]]$E)
  infected[[j]] <- rowSums(out_mle[[j]]$I)
  recovered[[j]] <- rowSums(out_mle[[j]]$R) 
  recovered1[[j]] <- rowSums(out_mle[[j]]$R_1w)
  recovered2[[j]] <- rowSums(out_mle[[j]]$R_2w)
  recovered3[[j]] <- rowSums(out_mle[[j]]$R_3w)

  # plot for quick check of fit
  plot(daily_cases[[j]]~times, type = "l")
  points(times, case_data_sub$inc, pch = 16, col = "red")
  
} # end of for loop over breakpoints

  todays_date <- Sys.Date()
  # path_out <- "/rivm/s/ainsliek/code/vacamole/inst/extdata/results/model_fits/"
# save outputs
  if(save_output_to_file){
    saveRDS(mles, file = paste0(path_out, "mles_from_fits_", todays_date, ".rds"))
    saveRDS(beta_draws, file = paste0(path_out, "beta_draws_from_fits_", todays_date, ".rds"))
    names(out_mle) <- paste0("end_date_", breakpoints$date) # name list elements for easier indexing
    saveRDS(out_mle, file = paste0(path_out, "output_from_fits_", todays_date, ".rds"))
  }
  rtn <- list(cases = unique(unlist(daily_cases)),
              ml_est = mles,
              beta_draws = beta_draws,
              out_mle = out_mle,
              susceptibles = unique(unlist(susceptibles)),
              recovered = unique(unlist(recovered)),
              recovered1 = unique(unlist(recovered1)),
              recovered2 = unique(unlist(recovered2)),
              recovered3 = unique(unlist(recovered3))
              )
  return(rtn)
} # end of function


# ----------------------------------------------------
# run simulations for mle, lower, and upper bounds 
# of beta
# ----------------------------------------------------
# fit_date <- "2021-10-01"
# beta_mles <- readRDS(paste0(path_out,"mles_from_fits_",fit_date,".rds"))
# beta_mles_list <- split(beta_mles, seq(nrow(beta_mles)))
# beta_draws <- readRDS(paste0(path_out,"beta_draws_from_fits_",fit_date,".rds"))

# run for 200 contact matrices
# mle_run <- model_run_wrapper(breakpoints = breakpoints, beta_values = beta_mles_list, init_conditions = init, params = params)
# ci_run  <- model_run_wrapper(breakpoints = breakpoints, beta_values = beta_draws, init_conditions = init, params = params, mle = FALSE)
# ci_out <- list()
# for (i in 1:n_bp){
#   ci_out[[i]] <- do.call("rbind", ci_run[[i]])
# }
# ci_out_wide <- do.call("cbind", ci_out)
# matplot(t(ci_out_wide), type = "l")
# 
# # try getting quantiles
# bounds <- apply(ci_out_wide, 2, quantile, probs = c(0.025, 0.975))
# matplot(t(bounds), type = "l")
# 
# # save outputs -------------------------------------
# # --------------------------------------------------
# #  combine all piecewise results to plot together
# cases_mle <- unique(unlist(mle_run))
# cases_lower <- unique(bounds[1,])
# cases_upper <- unique(bounds[2,])
# times_all <- 1:length(cases_mle)
# 
# model_fit <- data.frame(time = times_all, date = case_data$date, real = case_data$inc, mle = cases_mle, lower = cases_lower, upper = cases_upper)
# if(save_output_to_file){
#   saveRDS(model_fit, file = paste0(path_out, "model_fit_df_", todays_date, ".rds"))
# }
# return(model_fit)
# }
# --------------------------------------------------