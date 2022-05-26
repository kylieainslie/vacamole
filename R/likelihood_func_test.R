# --------------------------------------------------
# likelihood function for model fit to case data
# --------------------------------------------------
#' Likelihood function for fitting SEIR model to daily case data
#' @param x vector of parameters to fit
#' @param t vector of time points
#' @param data real data to fit model to
#' @param params list of parameters values for input into SEIR model
#' @param init named vector of initial conditions for each model compartment
#' @param stochastic logical, if TRUE, a stochastic model is fit to the data
#' otherwise a deterministic model is used.
#' @keywords vacamole
#' @importFrom stats dnbinom
#' @importFrom rARPACK eigs
#' @export
likelihood_func_test <- function(x,
                             t,
                             data,
                             params,
                             init,
                             model_func,
                             ...) {
  
  r0 <- x[1]
  S_diag <- diag(init[c(2:10)])
  rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  params$beta <- (r0 / rho) * params$gamma
  
  seir_out <- lsoda(init, t, model_func, params, rtol = 0.00001, hmax = 0.02) # hmax = 0.02, age_struct_seir_ode_test
  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output2(seir_out)
  daily_cases <- (params$sigma * rowSums(out$E)) * params$p_report
  daily_cases <- ifelse(daily_cases == 0, 0.0001, daily_cases) # prevent likelihood function function from being Inf
  
  inc_obs <- data$inc
  
  # lik <- sum(dpois(x = inc_obs,lambda = incidence,log=TRUE))
  alpha <- x[2]
  #size <- daily_cases * (alpha / (1 - alpha))
  lik <- -sum(stats::dnbinom(x = inc_obs, mu = daily_cases, size = alpha, log = TRUE))
  
  #print(x)
  #print(lik)
  
  lik
}
