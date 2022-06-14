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
                             ...) {
  
  r0 <- x[1]
  S_diag <- diag(init[c(2:10)])
  rho <- as.numeric(eigs(S_diag %*% params$contact_mat, 1)$values)
  params$beta <- (r0 / rho) * mean(params$gamma)
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init, times, age_struct_seir_ode_test, params, method = rk45, rtol = 1e-08, hmax = 0.02) # , rtol = 1e-08, hmax = 0.02
  out <- as.data.frame(seir_out)
  out1 <- postprocess_age_struct_model_output2(out)
  daily_cases <- rowSums(params$sigma * out1$E * params$p_report)
  daily_cases <- ifelse(daily_cases == 0, 0.0001, daily_cases) # prevent likelihood function function from being Inf
  
  inc_obs <- data$inc
  
  # lik <- sum(dpois(x = inc_obs,lambda = incidence,log=TRUE))
  alpha <- x[2]
  #size <- daily_cases * (alpha / (1 - alpha))
  lik <- -sum(stats::dnbinom(x = inc_obs, mu = daily_cases, size = alpha, log = TRUE))
  
  # print(x)
  # print(lik)
  
  lik
}
