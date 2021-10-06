# --------------------------------------------------
# likelihood function for model fit to case data
# --------------------------------------------------
#' @param x
#' @param contact_matrix
#' @param t
#' @param data
#' @param params
#' @param init
#' @param stochastic
#' @keywords vacamole
#' @export
likelihood_func <- function(x,
                            contact_matrix,
                            t,
                            data,
                            params,
                            init,
                            stochastic = FALSE) {
  #params$beta <- x[1] # pars["beta"]
  
  # params$beta1 <- beta1 #pars["beta1"]
  r0 <- x[1]
  S_diag <- diag(init[c(2:10)])
  rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  params$beta <- (r0 / rho) * params$gamma
  
  if (stochastic){
    seir_out <- stochastic_age_struct_seir_ode(times = t,init = init, params = params)
    out <- apply(seir_out, 3, rowSums)
    daily_cases <- params$sigma * (out[,"E"] + out[,"Ev_1d"] + out[,"Ev_2d"]) * params$p_report
    daily_cases <- ifelse(daily_cases == 0, 0.0001, daily_cases) # prevent likelihood function function from being Inf
  } else {
    seir_out <- lsoda(init,t,age_struct_seir_ode,params) #
    seir_out <- as.data.frame(seir_out)
    out <- postprocess_age_struct_model_output(seir_out)
    daily_cases <- (params$sigma * rowSums(out$E + out$Ev_1d + out$Ev_2d)) * params$p_report
    daily_cases <- ifelse(daily_cases == 0, 0.0001, daily_cases) # prevent likelihood function function from being Inf
  }
  
  inc_obs <- data$inc
  
  # lik <- sum(dpois(x = inc_obs,lambda = incidence,log=TRUE))
  alpha <- x[2]
  size <- daily_cases * (alpha/(1-alpha))
  lik <- -sum(dnbinom(x = inc_obs, mu = daily_cases, size = size, log = TRUE))
  lik
}