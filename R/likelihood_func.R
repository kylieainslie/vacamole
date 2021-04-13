# likelihood function for SEIR model to fit to real data
dat <- osiris1
tstep <- 1
weeks <- floor(dim(dat)[1]/7)
times <- seq(0,(weeks*7)-1,by=tstep)

likelihood_func <- function(pars){
  params$beta <- pars[1]
  #gamma <- pars[2]
  seir_out <- lsoda(init,times,age_struct_seir_ode,params) #
  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output(seir_out)
  
  # results <- as.data.frame(deSolve::ode(y=c(S=S0,I=I0,R=R0),
  #                                       times=t, func=SIR_odes,
  #                                       parms=c(beta,gamma)))
  # incidence <- c(0,-diff(results$S))
  daily_cases <- params$sigma * rowSums(out$E + out$Ev_1d + out$Ev_2d) * params$p_report
  
  ## Get weekly incidence
  weekly_cases <- colSums(matrix(daily_cases,nrow=7))
  
  lik <- sum(dpois(dat$inc,weekly_cases,log=TRUE))
  lik
}
real_pars <- c(params$beta)
test_pars <- c(0.06)
## Likelihood of true pars
print(likelihood_func(real_pars))
## Likelihood of test pars
print(likelihood_func(test_pars))
