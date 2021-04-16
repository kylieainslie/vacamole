# Fit model to OSIRIS data using MCMC
library(lazymcmc)
# read in OSIRIS data
osiris <- readRDS("inst/extdata/data/osiris_20210409.rds")
osiris1 <- osiris %>%
  group_by(date) %>%
  summarise_at(.vars = "n", .funs = "sum") %>%
  filter(date >= as.Date("2021-01-31")) %>%
  rename(inc = n)
  # complete(., date = full_seq(date, 1)) %>%
  # mutate(n = ifelse(is.na(n), 0, n))

# specify params for MCMC
parTab <- data.frame(names=c("beta"),
                     values=c(0.0006),
                     fixed=c(0),
                     steps=c(0.0001),
                     lower_bound=c(0),
                     upper_bound=c(1))

mcmcPars <- c("iterations"=2000,"popt"=0.44,"opt_freq"=100,
              "thin"=1,"adaptive_period"=1000,"save_block"=100)


## Putting model solving code in a function for later use
solve_model <- function(pars, t, params, age_struct_seir_ode){

  params$beta <- pars["beta"]
  seir_out <- lsoda(init,t,age_struct_seir_ode,params) #
  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output(seir_out)
  
  daily_cases <- params$sigma * rowSums(out$E + out$Ev_1d + out$Ev_2d) * params$p_report
  
  return(daily_cases)
}

## We need to put our likelihood function in a closure environment
## It's important to write the function in this form!
create_lik <- function(parTab, data, PRIOR_FUNC,...){
  par_names <- parTab$names
  
  ## Extract observed incidence
  inc <- data$inc
  
  ## Get times to solve model over
  tstep <- 1
  weeks <- floor(dim(data)[1]/7)
  t <- seq(0,(weeks*7)-1,by=tstep)
  
  ## We will pass S0, I0, R0 and SIR_odes 
  #N <- S0 + I0 + R0
  ## using the `...` bit.
  likelihood_func <- function(pars){
    names(pars) <- par_names
    incidence <- solve_model(pars, t, params, age_struct_seir_ode)
    ## Get weekly incidence
    #weekly_incidence <- colSums(matrix(incidence,nrow=7/tstep))
    
    lik <- sum(dpois(x = inc,lambda = incidence,log=TRUE))
    if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
    lik
    
  }
}

## Starting points, chosen pseudo at random
## seeding chains for SIR models is hard with deSolve,
## so I've chosen points near the true values.
startTab <- parTab
startTab$values <- c(0.0006)
tstep <- 1
weeks <- floor(dim(osiris1)[1]/7)
t <- seq(0,(weeks*7)-1,by=tstep)
osiris2 <- osiris1[t+1,]

output <- run_MCMC(parTab=startTab, data=osiris2, mcmcPars=mcmcPars, filename="SEIR_fit_no_prior",
                   CREATE_POSTERIOR_FUNC = create_lik, mvrPars = NULL, PRIOR_FUNC=NULL,
                   params = params, age_struct_seir_ode = age_struct_seir_ode)
#chain <- read.csv(output$file)
chain <- read.csv("SEIR_fit_no_prior_univariate_chain.csv")
plot(coda::as.mcmc(chain[,c("beta")]))

chain1 <- chain[chain$sampno >= 200, ] # mcmcPars["adaptive_period"] ,2:(ncol(chain)-1)

## Prior function
# my_prior <- function(pars){
#   ## Diffuse gaussian priors on both parameters
#   beta_prior <- dgamma(pars[1],0.001,1000)
#   #gamma_prior <- dnorm(pars[2],0.2,100,1)
#   return(beta_prior)
# }

## Use the previous chain to get a good starting
## covariance matrix
# startTab$values <- get_best_pars(chain)
# chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
# covMat <- cov(chain)
# mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)
# 
# ## Remove restrictions on parameter bounds -
# ## not important with the multivariate sampler
# startTab$lower_bound <- -Inf
# startTab$upper_bound <- Inf
# 
# ## Run with mvrPars
# mcmcPars <- c("iterations"=20000,"popt"=0.234,"opt_freq"=500,
#               "thin"=1,"adaptive_period"=7500,"save_block"=100)
# 
# output2 <- run_MCMC(parTab=startTab, data=dat, mcmcPars=mcmcPars, filename="SIR_fitting",
#                     CREATE_POSTERIOR_FUNC = create_lik, mvrPars = mvrPars, PRIOR_FUNC=my_prior,
#                     S0=999,I0=1,R0=0,SIR_odes=SIR_odes)
# chain2 <- read.csv(output2$file)
# chain2 <- chain2[chain2$sampno >= mcmcPars["adaptive_period"],]
# plot(coda::as.mcmc(chain2[,c("R0","gamma")]))

## Get the maximum likelihood parameters and estimate
## the model trajectory for each week
best_pars <- get_best_pars(chain1)
times <- seq(1,weeks*7,by=1)
best_trajectory <- solve_model(best_pars, times, params, age_struct_seir_ode)
best_incidence <- colSums(matrix(best_trajectory,nrow=7))
best_dat <- data.frame(t=seq(1,weeks,by=1), y=best_incidence)

## Bit of code to generate prediction intervals
n_samps <- 100 #100
trajectories <- matrix(nrow=n_samps,ncol=weeks)
samps <- sample(nrow(chain1),n_samps)
for(i in 1:n_samps){
  print(i)
  pars <- as.numeric(chain1[samps[i],c("beta")])
  names(pars) <- c("beta")
  tmp <- solve_model(pars, t, params, age_struct_seir_ode)
  trajectories[i,] <- colSums(matrix(tmp, nrow=7))
}
bounds <- apply(trajectories,2,function(x) quantile(x, c(0.025,0.975)))

osiris_weekly_inc <- colSums(matrix(osiris2$inc,nrow=7))
plot(osiris_weekly_inc~seq(1,weeks,by=1),col="red",pch=16, 
     xlab="Time (weeks)",ylab="Incidence", ylim=c(0,42000)) 
lines(best_dat,col="blue")
lines(bounds[1,],col="blue",lty=2,lwd=0.5)
lines(bounds[2,],col="blue",lty=2,lwd=0.5)
legend("topright",c("Data","Model","95% credible intervals"),
       col=c("red","blue","blue"),lty=c(0,1,2),pch=c(16,NA,NA))

