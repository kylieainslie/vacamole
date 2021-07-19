# Fit model to OSIRIS data using MCMC

# load packages
library(dplyr)
library(ggplot2)
library(lazymcmc)

# read in OSIRIS data -----------------------------------------------
osiris <- readRDS("inst/extdata/data/Osiris_Data_20210715_1054.rds")

osiris1 <- osiris %>%
  filter(!is.na(date)) %>%
  mutate(roll_avg = zoo::rollmean(inc, k = 7, fill = 0))

osiris2 <- osiris1 %>%
  filter(date < as.Date("2020-03-16")) # date of first lockdown

# plot data ---------------------------------------------------------
p <- ggplot(osiris1, aes(x = date, y = inc)) +
  geom_line() +
  geom_line(aes(x = date, y = roll_avg, color = "red")) +
  theme(panel.background = element_blank())
p

# specify params for MCMC -------------------------------------------
parTab <- data.frame(names=c("r0", "alpha"),
                     values=c(2.3, 0.1),
                     fixed=c(0,0),
                     steps=c(0.1,0.1),
                     lower_bound=c(1, 0),
                     upper_bound=c(10, 1)
                     )

mcmcPars <- c("iterations"=3000,"popt"=0.44,"opt_freq"=100,
              "thin"=1,"adaptive_period"=1000,"save_block"=50)


## Putting model solving code in a function for later use
solve_model <- function(pars, t, params, init){
  
  r0 <- pars["r0"]
  S_diag <- diag(init[c(2:10)])
  rho <- as.numeric(eigs(S_diag %*% params$c_start, 1)$values)
  
  params$beta <- (r0 / rho) * params$gamma
  seir_out <- stochastic_age_struct_seir_ode(times = t,init = init, params = params)
  out <- apply(seir_out, 3, rowSums)
  
  incidence <- c(0,-diff(out[,"S"]))
  daily_cases <- incidence * params$p_report
    #params$sigma * (out[,"E"] + out[,"Ev_1d"] + out[,"Ev_2d"]) * params$p_report
  
  return(daily_cases)
}

## We need to put our likelihood function in a closure environment
## It's important to write the function in this form!
create_lik <- function(parTab, data, PRIOR_FUNC,...){
  par_names <- parTab$names
  
  ## Extract observed incidence
  inc_obs <- data$inc
  
  ## Get times to solve model over
  # tstep <- 1
  # weeks <- floor(dim(data)[1]/7)
  # t <- seq(0,(weeks*7)-1,by=tstep)
  # 
  ## We will pass S0, I0, R0 and SIR_odes 
  #N <- S0 + I0 + R0
  ## using the `...` bit.
  likelihood_func <- function(pars){
    names(pars) <- par_names
    incidence <- solve_model(pars, t, params, init)
    ## Get weekly incidence
    #weekly_incidence <- colSums(matrix(incidence,nrow=7/tstep))
    
    # lik <- sum(dpois(x = inc_obs,lambda = incidence,log=TRUE))
    alpha <- pars["alpha"]
    size <- 1 / alpha
    prob <- size / (incidence + size)
    lik <- sum(dnbinom(x = inc_obs, prob = prob, size = size, log = TRUE))
    
    if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
    lik
    
  }
}

## Prior function
my_prior <- function(pars){
  # diffuse gaussian priors on both parameters
  r0_prior <- dnorm(pars[1],2.3,100,1)
  alpha_prior <- dnorm(pars[2],0.5,100,1)
  return(r0_prior + alpha_prior)
}

## Specify up model inputs
# Create list of parameter values for input into model solver
params <- list(dt = 1/6,
               beta = 0.0004,             # transmission rate
               beta1 = 0.14,              # amplitude of seasonal forcing
               gamma = g,                 # 1/gamma = infectious period
               sigma = s,                 # 1/sigma = latent period
               epsilon = 0.01,            # import case
               N = n_vec,                 # Population (no need to change)
               h = h,                     # Rate from infection to hospital admission/ time from infection to hosp admission
               i1 = i1,
               i2 = i2,
               d = d, 
               d_ic = d_ic,
               d_hic = d_hic,
               r = r,
               r_ic = r_ic,
               p_report = 1/3,
               c_start = t1,
               c_lockdown = t2,
               c_relaxed = t4,
               c_very_relaxed = t3,
               c_normal = t1,
               keep_cm_fixed = TRUE,
               vac_inputs = NULL,
               use_cases = TRUE,                           # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 35.7/100000 * sum(n_vec),      #35.7  # 20 for IC admissions
               no_vac = FALSE,
               t_calendar_start = yday(as.Date("2020-01-01")),   # calendar start date (ex: if model starts on 31 Jan, then t_calendar_start = 31)
               breakpoints = NULL  # breakpoints - start_date    # time points when parameters can change (if NULL, then beta is constant over time)
)

# initial values
empty_state <- c(rep(0, 9))
init <- c(
  t = 0,
  S = n_vec - 1,
  Shold_1d = empty_state,
  Sv_1d = empty_state,
  Shold_2d = empty_state,
  Sv_2d = empty_state,
  E = empty_state,
  Ev_1d = empty_state,
  Ev_2d = empty_state,
  I = c(rep(1,9)),
  Iv_1d = empty_state,
  Iv_2d = empty_state,
  H = empty_state,
  Hv_1d = empty_state,
  Hv_2d = empty_state,
  H_IC = empty_state,
  H_ICv_1d = empty_state,
  H_ICv_2d = empty_state,
  IC = empty_state,
  ICv_1d = empty_state,
  ICv_2d = empty_state,
  D = empty_state,
  R = empty_state,
  Rv_1d = empty_state,
  Rv_2d = empty_state
)

## Starting points, chosen pseudo at random
# starting value for beta based on R0 at beginning of epidemic
startTab <- parTab
startTab$values <- c(3, 0.5)
t <- seq(0,dim(osiris2)[1]-1, by=params$dt)
#osiris2 <- osiris1[t+1,]

output <- run_MCMC(parTab=startTab, data=osiris2, mcmcPars=mcmcPars, filename="fit_to_OSIRIS_no_prior",
                   CREATE_POSTERIOR_FUNC = create_lik, mvrPars = NULL, PRIOR_FUNC=NULL,
                   params = params, init = init)
chain <- read.csv(output$file)
#chain <- read.csv("fit_to_OSIRIS_no_prior_univariate_chain.csv")
plot(coda::as.mcmc(chain[,c("r0", "alpha")]))

chain1 <- chain[chain$sampno >= mcmcPars["adaptive_period"], ] #  ,2:(ncol(chain)-1)

## Use the previous chain to get a good starting
## covariance matrix
startTab$values <- get_best_pars(chain)
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)
# 
## Remove restrictions on parameter bounds -
## not important with the multivariate sampler
startTab$lower_bound <- -Inf
startTab$upper_bound <- Inf
# 
## Run with mvrPars
mcmcPars <- c("iterations"=20000,"popt"=0.234,"opt_freq"=500,
              "thin"=1,"adaptive_period"=7500,"save_block"=100)
# 
output2 <- run_MCMC(parTab=startTab, data=dat, mcmcPars=mcmcPars, filename="SIR_fitting",
                    CREATE_POSTERIOR_FUNC = create_lik, mvrPars = mvrPars, PRIOR_FUNC=my_prior,
                    S0=999,I0=1,R0=0,SIR_odes=SIR_odes)
chain2 <- read.csv(output2$file)
chain2 <- chain2[chain2$sampno >= mcmcPars["adaptive_period"],]
plot(coda::as.mcmc(chain2[,c("R0","gamma")]))

## Get the maximum likelihood parameters and estimate
## the model trajectory for each day
best_pars <- get_best_pars(chain1)
best_trajectory <- solve_model(best_pars, t, params, age_struct_seir_ode)
#best_incidence <- colSums(matrix(best_trajectory,nrow=7))
best_dat <- data.frame(t = t, y = best_trajectory)

## Bit of code to generate prediction intervals
n_samps <- 100
trajectories <- matrix(nrow = n_samps,ncol = length(t))
samps <- sample(nrow(chain1),n_samps)
for(i in 1:n_samps){
  print(i)
  pars <- as.numeric(chain1[samps[i],c("beta")])
  names(pars) <- c("beta")
  trajectories[i,] <- solve_model(pars, t, params, age_struct_seir_ode)
}
bounds <- apply(trajectories,2,function(x) quantile(x, c(0.025,0.975)))

plot(osiris1$inc ~ t, col = "red", pch = 16, 
     xlab = "Time (days)",ylab = "Daily Cases", ylim = c(0,7000)) 
lines(best_dat,col = "blue")
lines(bounds[1,], col = "blue", lty = 2, lwd = 0.5)
lines(bounds[2,], col = "blue", lty = 2, lwd = 0.5)
legend("topright", c("Data","Model","95% credible intervals"),
       col = c("red","blue","blue"), lty = c(0,1,2), pch = c(16,NA,NA))

