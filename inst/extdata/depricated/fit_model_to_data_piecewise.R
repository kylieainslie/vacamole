# Fit model to data: piecewise -------------------------------------

# Load packages ----------------------------------------------------
library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
library(rARPACK)
library(readr)
library(lubridate)
library(FSA)

# Source functions ------------------------------------------------------
source("R/age_struct_seir_ode.R")
source("R/postprocess_age_struct_model_output.R")
source("R/get_foi.R")
source("R/get_beta.R")
source("R/get_vac_rate.R")
source("R/get_vac_rate_2.R")

# read in hospital admission counts -------------------------------------
nice <- readRDS("inst/extdata/data/nice_admissions_20210212.rds") %>%
  filter(AdmissionDate >= "2020-01-01")
# make dates consecutive
nice1 <- complete(nice, AdmissionDate = full_seq(AdmissionDate, 1)) %>%
  mutate(n = ifelse(is.na(n), 0, n),
         time = seq(1,dim(nice1)[1], by = 1))

# list of dates at different intervention points ------------------------
t_vec <- lubridate::yday(as.Date(c("2020-01-01", "2020-03-28", 
                                   "2020-06-01", "2020-09-29", 
                                   "2020-12-15", "2021-02-11"))) 
t_vec[length(t_vec)] <- t_vec[length(t_vec)] + 366 #2020 was a leap year

# make separate data sets for each time interval ------------------------
data_sets <- list()
for (i in 1:(length(t_vec)-1)){
  if (i == length(t_vec) - 1){
    times <- seq(t_vec[i], t_vec[i+1], by = 1) 
  } else {
    times <- seq(t_vec[i], t_vec[i+1] - 1, by = 1) 
  }
  df <- nice1 %>%
    filter(time %in% times)
  data_sets[[i]] <- df
}

# probabilities ---------------------------------------------------------
p_infection2admission <- c(0.003470, 0.000377, 0.000949, 0.003880, 0.008420, 
                           0.016500,0.025100, 0.049400, 0.046300)
p_admission2death <- c(0.00191, 0.00433, 0.00976, 0.02190, 0.02500, 0.04010,
                       0.10600, 0.22900, 0.31100)

# contact matrix --------------------------------------------------------
contact_matrices_all <- readRDS("inst/extdata/data/contact_matrices_for_model_input.rds")
c1 <- as.matrix(contact_matrices_all$baseline[,-1])
c2 <- as.matrix(contact_matrices_all$april2020[,-1])
c3 <- as.matrix(contact_matrices_all$june2020[,-1])
c4 <- as.matrix(contact_matrices_all$september2020[,-1])
# age distribution and pop size -----------------------------------------
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.04622194)
n <- 17407585                                 # Dutch population size
n_vec <- n * age_dist

empty_state <- c(rep(0,9))


# use nls to fit non-linear least squares ------------------------------
# ----------------------------------------------------------------------
# fit initial epidemic timepoint (before any measures are introduced)
model_fit <- nls(n ~ fit_to_data_wrapper_init(AdmissionDate, 
                                              params, 
                                              times_vec = time,
                                              s = 0.5,
                                              g = 0.5),
                 data = data_sets[[1]], # where x and y are
                 start = list(params = c(3, 0.03)), #initial conditions 
                 lower = c(2, 0), # lower bound of parameter values 
                 upper = c(4, 10), # upper bound of parameter values 
                 trace = T,
                 algorithm = "port", # must be used to specify bounds
                 control = list(maxiter = 1000000, 
                                warnOnly=T)
)
summary(model_fit)
# plot -----------------------------------------------------------------
plot(n~AdmissionDate,data=data_sets[[1]],type = "l")
lines(x = data_sets[[1]]$AdmissionDate, y = model_fit$m$fitted(), 
      type = "l", col = "red")

# ----------------------------------------------------------------------
# fit other time intervals using values from initial fit
tmp <- get_beta(R0 = coef(model_fit)[1], contact_matrix = c1, N = n_vec, sigma = 0.5, 
                gamma = 0.5) 
beta <- tmp$beta
i = 5
fit <- nls(n ~ fit_to_data_wrapper(AdmissionDate, 
                                     params, 
                                     contact_mat = c1,
                                     beta = beta,
                                     times_vec = time,
                                     s = 0.5,
                                     g = 0.5
                                     ),
             data = data_sets[[5]], # where x and y are
             start = list(params = c(18742, 0.4)), #initial conditions 
             lower = c(0, 0), # lower bound of parameter values 
             upper = c(10e6, 2), # upper bound of parameter values 
             trace = T,
             algorithm = "port", # must be used to specify bounds
             control = list(maxiter = 1000000, warnOnly=T)
            )

# smooth data
loess_fit <- loess(n ~ time, data=data_sets[[5]], span=0.15)
smoothed <- predict(loess_fit) 
plot(n~AdmissionDate,data=data_sets[[5]],type = "l")
lines(smoothed, x = data_sets[[5]]$AdmissionDate, col = "blue")
lines(x = data_sets[[5]]$AdmissionDate, y = fit$m$fitted(), 
      type = "l", col = "red")

