# Simulate SEIR compartmental model with 2 vaccine doses

# Load packages ----------------------------------------------------
library(deSolve)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(readxl)
source("R/seir_ode.R")

# Read in vac schedule ---------------------------------------------
vac_schedule_p <- read_xlsx("inst/extdata/data/cum_upt_A.xlsx", sheet = 1)
vac_schedule_m <- read_xlsx("inst/extdata/data/cum_upt_A.xlsx", sheet = 2)
vac_schedule_az <- read_xlsx("inst/extdata/data/cum_upt_A.xlsx", sheet = 3)

# Specify parameter values -----------------------------------------
ve_est <- list(Pfizer  = c(0.926, 0.948),
               Moderna = c(0.896, 0.941),
               AZ_60 = c(0.583, 0.621),
               AZ_30 = c(0.3, 0.3),
               AZ_10 = c(0.1, 0.1))

days_to_protection <- list(Pfizer  = c(14, 7),
                           Moderna = c(14, 14),
                           AstraZeneca = c(21, 14))

vac_start_day <- list(Pfizer = c(52, 94),
                      AstraZeneca = c(18, 92),
                      Pfizer_2 = c(88, 130),
                      AstraZeneca_2 = c(25, 109))

initial_inputs <- list(y60_64 = list(N = 1117798,
                                     E = 1902 * .5288,
                                     I = 1941 * 0.5288,
                                     R = 250590 * 0.5288),
                       y65_69 = list(N = 996048,
                                     E = 1902 * .4712,
                                     I = 1941 * 0.4712,
                                     R = 250590 * 0.4712),
                       all = list(N = 2113846,
                                  E = 1902,
                                  I = 1941,
                                  R = 250590))

output_tag <- "new_AZ_30_cfoi_60_64"
inputs <- initial_inputs$y60_64
# dose1_input <- as.numeric(vac_schedule_az$az_d1_7[which(as.Date(vac_schedule_az$date) > as.Date("2021-01-20"))])
# dose2_input <- as.numeric(vac_schedule_az$az_d2_7[which(as.Date(vac_schedule_az$date) > as.Date("2021-01-20"))])
# Input parameters:
params <- list(beta = 0.61, 
            gamma = 0.5,                   # R0 = beta/gamma
            sigma = 0.5,                   # 1/sigma = latent period
            N = inputs$N,                  # Population (no need to change)
            vac_per_day = 25000,     # Number of vaccines per day (dose 1)
            vac_per_day2 = 25000,    # Number of vaccines per day (dose 2)
            tv = vac_start_day$AstraZeneca_2[1],                     # Time vaccination starts (dose 1)
            tv2 = vac_start_day$AstraZeneca_2[2],                    # Time vaccination starts (dose 2)
            delay = days_to_protection$AstraZeneca[1],   # Delay from vaccination to protection (days)
            delay2 = days_to_protection$AstraZeneca[2],  # Delay for dose 2
            eta = 1- ve_est$AZ_30[1],        # 1 - VE (dose 1)
            eta2 = 1- ve_est$AZ_30[2],       # 1 - VE (dose 2)
            uptake = 0.85,                 # Proportion of population able and willing to be vaccinated
            h = 0.0251,                    # Rate from infection to hospital admission
            d = 0.106,                     # Rate from admission to death
            r = 0.0206,                    # Rate from admission to recovery
            constant_foi = TRUE,
            init_inf = inputs$I,
            vac_input_perc = TRUE          # is vaccine distribution a percentage of the population?
)                  

# Specify initial values -------------------------------------------
times <- seq(0,200,length.out = 201)       # Vector of times
timeInt <- times[2]-times[1]             # Time interval (for technical reasons)
init <- c(t = times[1],                  # Initial conditions
          S = inputs$N-inputs$E-inputs$I-inputs$R,
          Shold = 0,
          Sv = 0,
          Shold2 = 0,
          Sv2 = 0,
          E = inputs$E,
          Ev = 0,
          Ev2 = 0,
          I = inputs$I,
          Iv = 0,
          Iv2 = 0,
          H = 50,
          D = 0,
          R = inputs$R,
          Rv = 0,
          Rv2 = 0)                      

# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times,seir_ode,params)

# Summarise results ------------------------------------------------
beta <- params$beta * timeInt
eta <- params$eta
eta2 <- params$eta2
N <- params$N
h <- params$h
gamma <- params$gamma


lambda <- beta * rowSums(seir_out[,c(11:13)])/N
time <- seir_out[,1]
S <- seir_out[,3]
Shold <- seir_out[,4]
Sv <- seir_out[,5]
Shold2 <- seir_out[,6]
Sv2 <- seir_out[,7]
inc <- (S + Shold + (eta * (Sv + Shold2)) + (eta2 * Sv2)) * lambda
I <- seir_out[,11]
Iv <- seir_out[,12]
Iv2 <- seir_out[,13]
hosp <- h * (I + Iv + Iv2)
hosp2 <- inc * params$h

# Create object for plotting:
df <- data.frame(time = time, 
                 incidence = inc, 
                 hosp_admissions = hosp,
                 hosp_admissions_from_inc = hosp2)
saveRDS(df, file = paste0("inst/extdata/results/",output_tag,"_output.rds"))

# Calculate summary data:
Value <- c(time[which.max(inc)],max(inc),max(hosp),sum(inc),sum(hosp))
Value <- round(Value, digits=2)
Value <- data.frame(Value)
rownames(Value) <- c("Peak time","Peak inc.","Peak hosp.","Cum. inc.","Cum. hosp.")
print(Value)

# Make plot:
colors <- c("Incidence" = "red", "Hospital Admissions" = "Green")

g <- ggplot(df, aes(time)) +
  geom_line(aes(y=incidence, color="Incidence"), lwd=1) +
  geom_line(aes(y=hosp_admissions, color="Hospital Admissions"),lwd=1) +
  geom_vline(xintercept = params[[7]],linetype="dashed", color = "blue") +
  geom_vline(xintercept = params[[8]],linetype="dotted", color = "blue") +
  labs(y = "", x = "Time (days)", color = "Legend") +
  scale_color_manual(values = colors) +
  #annotation_custom(tableGrob(Value), xmin=125, xmax=200, ymin=25, ymax=100) +
  theme(legend.position = "bottom",
        panel.background = element_blank()
  )
plot(g)

