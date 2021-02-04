# Simulate SEIR compartmental model with 2 vaccine doses

# Load packages ----------------------------------------------------
library(deSolve)
library(reshape2)
library(ggplot2)
library(gridExtra)
source("R/seir_ode.R")

# Specify parameter values -----------------------------------------
# Input parameters:
params <- list(beta = 0.61, 
               gamma = 0.5,                   # R0 = beta/gamma
               sigma = 0.5,                   # 1/sigma = latent period
               N = 100000,                   # Population (no need to change)
               vac_per_day = 25000,           # Number of vaccines per day (dose 1)
               vac_per_day2 = 25000,          # Number of vaccines per day (dose 2)
               tv = 18,                       # Time vaccination starts (dose 1)
               tv2 = 50,                      # Time vaccination starts (dose 2)
               delay = 14,                    # Delay from vaccination to protection (days)
               delay2 = 14,                   # Delay for dose 2
               eta = 1 - 0.7,                 # 1 - VE (dose 1)
               eta2 = 1 - 0.9,                # 1 - VE (dose 2)
               uptake = 0.85,                 # Proportion of population able and willing to be vaccinated
               h = 0.0251,                    # Rate from infection to hospital admission
               d = 0.106,                     # Rate from admission to death
               r = 0.0206,                    # Rate from admission to recovery
               constant_foi = FALSE,
               init_inf = 0,
               vac_input_perc = TRUE          # is vaccine distribution a percentage of the population?
)                  

# Specify initial values -------------------------------------------
times <- seq(0,200,length.out = 201)       # Vector of times
timeInt <- times[2]-times[1]             # Time interval (for technical reasons)
init <- c(t = times[1],                  # Initial conditions
          S = params$N,
          Shold = 0,
          Sv = 0,
          Shold2 = 0,
          Sv2 = 0,
          E = 0,
          Ev = 0,
          Ev2 = 0,
          I = 20,
          Iv = 0,
          Iv2 = 0,
          H = 0,
          Hv = 0,
          Hv2 = 0,
          D = 0,
          R = 0,
          Rv = 0,
          Rv2 = 0
)                      

# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times,seir_ode,params)
seir_out <- as.data.frame(seir_out)
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

# Create object for plotting:
df <- data.frame(time = time, 
                 incidence = inc, 
                 hosp_admissions = hosp)
#saveRDS(df, file = paste0("inst/extdata/results/",output_tag,"_output.rds"))

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
