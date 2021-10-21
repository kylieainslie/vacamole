# -----------------------------------------------------------------
# Convert contact matrices for model input
#------------------------------------------------------------------
library(plyr)
library(dplyr)
library(tidyr)

# read in all contact matrices ------------------------------------
contact_matrices_pienter3 <- readRDS("inst/extdata/data/contact_matrices/raw/Contactpatterns_Pienter3_10y.rds") # April 2017
contact_matrices_pico1 <- readRDS("inst/extdata/data/contact_matrices/raw/Contactpatterns_PICO1_10y.rds")       # April 2020
contact_matrices_pico2 <- readRDS("inst/extdata/data/contact_matrices/raw/Contactpatterns_PICO2_10y.rds")       # June 2020
contact_matrices_pico3 <- readRDS("inst/extdata/data/contact_matrices/raw/Contactpatterns_PICO3_10y.rds")       # September 2020
contact_matrices_pico4 <- readRDS("inst/extdata/data/contact_matrices/raw/Contactpatterns_PICO4_10y.rds")       # February 2021
contact_matrices_pico5 <- readRDS("inst/extdata/data/contact_matrices/raw/Contactpatterns_PICO5_10y.rds")       # June 2021

# need to run this before running convert_contact_matrices()------
age_dist <- c(
  0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463,
  0.14514332, 0.12092904, 0.08807406, 0.04622194
)
n <- 17407585 # Dutch population size
n_vec <- n * age_dist
N_diag <- diag(1 / n_vec)

# relative susceptibility/infectiousness by age group
# the 0-9 year age group is the reference group
rel_trans <- c(1.000, 3.051, 5.751, 3.538, 3.705, 4.365, 5.688, 5.324, 7.211)

# convert contact matrices ----------------------------------------
baseline_2017 <- convert_contact_matrices(contact_matrices_pienter3)
april_2020 <- convert_contact_matrices(contact_matrices_pico1)
june_2020 <- convert_contact_matrices(contact_matrices_pico2)
september_2020 <- convert_contact_matrices(contact_matrices_pico3)
february_2021 <- convert_contact_matrices(contact_matrices_pico4)
june_2021 <- convert_contact_matrices(contact_matrices_pico5)

# save output -----------------------------------------------------
sve_path <- "inst/extdata/data/contact_matrices/converted/"
saveRDS(baseline_2017, file = paste0(save_path,"contact_matrices_baseline_2017.rds"))
saveRDS(april_2020, file = paste0(save_path,"contact_matrices_april_2020.rds"))
saveRDS(june_2020, file = paste0(save_path,"contact_matrices/contact_matrices_june_2020.rds"))
saveRDS(september_2020, file = paste0(save_path,"contact_matrices_september_2020.rds"))
saveRDS(february_2021, file = paste0(save_path,"contact_matrices_february_2021.rds"))
saveRDS(june_2021, file = paste0(save_path,"contact_matrices_june_2021.rds"))

