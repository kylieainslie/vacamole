# Some wrangling of contact matrices
library(plyr)
library(dplyr)
library(tidyr)

# read in all contact matrices ------------------------------------
contact_matrices_pienter3 <- readRDS("inst/extdata/data/contact_matrices/Contactpatterns_Pienter3_10y.rds") # september 2020
contact_matrices_pico1 <- readRDS("inst/extdata/data/contact_matrices/Contactpatterns_PICO1_10y.rds") # september 2020
contact_matrices_pico2 <- readRDS("inst/extdata/data/contact_matrices/Contactpatterns_PICO2_10y.rds") # september 2020
contact_matrices_pico3 <- readRDS("inst/extdata/data/contact_matrices/Contactpatterns_PICO3_10y.rds") # september 2020
contact_matrices_pico4 <- readRDS("inst/extdata/data/contact_matrices/Contactpatterns_PICO4_10y.rds") # february 2021
contact_matrices_pico5 <- readRDS("inst/extdata/data/contact_matrices/Contactpatterns_PICO5_10y.rds") # june 2021

# need to run this before running convert_contact_matrices()------
age_dist <- c(
  0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463,
  0.14514332, 0.12092904, 0.08807406, 0.04622194
)
n <- 17407585 # Dutch population size
n_vec <- n * age_dist
N_diag <- diag(1 / n_vec)

rel_trans <- c(1.000, 3.051, 5.751, 3.538, 3.705, 4.365, 5.688, 5.324, 7.211)
get_transmission_matrix <- function(x, contact_mat) {
  # multiply by relative susc/inf
  tmp <- sweep(contact_mat, 1, x, "*") # rows
  rtn <- sweep(tmp, 2, x, "*") # columns
  
  # output
  return(rtn)
}
# -----------------------------------------------------------------

# Define function to convert contact matrices --------------------
convert_contact_matrices <- function(x){
  # create empty list to store converted matrices
  rtn <- list()
  # loop over different realisations of contact matrices
  var_names <- paste0("c_smt.", 1:200)
  
  for (i in 1:200){
    
    tmp <- x %>%
      select(part_age, cnt_age, var_names[i]) %>%
      mutate(contact_type = c(rep("all", 81),
                              rep("community", 81),
                              rep("household", 81))) %>%
      filter(contact_type == "all") %>%
      select(-contact_type) %>%
      pivot_wider(., names_from = cnt_age, values_from = var_names[i]) %>%
      select(-part_age)
    
    # convert to transmission matrix
    tmp1 <- as.matrix(tmp) %*% N_diag
    tmp2 <- get_transmission_matrix(rel_trans, tmp1)
    
    rtn[[i]] <- tmp2
    
  }
  
  # add mean to list of matrices
  tmp_mean <- aaply(laply(rtn, as.matrix), c(2, 3), mean)
  #tmp1_mean <- as.matrix(tmp_mean) %*% N_diag
  #tmp2_mean <- get_transmission_matrix(rel_trans, tmp1_mean)
  rtn$mean <- tmp_mean
  
  return(rtn)
}
# -----------------------------------------------------------------

# convert contact matrices and save to file
baseline_2017 <- convert_contact_matrices(contact_matrices_pienter3)
saveRDS(baseline_2017, file = "inst/extdata/data/contact_matrices/contact_matrices_baseline_2017.rds")

april_2020 <- convert_contact_matrices(contact_matrices_pico1)
saveRDS(april_2020, file = "inst/extdata/data/contact_matrices/contact_matrices_april_2020.rds")

june_2020 <- convert_contact_matrices(contact_matrices_pico2)
saveRDS(june_2020, file = "inst/extdata/data/contact_matrices/contact_matrices_june_2020.rds")

september_2020 <- convert_contact_matrices(contact_matrices_pico3)
saveRDS(september_2020, file = "inst/extdata/data/contact_matrices/contact_matrices_september_2020.rds")

february_2021 <- convert_contact_matrices(contact_matrices_pico4)
saveRDS(february_2021, file = "inst/extdata/data/contact_matrices/contact_matrices_february_2021.rds")

june_2021 <- convert_contact_matrices(contact_matrices_pico5)
saveRDS(june_2021, file = "inst/extdata/data/contact_matrices/contact_matrices_june_2021.rds")
# -----------------------------------------------------------------

