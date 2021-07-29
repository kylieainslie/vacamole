# Some wrangling of contact matrices
library(plyr)
library(dplyr)
library(tidyr)
# read in all contact matrices
#contact_matrices_pico3 <- read.delim("inst/extdata/data/S2_contact_matrices_withPico3_10y.tsv", header=TRUE, allowEscapes=FALSE, sep="\t")
contact_matrices_pico4 <- readRDS("inst/extdata/data/Contactpatterns_PICO4_10y.rds") # february 2021
contact_matrices_pico5 <- readRDS("inst/extdata/data/Contactpatterns_PICO5_10y.rds") # june 2021

var_names <- paste0("c_smt.", 1:200)

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

# February 2021
february_2021 <- list()

for (i in 1:200){
  
  tmp <- contact_matrices_pico4 %>%
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
  
  february_2021[[i]] <- tmp

}

# add mean to list of matrices
tmp_mean <- aaply(laply(february_2021, as.matrix), c(2, 3), mean)
tmp1_mean <- as.matrix(tmp_mean) %*% N_diag
tmp2_mean <- get_transmission_matrix(rel_trans, tmp1_mean)
february_2021$mean <- tmp2_mean
 
# save contact matrices to file
saveRDS(february_2021, file = "inst/extdata/data/contact_matrices_february_2021.rds")

# June 2021
june_2021 <- list()

for (i in 1:200){
  
  tmp <- contact_matrices_pico5 %>%
    select(part_age, cnt_age, var_names[i]) %>%
    mutate(contact_type = c(rep("all", 81),
                            rep("community", 81),
                            rep("household", 81))) %>%
    filter(contact_type == "all") %>%
    select(-contact_type) %>%
    pivot_wider(., names_from = cnt_age, values_from = var_names[i]) %>%
    select(-part_age)
  
  june_2021[[i]] <- tmp
  
}

# add mean to list of matrices
june_2021$mean <- aaply(laply(june_2021, as.matrix), c(2, 3), mean)

# save contact matrices to file
saveRDS(june_2021, file = "inst/extdata/data/contact_matrices_june_2021.rds")


### old code ### -----------------------------------------------
# --------------------------------------------------------------

# demographic info to make contact matrices symmetrical
# population distributions
# demo <- data.frame(age_group = c("[0,10)", "[10,20)", "[20,30)", "[30,40)","[40,50)", "[50,60)",
#                                  "[60,70)", "[70,80)", "[80,Inf]"), 
#                    frac_2017 = c(0.105, 0.118, 0.126, 0.120, 0.138, 0.145, 0.122, 0.0808, 0.0447), 
#                    frac_2019 = c(0.103, 0.116, 0.127, 0.122, 0.131, 0.145, 0.121, 0.0881, 0.0462),
#                    frac_2020 = c(0.102, 0.115, 0.128, 0.123, 0.127, 0.145, 0.121, 0.0904, 0.0472))

# baseline (2017)
# baseline <- contact_matrices_all %>%
#   filter(survey == "baseline") %>%
#   filter(contact_type == "all") %>%
#   mutate(f_pop = rep(demo$frac_2017,9),
#          c_est = round(m_est/f_pop,1)) 
# 
# baseline_asym <- baseline %>%
#   select(-survey, -contact_type, -c_est, -f_pop) %>%
#   pivot_wider(., names_from = cont_age, values_from = m_est)
# 
# baseline_sym <- baseline %>%
#   select(-survey, -contact_type, -m_est, -f_pop) %>%
#   pivot_wider(., names_from = cont_age, values_from = c_est)
# 
# # April 2020 (lockdown)
# april2020 <- contact_matrices_all %>%
#   filter(survey == "April 2020") %>%
#   filter(contact_type == "all") %>%
#   mutate(f_pop = rep(demo$frac_2019,9),
#          c_est = round(m_est/f_pop,1))
# 
# april2020_asym <- april2020 %>%
#   select(-survey, -contact_type, -c_est, -f_pop) %>%
#   pivot_wider(., names_from = cont_age, values_from = m_est)
# 
# april2020_sym <- april2020 %>%
#   select(-survey, -contact_type, -m_est, -f_pop) %>%
#   pivot_wider(., names_from = cont_age, values_from = c_est)
# 
# # June 2020 
# june2020 <- contact_matrices_all %>%
#   filter(survey == "June 2020") %>%
#   filter(contact_type == "all") %>%
#   mutate(f_pop = rep(demo$frac_2019,9),
#          c_est = round(m_est/f_pop,1))
# 
# june2020_asym <- june2020 %>%
#   select(-survey, -contact_type, -c_est, -f_pop) %>%
#   pivot_wider(., names_from = cont_age, values_from = m_est)
# 
# june2020_sym <- june2020 %>%
#   select(-survey, -contact_type, -m_est, -f_pop) %>%
#   pivot_wider(., names_from = cont_age, values_from = c_est)
# 
# # September 2020 
# september2020 <- contact_matrices_all %>%
#   filter(survey == "September 2020") %>%
#   filter(contact_type == "all") %>%
#   mutate(f_pop = rep(demo$frac_2019,9),
#          c_est = round(m_est/f_pop,1)) 
# 
# september2020_asym <- september2020 %>%
#   select(-survey, -contact_type, -c_est, -f_pop) %>%
#   pivot_wider(., names_from = cont_age, values_from = m_est)
# 
# september2020_sym <- september2020 %>%
#   select(-survey, -contact_type, -m_est, -f_pop) %>%
#   pivot_wider(., names_from = cont_age, values_from = c_est)

# put in a list and write to rds
# cm <- list(baseline_sym = baseline_sym,
#            baseline_asym = baseline_asym,
#            april2020_sym = april2020_sym,
#            april2020_asym = april2020_asym,
#            june2020_sym = june2020_sym,
#            june2020_asym = june2020_asym,
#            september2020_sym = september2020_sym,
#            september2020_asym = september2020_asym,
#            february2021_sym = february2021_sym,
#            february2021_asym = february2021_asym
# )
# saveRDS(cm,"inst/extdata/data/contact_matrices_for_model_input.rds")

