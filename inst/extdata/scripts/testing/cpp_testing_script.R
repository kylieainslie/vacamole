# script for testing C++ functions

# source funtions from C++ file
sourceCpp("src/assign_hh.cpp")

hh_info_dat <- data.frame(hh_size = c(1,2,2,3,3,4,4,5,10), 
                          hh_type = c("Single", 
                                      "Couple","Single, 1 Child",
                                      "Couple, 1 Child", "Single, 2 Children",
                                      "Couple, 2 Children", "Single, 3 Children",
                                      "Couple, 3 Children","Multi-person (other)"),
                          hh_prob = c(0.385, 0.283, 0.045, 0.095, 0.022, 0.114, 0.007, 0.045, 0.005))

# simulate a bunch of populations to check averages
sims <- 100
pop_avg <- 0
hh_prop_avg <- c(rep(0, 6))

for(i in 1:sims){
  # simulate population
  my_households <- assign_hh(pop_size = 1000, hh_info_mat = hh_info_dat)

  n_hh <- dim(my_households)[1]
  
  # check proportion of households of each size
  if(length(table(my_households$hh_size)) == 5){
    tab_hh_sizes <- c(table(my_households$hh_size), 0)
  } else{tab_hh_sizes <- table(my_households$hh_size)}
  
  hh_prop_avg <- hh_prop_avg + (tab_hh_sizes / n_hh)

  # check the number of people
  pop_avg <- pop_avg + sum(my_households$hh_size)

}

pop_avg/sims
hh_prop_avg/sims

# source funtions from C++ file
sourceCpp("src/initialize_pop.cpp")

# simulate populations with households and check age distributions
n <- 1000
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 0.14514332, 0.12092904,
              0.08807406, 0.04622194)
sims <- 1000
pop_avg <- 0
age_group_prop_avg <- c(rep(0, 9))

for(i in 1:sims){
  # simulate population
  my_pops <- initialize_pop(pop_size = n, age_dist = age_dist, hh_info = hh_info_dat, n_postcodes = 10)
  
  n_indiv <- dim(my_pops)[1]
  
  # check proportion of households of each size
  # if(length(table(my_households$hh_size)) == 5){
  #   tab_hh_sizes <- c(table(my_households$hh_size), 0)
  # } else{tab_hh_sizes <- table(my_households$hh_size)}
  # 
  tab_age_group_prop <- table(my_pops[,"Age_Group"]) / n_indiv
  age_group_prop_avg <- age_group_prop_avg + tab_age_group_prop
  
  # check the number of people
  pop_avg <- pop_avg + n_indiv
  
}

pop_avg/sims
age_group_prop_avg/sims
age_dist
