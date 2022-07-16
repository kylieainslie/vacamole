# define VE waning function --------------------------------------------------
# it is based on the weighted overage of the VE, daily vaccination rate, and
# time since vaccination

calc_waning <- function(prop, time_points, k, t0){

  tmp <- rep(NA, length(time_points))

  # loop over time points to calculate weighted average of waning based on 
  # proportion of people vaccinated within time point 0 and the current time point
  for(t in 1:length(time_points)){
    if(is.na(time_points[t])){
      next  # go to next time point if t = NaN (because vaccination hasn't started yet)
    } 
    # calculate waning function for time_points
    t_vec <- 1:t
    wane_func_vals = (1 / (1 + exp(-k * (t_vec - t0))))
    # recalculate proportion of people vaccinated at each time point
    prop_recalc <- prop[t:1]/sum(prop[t:1])
    # calculate waning
    tmp[t] <- sum(prop_recalc * wane_func_vals)
  }
  
  # for use in group_modify(), output needs to be a dataframe
  data.frame(t = time_points, w = tmp)
}
