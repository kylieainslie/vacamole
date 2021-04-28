# write wrapper function for seir model code that outputs only hospital
# admission counts ------------------------------------------------------
fit_to_data_wrapper_init <- function(x, params, times_vec, s, g){ 
  
  r0 <- params[1]
  init_i <- params[2]
  #s <- params[3]
  #g <- params[4]
  
  tmp <- get_beta(R0 = r0, contact_matrix = c1, N = n_vec, sigma = s, 
                  gamma = s) 
  beta <- tmp$beta
  
  # loop over time periods -------------------------------------------
  params <- list(beta = beta,                    # transmission rate
                 gamma = g,                      # 1/gamma = infectious period
                 sigma = s,                      # 1/sigma = latent period
                 N = n_vec,                      # Population (no need to change)
                 h = p_infection2admission,      # Rate from infection to hospital admission
                 d = p_admission2death,          # Rate from admission to death
                 r = 0.0206,                     # Rate from admission to recovery
                 c_start = c1,
                 no_vac = TRUE
  )
  
  # Specify initial values -------------------------------------------
  times <- times_vec
  #print(length(times_vec))
  timeInt <- times[2]-times[1]             
  init <- c(t = times[1],                  
            S = params$N - c(rep(init_i/9, 9)),
            Shold_1d = empty_state,
            Sv_1d = empty_state,
            Shold_2d = empty_state,
            Sv_2d = empty_state,
            E = empty_state,
            Ev_1d = empty_state,
            Ev_2d = empty_state,
            I = c(rep(init_i/9, 9)),
            Iv_1d = empty_state,
            Iv_2d = empty_state,
            H = empty_state,
            Hv_1d = empty_state,
            Hv_2d = empty_state,
            D = empty_state,
            R = empty_state,
            Rv_1d = empty_state,
            Rv_2d = empty_state
  )                      
  
  # Solve model ------------------------------------------------------
  seir_out <- lsoda(init,times,age_struct_seir_ode,params)
  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output(seir_out)
  hosp_by_age_group <- sweep((out$I + out$Iv_1d + out$Iv_2d), 2, params$h, "*")
  rtn <- rowSums(hosp_by_age_group)
  #print(length(rtn))
  
  return(rtn)
}

# 
fit_to_data_wrapper <- function(x, params, contact_mat, beta, times_vec, 
                                s, g){ 
  
  init_i <- params[1]
  delta <- params[2]
  #print(delta)
  # tmp <- get_beta(R0 = r0, contact_matrix = c1, N = n_vec, sigma = s, 
  #                 gamma = s) 
  # beta <- tmp$beta
  b <- beta * delta
  
  # loop over time periods -------------------------------------------
  params <- list(beta = b,            # transmission rate
                 gamma = g,                      # 1/gamma = infectious period
                 sigma = s,                      # 1/sigma = latent period
                 N = n_vec,                      # Population (no need to change)
                 h = p_infection2admission,      # Rate from infection to hospital admission
                 d = p_admission2death,          # Rate from admission to death
                 r = 0.0206,                     # Rate from admission to recovery
                 c_start = contact_mat,
                 no_vac = TRUE
  )
  
  # Specify initial values -------------------------------------------
  times <- times_vec
  timeInt <- times[2]-times[1]             
  init <- c(t = times[1],                  
            S = params$N - c(rep(init_i/9, 9)),
            Shold_1d = empty_state,
            Sv_1d = empty_state,
            Shold_2d = empty_state,
            Sv_2d = empty_state,
            E = empty_state,
            Ev_1d = empty_state,
            Ev_2d = empty_state,
            I = c(rep(init_i/9, 9)),
            Iv_1d = empty_state,
            Iv_2d = empty_state,
            H = empty_state,
            Hv_1d = empty_state,
            Hv_2d = empty_state,
            D = empty_state,
            R = empty_state,
            Rv_1d = empty_state,
            Rv_2d = empty_state
  )                      
  
  # Solve model ------------------------------------------------------
  seir_out <- lsoda(init,times,age_struct_seir_ode,params)
  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output(seir_out)
  hosp_by_age_group <- sweep((out$I + out$Iv_1d + out$Iv_2d), 2, params$h, "*")
  rtn <- rowSums(hosp_by_age_group)
  
  
  return(rtn)
}
