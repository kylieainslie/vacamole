#' Stochastic Age-structured SEIR ODE model of vaccination with 2 doses and delay to protection
#' @param times vector of times
#' @param init list of inititial states
#' @param params list of parameter values
#' @return List of summary results
#' @keywords vacamole
#' @export
# Define model -----------------------------------------------------
stochastic_age_struct_seir_ode <- function(times,init,params){
  with(as.list(c(params,init)), {
    #print(t)
    # define initial state vectors from input ----------------------
    S = c(S1, S2, S3, S4, S5, S6, S7, S8, S9)
    Shold_1d = c(Shold_1d1, Shold_1d2, Shold_1d3, Shold_1d4, Shold_1d5, Shold_1d6, 
                 Shold_1d7, Shold_1d8, Shold_1d9)
    Sv_1d = c(Sv_1d1, Sv_1d2, Sv_1d3, Sv_1d4, Sv_1d5, Sv_1d6, Sv_1d7, Sv_1d8, Sv_1d9)
    Shold_2d = c(Shold_2d1, Shold_2d2, Shold_2d3, Shold_2d4, Shold_2d5, Shold_2d6, 
                 Shold_2d7, Shold_2d8, Shold_2d9)
    Sv_2d = c(Sv_2d1, Sv_2d2, Sv_2d3, Sv_2d4, Sv_2d5, Sv_2d6, Sv_2d7, Sv_2d8, Sv_2d9)
    E = c(E1, E2, E3, E4, E5, E6, E7, E8, E9)
    Ev_1d = c(Ev_1d1, Ev_1d2, Ev_1d3, Ev_1d4, Ev_1d5, Ev_1d6, Ev_1d7, Ev_1d8, Ev_1d9)
    Ev_2d = c(Ev_2d1, Ev_2d2, Ev_2d3, Ev_2d4, Ev_2d5, Ev_2d6, Ev_2d7, Ev_2d8, Ev_2d9)
    I = c(I1, I2, I3, I4, I5, I6, I7, I8, I9)
    Iv_1d = c(Iv_1d1, Iv_1d2, Iv_1d3, Iv_1d4, Iv_1d5, Iv_1d6, Iv_1d7, Iv_1d8, Iv_1d9)
    Iv_2d = c(Iv_2d1, Iv_2d2, Iv_2d3, Iv_2d4, Iv_2d5, Iv_2d6, Iv_2d7, Iv_2d8, Iv_2d9)
    H = c(H1, H2, H3, H4, H5, H6, H7, H8, H9)
    Hv_1d = c(Hv_1d1, Hv_1d2, Hv_1d3, Hv_1d4, Hv_1d5, Hv_1d6, Hv_1d7, Hv_1d8, Hv_1d9)
    Hv_2d = c(Hv_2d1, Hv_2d2, Hv_2d3, Hv_2d4, Hv_2d5, Hv_2d6, Hv_2d7, Hv_2d8, Hv_2d9)
    H_IC = c(H_IC1, H_IC2, H_IC3, H_IC4, H_IC5, H_IC6, H_IC7, H_IC8, H_IC9)
    H_ICv_1d = c(H_ICv_1d1, H_ICv_1d2, H_ICv_1d3, H_ICv_1d4, H_ICv_1d5, H_ICv_1d6, H_ICv_1d7, H_ICv_1d8, H_ICv_1d9)
    H_ICv_2d = c(H_ICv_2d1, H_ICv_2d2, H_ICv_2d3, H_ICv_2d4, H_ICv_2d5, H_ICv_2d6, H_ICv_2d7, H_ICv_2d8, H_ICv_2d9)
    IC = c(IC1, IC2, IC3, IC4, IC5, IC6, IC7, IC8, IC9)
    ICv_1d = c(ICv_1d1, ICv_1d2, ICv_1d3, ICv_1d4, ICv_1d5, ICv_1d6, ICv_1d7, ICv_1d8, ICv_1d9)
    ICv_2d = c(ICv_2d1, ICv_2d2, ICv_2d3, ICv_2d4, ICv_2d5, ICv_2d6, ICv_2d7, ICv_2d8, ICv_2d9)
    D = c(D1, D2, D3, D4, D5, D6, D7, D8, D9)
    R = c(R1, R2, R3, R4, R5, R6, R7, R8, R9)
    Rv_1d = c(Rv_1d1, Rv_1d2, Rv_1d3, Rv_1d4, Rv_1d5, Rv_1d6, Rv_1d7, Rv_1d8, Rv_1d9)
    Rv_2d = c(Rv_2d1, Rv_2d2, Rv_2d3, Rv_2d4, Rv_2d5, Rv_2d6, Rv_2d7, Rv_2d8, Rv_2d9)
    
    # determine vaccination rate -----------------------------------
    time_point <- floor(times) + 1
    
    alpha <- vac_inputs$alpha_dose1[time_point,]
    alpha2 <- vac_inputs$alpha_dose2[time_point,]
    eta <- vac_inputs$eta_dose1[time_point,]
    eta2 <- vac_inputs$eta_dose2[time_point,]
    delay <- vac_inputs$delay_dose1[time_point,]
    delay2 <- vac_inputs$delay_dose2[time_point,]
    eta_hosp <- vac_inputs$eta_hosp_dose1[time_point,]
    eta_hosp2 <- vac_inputs$eta_hosp_dose2[time_point,]
    eta_trans <- vac_inputs$eta_trans_dose1[time_point,]
    eta_trans2 <- vac_inputs$eta_trans_dose2[time_point,]
    
    # determine contact matrix based on criteria --------------------
    ic_admin <- sum(i1 * (H + Hv_1d + Hv_2d))
    
    cases <- sum(sigma * (E + Ev_1d + Ev_2d) * p_report)
    criteria <- (use_cases) * cases + (!use_cases) * ic_admin 
    # initialise flags
    if(times == 0){
      flag_relaxed <- 0
      flag_very_relaxed <- 0
      flag_normal <- 0
    }
    
    # determine contact matrix to use based on criteria
    tmp2 <- choose_contact_matrix(params, times, criteria, flag_relaxed, 
                                  flag_very_relaxed, flag_normal, keep_fixed = keep_cm_fixed)
    contact_mat <- tmp2$contact_matrix
    flag_relaxed <- tmp2$flag_relaxed
    flag_very_relaxed <- tmp2$flag_very_relaxed
    flag_normal <- tmp2$flag_normal
    
    # determine force of infection ----------------------------------
    calendar_day <- t_calendar_start + times
    beta_t <- beta * (1 + beta1 * cos(2 * pi * calendar_day/365.24))
    # lambda <- beta * (contact_mat %*% (I + (eta_trans * Iv_1d) + (eta_trans2 * Iv_2d)))
    lambda <- beta_t * (contact_mat %*% (I + (eta_trans * Iv_1d) + (eta_trans2 * Iv_2d)))
    # ---------------------------------------------------------------
    ### probabilities of transitioning
    # from S
    p_S_ <- as.numeric(1 - exp(-lambda-alpha))            # total probability of moving from S
    p_S_Shold1 <- as.numeric(alpha/(lambda + alpha))      # relative probability of moving from S -> Shold_1d
    p_S_E <- as.numeric(lambda/(lambda + alpha))          # relative probability of moving from S -> E
    # from Shold_1d
    p_Shold1_ <- as.numeric(1 - exp(-lambda-(1/delay)))               # total probability of moving from Shold_1d
    p_Shold1_E <- as.numeric(lambda/(lambda + (1/delay)))             # relative probability of moving from Shold_1d -> E
    p_Shold1_Sv1 <- as.numeric((1/delay)/(lambda + (1/delay)))        # relative probability of moving from Shold_1d -> Sv1
    # from Sv_1d
    p_Sv1_ <- as.numeric(1 - exp(-eta*lambda-alpha2))                 # total probability of moving from Sv1
    p_Sv1_Shold2 <- as.numeric(alpha2/(eta*lambda + alpha2))          # relative probability of moving from Sv1 -> Shold_2d
    p_Sv1_Ev1 <- as.numeric(eta*lambda/(eta*lambda + alpha2))         # relative probability of moving from Sv1 -> Ev1
    # from Shold_2d
    p_Shold2_ <- as.numeric(1 - exp(-eta*lambda-(1/delay2)))          # total probability of moving from Shold_2s
    p_Shold2_Sv2 <- as.numeric((1/delay2)/(eta*lambda + (1/delay2)))  # relative probability of moving from Shold_2d -> Sv2
    p_Shold2_Ev1 <- as.numeric(eta*lambda/(eta*lambda + (1/delay2)))  # relative probability of moving from Shold_2d -> Ev1
    # from Sv_2d
    p_Sv2_Ev2 <- as.numeric(1 - exp(-eta2*lambda))                    # probability of moving Sv2 -> E
    # from E
    p_E_I <- c(rep(1 - exp(-sigma),9))                                # probability of moving E -> I (or Ev1 -> Iv1 or Ev2 -> Iv2)
    # from I
    p_I_ <- 1 - exp(-gamma-h)                               # total probability of moving from I
    p_I_R <- gamma/(gamma + h)                            # relative probability of moving from I -> R
    p_I_H <- h / (gamma + h)                              # relative probability of moving from I -> H
    # from Iv_1d
    p_Iv1_ <- as.numeric(1 - exp(-gamma-eta_hosp*h))                  # total probability of moving from Iv1
    p_Iv1_Hv1 <- as.numeric(eta_hosp*h/(eta_hosp*h + gamma))          # relative probability of moving from Iv1 -> Hv1
    p_Iv1_Rv1 <- as.numeric(gamma/(eta_hosp*h + gamma))               # relative probability of moving from Iv1 -> Rv1
    # from Iv_2d
    p_Iv2_ <- as.numeric(1 - exp(-gamma-eta_hosp2*h))                 # total probability of moving from Iv2
    p_Iv2_Hv2 <- as.numeric(eta_hosp2*h/(eta_hosp2*h + gamma))        # relative probability of moving from Iv2 -> Hv2
    p_Iv2_Rv2 <- as.numeric(gamma/(eta_hosp2*h + gamma))              # relative probability of moving from Iv2 -> Rv2
    # from H
    p_H_ <- 1 - exp(-i1-d-r)                              # total probability of moving from H (or Hv1 or Hv2)
    p_H_IC <- i1/(i1+d+r)                                 # relative probability of moving from H -> IC
    p_H_D <- d/(i1+d+r)                                   # relative probability of moving from H -> D   
    p_H_R <- r/(i1+d+r)                                     # relative probability of moving from H -> R
    # from IC
    p_IC_ <- 1 - exp(-i2-d_ic)                            # total probability of moving from IC
    p_IC_HIC <- i2/(i2+d_ic)                              # relative probability of moving IC -> HIC
    p_IC_D <- d_ic/(i2+d_ic)                             # relative probability of moving from IC -> D
    # from H_IC
    p_HIC_ <- 1 - exp(-d_hic-r_ic)                        # total probability of moving from HIC
    p_HIC_D <- d_hic/(d_hic+r_ic)                         # relative probability of moving from HIC -> D
    p_HIC_R <- r_ic/(d_hic+r_ic)                          # relative probability of moving from HIC -> R
    
    ### number of individuals transitioning between compartments
    # S
    print("S")
    n_S_ <- mapply(FUN = rbinom, n = 1, size = round(S), prob = p_S_)
    x_S_ <- cbind(n_S_, p_S_Shold1, p_S_E)
    n_S_Shold1_E <- apply(x_S_, 1, my_rmultinom)
    # Shold1
    print("Shold1")
    n_Shold1_ <- mapply(FUN = rbinom, n = 1, size = round(Shold_1d), prob = p_Shold1_)
    x_Shold1_ <- cbind(n_Shold1_, p_Shold1_E, p_Shold1_Sv1)
    n_Shold1_E_Sv1 <- apply(x_Shold1_, 1, my_rmultinom)
    # Sv1
    print("Sv1")
    n_Sv1_ <- mapply(FUN = rbinom, n = 1, size = round(Sv_1d), prob = p_Sv1_)
    x_Sv1_ <- cbind(n_Sv1_, p_Sv1_Shold2, p_Sv1_Ev1)
    n_Sv1_Shold2_Ev1 <- apply(x_Sv1_, 1, my_rmultinom)
    # Shold2
    print("Shold2")
    n_Shold2_ <- mapply(FUN = rbinom, n = 1, size = round(Shold_2d), prob = p_Shold2_)
    x_Shold2_ <- cbind(n_Shold2_, p_Shold2_Ev1, p_Shold2_Sv2)
    n_Shold2_Ev1_Sv2 <- apply(x_Shold2_, 1, my_rmultinom)
    # Sv2
    print("Sv2")
    n_Sv2_Ev2 <- mapply(FUN = rbinom, n = 1, size = round(Sv_2d), prob = p_Sv2_Ev2)
    # E
    print("E")
    n_E_I <- mapply(FUN = rbinom, n = 1, size = E, prob = p_E_I)
    n_Ev1_Iv1 <- mapply(FUN = rbinom, n = 1, size = Ev_1d, prob = p_E_I)
    n_Ev2_Iv2 <- mapply(FUN = rbinom, n = 1, size = Ev_2d, prob = p_E_I)
    # I
    print("I")
    n_I_ <- mapply(FUN = rbinom, n = 1, size = round(I), prob = p_I_)
    x_I_ <- cbind(n_I_, p_I_R, p_I_H)
    n_I_R_H <- apply(x_I_, 1, my_rmultinom)
    # Iv1
    print("Iv1")
    n_Iv1_ <- mapply(FUN = rbinom, n = 1, size = round(Iv_1d), prob = p_Iv1_)
    x_Iv1_ <- cbind(n_Iv1_, p_Iv1_Rv1, p_Iv1_Hv1)
    n_Iv1_Rv1_Hv1 <- apply(x_Iv1_, 1, my_rmultinom)
    # Iv2
    print("Iv2")
    n_Iv2_ <- mapply(FUN = rbinom, n = 1, size = round(Iv_2d), prob = p_Iv2_)
    x_Iv2_ <- cbind(n_Iv2_, p_Iv2_Rv2, p_Iv2_Hv2)
    n_Iv2_Rv2_Hv2 <- apply(x_Iv2_, 1, my_rmultinom)
    # H
    print("H")
    n_H_ <- mapply(FUN = rbinom, n = 1, size = round(H), prob = p_H_)
    x_H_ <- cbind(n_H_, p_H_IC, p_H_D, p_H_R)
    n_H_IC_D_R <- apply(x_H_, 1, my_rmultinom)
    # Hv1
    print("Hv1")
    n_Hv1_ <- mapply(FUN = rbinom, n = 1, size = round(Hv_1d), prob = p_H_)
    x_Hv1_ <- cbind(n_Hv1_, p_H_IC, p_H_D, p_H_R)
    n_Hv1_ICv1_D_Rv1 <- apply(x_Hv1_, 1, my_rmultinom)
    # Hv2
    print("Hv2")
    n_Hv2_ <- mapply(FUN = rbinom, n = 1, size = round(Hv_2d), prob = p_H_)
    x_Hv2_ <- cbind(n_Hv2_, p_H_IC, p_H_D, p_H_R)
    n_Hv2_ICv2_D_Rv2 <- apply(x_Hv2_, 1, my_rmultinom)
    # IC
    print("IC")
    n_IC_ <- mapply(FUN = rbinom, n = 1, size = round(IC), prob = p_IC_)
    x_IC_ <- cbind(n_IC_, p_IC_HIC, p_IC_D)
    n_IC_HIC_D <- apply(x_IC_, 1, my_rmultinom)
    # ICv1
    print("ICv1")
    n_ICv1_ <- mapply(FUN = rbinom, n = 1, size = round(ICv_1d), prob = p_IC_)
    x_ICv1_ <- cbind(n_ICv1_, p_IC_HIC, p_IC_D)
    n_ICv1_HICv1_D <- apply(x_ICv1_, 1, my_rmultinom)
    # ICv2
    print("ICv2")
    n_ICv2_ <- mapply(FUN = rbinom, n = 1, size = round(ICv_2d), prob = p_IC_)
    x_ICv2_ <- cbind(n_ICv2_, p_IC_HIC, p_IC_D)
    n_ICv2_HICv2_D <- apply(x_ICv2_, 1, my_rmultinom)
    # H_IC
    print("H_IC")
    n_HIC_ <- mapply(FUN = rbinom, n = 1, size = round(H_IC), prob = p_HIC_)
    n_HIC_D_R <- rmultinom(1, size = n_HIC_, prob = c(p_HIC_D, p_HIC_R))
    # H_ICv1
    print("H_ICv1")
    n_HICv1_ <- mapply(FUN = rbinom, n = 1, size = round(H_ICv_1d), prob = p_HIC_)
    x_HICv1_ <- cbind(n_HICv1_, p_HIC_D, p_HIC_R)
    n_HICv1_D_Rv1 <- apply(x_HICv1_, 1, my_rmultinom)
    # H_ICv2
    print("H_ICv2")
    n_HICv2_ <- mapply(FUN = rbinom, n = 1, size = round(H_ICv_2d), prob = p_HIC_)
    x_HICv2_ <- cbind(n_HICv2_, p_HIC_D, p_HIC_R)
    n_HICv2_D_Rv2 <- apply(x_HICv2_, 1, my_rmultinom)
    
    ################################################################
    # ODEs:
    dS        <- S - n_S_Shold1_E[1,] - n_S_Shold1_E[2,]
    dShold_1d <- Shold_1d + n_S_Shold1_E[1,] - n_Shold1_E_Sv1[1,] - n_Shold1_E_Sv1[2,]
    dSv_1d    <- Sv_1d + n_Shold1_E_Sv1[2,] - n_Sv1_Shold2_Ev1[1,] - n_Sv1_Shold2_Ev1[2,]
    dShold_2d <- Shold_2d - n_Shold2_Ev1_Sv2[1,] - n_Shold2_Ev1_Sv2[2,]
    dSv_2d    <- Sv_2d - n_Sv2_Ev2
    dE        <- E + n_S_Shold1_E[2,] + n_Shold1_E_Sv1[1,] - n_E_I
    dEv_1d    <- Ev_1d + n_Sv1_Shold2_Ev1[2,] + n_Shold2_Ev1_Sv2[1,] - n_Ev1_Iv1
    dEv_2d    <- Ev_2d + n_Sv2_Ev2 - n_Ev2_Iv2
    dI        <- I + n_E_I - n_I_R_H[1,] - n_I_R_H[2,]
    dIv_1d    <- Iv_1d + n_Ev1_Iv1 - n_Iv1_Rv1_Hv1[1,] - n_Iv1_Rv1_Hv1[2,]  
    dIv_2d    <- Iv_2d + n_Ev2_Iv2 - n_Iv2_Rv2_Hv2[1,] - n_Iv2_Rv2_Hv2[2,]
    dH        <- H + n_I_R_H[2,] - n_H_IC_D_R[1,] - n_H_IC_D_R[2,] - n_H_IC_D_R[3,]
    dHv_1d    <- Hv_1d + n_Iv1_Rv1_Hv1[2,] - n_Hv1_ICv1_D_Rv1[1,] - n_Hv1_ICv1_D_Rv1[2,] - n_Hv1_ICv1_D_Rv1[3,]
    dHv_2d    <- Hv_2d + n_Iv1_Rv1_Hv1[2,] - n_Hv2_ICv2_D_Rv2[1,] - n_Hv2_ICv2_D_Rv2[2,]- n_Hv2_ICv2_D_Rv2[3,]
    dIC       <- IC + n_H_IC_D_R[1,] - n_IC_HIC_D[1,] - n_IC_HIC_D[2,]
    dICv_1d   <- ICv_1d + n_Hv1_ICv1_D_Rv1[1,] - n_ICv1_HICv1_D[1,] - n_ICv1_HICv1_D[2,]
    dICv_2d   <- ICv_2d + n_Hv2_ICv2_D_Rv2[1,] - n_ICv2_HICv2_D[1,] - n_ICv2_HICv2_D[2,]
    dH_IC     <- H_IC + n_IC_HIC_D[1,] - n_HIC_D_R[1,] - n_HIC_D_R[2,]
    dH_ICv_1d <- H_ICv_1d + n_ICv1_HICv1_D[1,] - n_HICv1_D_Rv1[1,] - n_HICv1_D_Rv1[2,]
    dH_ICv_2d <- H_ICv_2d + n_ICv2_HICv2_D[1,] - n_HICv2_D_Rv2[1,] - n_HICv2_D_Rv2[2,]
    dD        <- D + n_H_IC_D_R[2,] + n_Hv1_ICv1_D_Rv1[2,] + n_Hv2_ICv2_D_Rv2[2,] +
                     n_IC_HIC_D[2,] + n_ICv1_HICv1_D[2,] + n_ICv2_HICv2_D[2,] +
                     n_HIC_D_R[1,] + n_HICv1_D_Rv1[1,] + n_HICv2_D_Rv2[1,]
    dR        <- R + n_I_R_H[1,] + n_H_IC_D_R[3,] + n_HIC_D_R[2,]
    dRv_1d    <- Rv_1d + n_Iv1_Rv1_Hv1[1,] + n_Hv1_ICv1_D_Rv1[3,] + n_HICv1_D_Rv1[2,]
    dRv_2d    <- Rv_2d + n_Iv2_Rv2_Hv2[1,] + n_Hv2_ICv2_D_Rv2[3,] + n_HICv2_D_Rv2[2,]
    
    # assign variables to global environment, so they can be used for next iteration
    assign("flag_relaxed", flag_relaxed, envir = globalenv())
    assign("flag_very_relaxed", flag_very_relaxed, envir = globalenv())
    assign("flag_normal", flag_normal, envir = globalenv())
    ################################################################
    dt <- 1
    list(c(dt,dS,dShold_1d,dSv_1d,dShold_2d,dSv_2d,dE,dEv_1d,dEv_2d,
           dI,dIv_1d,dIv_2d,dH,dHv_1d,dHv_2d, dH_IC, dH_ICv_1d, dH_ICv_2d,
           dIC, dICv_1d, dICv_2d,dD,dR,dRv_1d,dRv_2d))
  })
}
