#' Age-structured SEIR ODE model of vaccination with 2 doses and delay to protection
#' @param times vector of times
#' @param init list of initial states
#' @param params list of parameter values
#' @return List of summary results
#' @keywords vacamole
#' @export
# Define model -----------------------------------------------------
age_struct_seir_ode <- function(times, init, params) {
  with(as.list(c(params, init)), {
    # print(t)
    # define initial state vectors from input ----------------------
    S <- c(S1, S2, S3, S4, S5, S6, S7, S8, S9)
    Shold_1d <- c(
      Shold_1d1, Shold_1d2, Shold_1d3, Shold_1d4, Shold_1d5, Shold_1d6,
      Shold_1d7, Shold_1d8, Shold_1d9
    )
    Sv_1d <- c(Sv_1d1, Sv_1d2, Sv_1d3, Sv_1d4, Sv_1d5, Sv_1d6, Sv_1d7, Sv_1d8, Sv_1d9)
    Shold_2d <- c(
      Shold_2d1, Shold_2d2, Shold_2d3, Shold_2d4, Shold_2d5, Shold_2d6,
      Shold_2d7, Shold_2d8, Shold_2d9
    )
    Sv_2d <- c(Sv_2d1, Sv_2d2, Sv_2d3, Sv_2d4, Sv_2d5, Sv_2d6, Sv_2d7, Sv_2d8, Sv_2d9)
    E <- c(E1, E2, E3, E4, E5, E6, E7, E8, E9)
    Ev_1d <- c(Ev_1d1, Ev_1d2, Ev_1d3, Ev_1d4, Ev_1d5, Ev_1d6, Ev_1d7, Ev_1d8, Ev_1d9)
    Ev_2d <- c(Ev_2d1, Ev_2d2, Ev_2d3, Ev_2d4, Ev_2d5, Ev_2d6, Ev_2d7, Ev_2d8, Ev_2d9)
    I <- c(I1, I2, I3, I4, I5, I6, I7, I8, I9)
    Iv_1d <- c(Iv_1d1, Iv_1d2, Iv_1d3, Iv_1d4, Iv_1d5, Iv_1d6, Iv_1d7, Iv_1d8, Iv_1d9)
    Iv_2d <- c(Iv_2d1, Iv_2d2, Iv_2d3, Iv_2d4, Iv_2d5, Iv_2d6, Iv_2d7, Iv_2d8, Iv_2d9)
    H <- c(H1, H2, H3, H4, H5, H6, H7, H8, H9)
    Hv_1d <- c(Hv_1d1, Hv_1d2, Hv_1d3, Hv_1d4, Hv_1d5, Hv_1d6, Hv_1d7, Hv_1d8, Hv_1d9)
    Hv_2d <- c(Hv_2d1, Hv_2d2, Hv_2d3, Hv_2d4, Hv_2d5, Hv_2d6, Hv_2d7, Hv_2d8, Hv_2d9)
    H_IC <- c(H_IC1, H_IC2, H_IC3, H_IC4, H_IC5, H_IC6, H_IC7, H_IC8, H_IC9)
    H_ICv_1d <- c(H_ICv_1d1, H_ICv_1d2, H_ICv_1d3, H_ICv_1d4, H_ICv_1d5, H_ICv_1d6, H_ICv_1d7, H_ICv_1d8, H_ICv_1d9)
    H_ICv_2d <- c(H_ICv_2d1, H_ICv_2d2, H_ICv_2d3, H_ICv_2d4, H_ICv_2d5, H_ICv_2d6, H_ICv_2d7, H_ICv_2d8, H_ICv_2d9)
    IC <- c(IC1, IC2, IC3, IC4, IC5, IC6, IC7, IC8, IC9)
    ICv_1d <- c(ICv_1d1, ICv_1d2, ICv_1d3, ICv_1d4, ICv_1d5, ICv_1d6, ICv_1d7, ICv_1d8, ICv_1d9)
    ICv_2d <- c(ICv_2d1, ICv_2d2, ICv_2d3, ICv_2d4, ICv_2d5, ICv_2d6, ICv_2d7, ICv_2d8, ICv_2d9)
    D <- c(D1, D2, D3, D4, D5, D6, D7, D8, D9)
    R <- c(R1, R2, R3, R4, R5, R6, R7, R8, R9)
    Rv_1d <- c(Rv_1d1, Rv_1d2, Rv_1d3, Rv_1d4, Rv_1d5, Rv_1d6, Rv_1d7, Rv_1d8, Rv_1d9)
    Rv_2d <- c(Rv_2d1, Rv_2d2, Rv_2d3, Rv_2d4, Rv_2d5, Rv_2d6, Rv_2d7, Rv_2d8, Rv_2d9)

    # determine vaccination rate -----------------------------------
    if (!no_vac) {
      index <- floor(times) + 1

      alpha <- vac_inputs$alpha_dose1[index, ]
      alpha2 <- vac_inputs$alpha_dose2[index, ]
      eta <- vac_inputs$eta_dose1[index, ]
      eta2 <- vac_inputs$eta_dose2[index, ]
      delay <- vac_inputs$delay_dose1[index, ]
      delay2 <- vac_inputs$delay_dose2[index, ]
      eta_hosp <- vac_inputs$eta_hosp_dose1[index, ]
      eta_hosp2 <- vac_inputs$eta_hosp_dose2[index, ]
      eta_trans <- vac_inputs$eta_trans_dose1[index, ]
      eta_trans2 <- vac_inputs$eta_trans_dose2[index, ]
    } else {
      alpha <- 0
      alpha2 <- 0
      eta <- 1
      eta2 <- 1
      delay <- 1
      delay2 <- 1
      eta_hosp <- 1
      eta_hosp2 <- 1
      eta_trans <- 1
      eta_trans2 <- 1
    }

    # determine contact matrix based on criteria --------------------
    ic_admin <- sum(i1 * (H + Hv_1d + Hv_2d))

    cases <- sum(sigma * (E + Ev_1d + Ev_2d) * p_report)
    criteria <- (use_cases) * cases + (!use_cases) * ic_admin

    # initialize flags
    if (times == 0 | params$keep_cm_fixed) {
      flag_relaxed <- 0
      flag_very_relaxed <- 0
      flag_normal <- 0
    }

    # determine contact matrix to use based on criteria
    tmp2 <- choose_contact_matrix(
      times = t,
      params = params,
      criteria = criteria,
      flag_relaxed = flag_relaxed,
      flag_very_relaxed = flag_very_relaxed,
      flag_normal = flag_normal,
      keep_fixed = keep_cm_fixed
    )

    contact_mat <- tmp2$contact_matrix
    flag_relaxed <- tmp2$flag_relaxed
    flag_very_relaxed <- tmp2$flag_very_relaxed
    flag_normal <- tmp2$flag_normal

    if (flag_normal > 0 & !is.null(beta_change)) {
      beta <- beta_change
    }

    # determine force of infection ----------------------------------
    calendar_day <- ifelse(times > 365, t_calendar_start + times - 365, t_calendar_start + times)
    beta_t <- beta * (1 + beta1 * cos(2 * pi * calendar_day / 365.24)) # incorporate seasonality in transmission rate
    lambda <- beta_t * (contact_mat %*% (I + (eta_trans * Iv_1d) + (eta_trans2 * Iv_2d)))
    lambda <- ifelse(lambda < 0, 0, lambda)
    # ---------------------------------------------------------------

    ################################################################
    # ODEs:
    dS <- -lambda * S - alpha * S
    dShold_1d <- alpha * S - (1 / delay) * Shold_1d - lambda * Shold_1d
    dSv_1d <- (1 / delay) * Shold_1d - eta * lambda * Sv_1d - alpha2 * Sv_1d
    dShold_2d <- alpha2 * Sv_1d - (1 / delay2) * Shold_2d - eta * lambda * Shold_2d
    dSv_2d <- (1 / delay2) * Shold_2d - eta2 * lambda * Sv_2d
    dE <- lambda * (S + Shold_1d) - sigma * E + epsilon
    dEv_1d <- eta * lambda * (Sv_1d + Shold_2d) - sigma * Ev_1d
    dEv_2d <- eta2 * lambda * Sv_2d - sigma * Ev_2d
    dI <- sigma * E - (gamma + h) * I
    dIv_1d <- sigma * Ev_1d - (gamma + eta_hosp * h) * Iv_1d
    dIv_2d <- sigma * Ev_2d - (gamma + eta_hosp2 * h) * Iv_2d
    dH <- h * I - (i1 + d + r) * H
    dHv_1d <- eta_hosp * h * Iv_1d - (i1 + d + r) * Hv_1d
    dHv_2d <- eta_hosp2 * h * Iv_2d - (i1 + d + r) * Hv_2d
    dH_IC <- i2 * IC - (r_ic + d_hic) * H_IC
    dH_ICv_1d <- i2 * ICv_1d - (r_ic + d_hic) * H_ICv_1d
    dH_ICv_2d <- i2 * ICv_2d - (r_ic + d_hic) * H_ICv_2d
    dIC <- i1 * H - (i2 + d_ic) * IC
    dICv_1d <- i1 * Hv_1d - (i2 + d_ic) * ICv_1d
    dICv_2d <- i1 * Hv_2d - (i2 + d_ic) * ICv_2d
    dD <- d * (H + Hv_1d + Hv_2d) + d_ic * (IC + ICv_1d + ICv_2d) + d_hic * (H_IC + H_ICv_1d + H_ICv_2d)
    dR <- gamma * I + r * H + r_ic * H_IC
    dRv_1d <- gamma * Iv_1d + r * Hv_1d + r_ic * H_ICv_1d
    dRv_2d <- gamma * Iv_2d + r * Hv_2d + r_ic * H_ICv_2d

    # assign variables to global environment, so they can be used for next iteration
    assign("flag_relaxed", flag_relaxed, envir = globalenv())
    assign("flag_very_relaxed", flag_very_relaxed, envir = globalenv())
    assign("flag_normal", flag_normal, envir = globalenv())
    ################################################################
    dt <- 1
    list(c(
      dt, dS, dShold_1d, dSv_1d, dShold_2d, dSv_2d, dE, dEv_1d, dEv_2d,
      dI, dIv_1d, dIv_2d, dH, dHv_1d, dHv_2d, dH_IC, dH_ICv_1d, dH_ICv_2d,
      dIC, dICv_1d, dICv_2d, dD, dR, dRv_1d, dRv_2d
    ))
  })
}
