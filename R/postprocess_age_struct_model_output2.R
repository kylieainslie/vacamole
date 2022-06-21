#TODO refer to manuscript or somewhere else for explanation of each compartment symbol
#' Postprocess output from age-structured SEIR ODE model of vaccination
#' with 2 doses and delay to protection
#' @param dat output from seir model as a data frame
#' @return List of summary results
#' @keywords vacamole
#' @importFrom stringr str_detect
#' @export
postprocess_age_struct_model_output2 <- function(dat) {
  # S ----------
  S <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "S[:digit:]")])
  Shold_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Shold_1d[:digit:]")])
  Sv_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Sv_1d[:digit:]")])
  Shold_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Shold_2d[:digit:]")])
  Sv_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Sv_2d[:digit:]")])
  Shold_3d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Shold_3d[:digit:]")])
  Sv_3d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Sv_3d[:digit:]")])
  Shold_4d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Shold_4d[:digit:]")])
  Sv_4d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Sv_4d[:digit:]")])
  Shold_5d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Shold_5d[:digit:]")])
  Sv_5d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Sv_5d[:digit:]")])
  
  # E ----------
  E <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "E[:digit:]")])
  Ev_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Ev_1d[:digit:]")])
  Ev_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Ev_2d[:digit:]")])
  Ev_3d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Ev_3d[:digit:]")])
  Ev_4d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Ev_4d[:digit:]")])
  Ev_5d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Ev_5d[:digit:]")])
  
  # I ----------
  I <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "I[:digit:]")])
  Iv_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Iv_1d[:digit:]")])
  Iv_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Iv_2d[:digit:]")])
  Iv_3d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Iv_3d[:digit:]")])
  Iv_4d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Iv_4d[:digit:]")])
  Iv_5d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Iv_5d[:digit:]")])
  
  # H ----------
  H <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H[:digit:]")])
  Hv_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Hv_1d[:digit:]")])
  Hv_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Hv_2d[:digit:]")])
  Hv_3d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Hv_3d[:digit:]")])
  Hv_4d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Hv_4d[:digit:]")])
  Hv_5d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Hv_5d[:digit:]")])
  
  # H_IC -------
  H_IC <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H_IC[:digit:]")])
  H_ICv_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H_ICv_1d[:digit:]")])
  H_ICv_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H_ICv_2d[:digit:]")])
  H_ICv_3d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H_ICv_3d[:digit:]")])
  H_ICv_4d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H_ICv_4d[:digit:]")])
  H_ICv_5d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H_ICv_5d[:digit:]")])
  
  # IC ---------
  IC <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "IC[:digit:]")]) %>%
    select(-names(dat)[str_detect(names(dat), pattern = "H_IC[:digit:]")])
  ICv_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "ICv_1d[:digit:]")]) %>%
    select(-names(dat)[str_detect(names(dat), pattern = "H_ICv_1d[:digit:]")])
  ICv_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "ICv_2d[:digit:]")]) %>%
    select(-names(dat)[str_detect(names(dat), pattern = "H_ICv_2d[:digit:]")])
  ICv_3d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "ICv_3d[:digit:]")]) %>%
    select(-names(dat)[str_detect(names(dat), pattern = "H_ICv_3d[:digit:]")])
  ICv_4d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "ICv_4d[:digit:]")]) %>%
    select(-names(dat)[str_detect(names(dat), pattern = "H_ICv_4d[:digit:]")])
  ICv_5d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "ICv_5d[:digit:]")]) %>%
    select(-names(dat)[str_detect(names(dat), pattern = "H_ICv_5d[:digit:]")])
  
  # D ----------
  D <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "D[:digit:]")])
  
  # R ----------
  R <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "R[:digit:]")])
  Rv_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_1d[:digit:]")])
  Rv_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_2d[:digit:]")])
  Rv_3d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_3d[:digit:]")])
  Rv_4d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_4d[:digit:]")])
  Rv_5d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_5d[:digit:]")])
  
  # R_1w ------
  R_1w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "R_1w[:digit:]")])
  Rv_1d_1w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_1d_1w[:digit:]")])
  Rv_2d_1w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_2d_1w[:digit:]")])
  Rv_3d_1w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_3d_1w[:digit:]")])
  Rv_4d_1w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_4d_1w[:digit:]")])
  Rv_5d_1w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_5d_1w[:digit:]")])
  
  # R_2w ------
  R_2w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "R_2w[:digit:]")])
  Rv_1d_2w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_1d_2w[:digit:]")])
  Rv_2d_2w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_2d_2w[:digit:]")])
  Rv_3d_2w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_3d_2w[:digit:]")])
  Rv_4d_2w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_4d_2w[:digit:]")])
  Rv_5d_2w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_5d_2w[:digit:]")])
  
  # R_3w ------
  R_3w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "R_3w[:digit:]")])
  Rv_1d_3w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_1d_3w[:digit:]")])
  Rv_2d_3w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_2d_3w[:digit:]")])
  Rv_3d_3w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_3d_3w[:digit:]")])
  Rv_4d_3w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_4d_3w[:digit:]")])
  Rv_5d_3w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_5d_3w[:digit:]")])
  
  # output ----
  rtn <- list(
    S = S,
    Shold_1d = Shold_1d,
    Sv_1d = Sv_1d,
    Shold_2d = Shold_2d,
    Sv_2d = Sv_2d,
    Shold_3d = Shold_3d,
    Sv_3d = Sv_3d,
    Shold_4d = Shold_4d,
    Sv_4d = Sv_4d,
    Shold_5d = Shold_5d,
    Sv_5d = Sv_5d,
    E = E,
    Ev_1d = Ev_1d,
    Ev_2d = Ev_2d,
    Ev_3d = Ev_3d,
    Ev_4d = Ev_4d,
    Ev_5d = Ev_5d,
    I = I,
    Iv_1d = Iv_1d,
    Iv_2d = Iv_2d,
    Iv_3d = Iv_3d,
    Iv_4d = Iv_4d,
    Iv_5d = Iv_5d,
    H = H,
    Hv_1d = Hv_1d,
    Hv_2d = Hv_2d,
    Hv_3d = Hv_3d,
    Hv_4d = Hv_4d,
    Hv_5d = Hv_5d,
    H_IC = H_IC,
    H_ICv_1d = H_ICv_1d,
    H_ICv_2d = H_ICv_2d,
    H_ICv_3d = H_ICv_3d,
    H_ICv_4d = H_ICv_4d,
    H_ICv_5d = H_ICv_5d,
    IC = IC,
    ICv_1d = ICv_1d,
    ICv_2d = ICv_2d,
    ICv_3d = ICv_3d,
    ICv_4d = ICv_4d,
    ICv_5d = ICv_5d,
    D = D,
    R = R,
    Rv_1d = Rv_1d,
    Rv_2d = Rv_2d,
    Rv_3d = Rv_3d,
    Rv_4d = Rv_4d,
    Rv_5d = Rv_5d,
    R_1w = R_1w,
    Rv_1d_1w = Rv_1d_1w,
    Rv_2d_1w = Rv_2d_1w,
    Rv_3d_1w = Rv_3d_1w,
    Rv_4d_1w = Rv_4d_1w,
    Rv_5d_1w = Rv_5d_1w,
    R_2w = R_2w,
    Rv_1d_2w = Rv_1d_2w,
    Rv_2d_2w = Rv_2d_2w,
    Rv_3d_2w = Rv_3d_2w,
    Rv_4d_2w = Rv_4d_2w,
    Rv_5d_2w = Rv_5d_2w,
    R_3w = R_3w,
    Rv_1d_3w = Rv_1d_3w,
    Rv_2d_3w = Rv_2d_3w,
    Rv_3d_3w = Rv_3d_3w,
    Rv_4d_3w = Rv_4d_3w,
    Rv_5d_3w = Rv_5d_3w
  )

  return(rtn)
}
