#TODO refer to manuscript or somewhere else for explanation of each compartment symbol
#' Postprocess output from age-structured SEIR ODE model of vaccination
#' with 2 doses and delay to protection
#' @param dat output from seir model as a data frame
#' @return List of summary results
#' @keywords vacamole
#' @importFrom stringr str_detect
#' @export
postprocess_age_struct_model_output <- function(dat) {
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
  # Shold_3d <- dat %>%
  #   select(names(dat)[str_detect(names(dat), pattern = "Shold_3d[:digit:]")])
  # Sv_3d <- dat %>%
  #   select(names(dat)[str_detect(names(dat), pattern = "Sv_3d[:digit:]")])
  E <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "E[:digit:]")])
  Ev_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Ev_1d[:digit:]")])
  Ev_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Ev_2d[:digit:]")])
  # Ev_3d <- dat %>%
  #   select(names(dat)[str_detect(names(dat), pattern = "Ev_3d[:digit:]")])
  I <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "I[:digit:]")])
  Iv_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Iv_1d[:digit:]")])
  Iv_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Iv_2d[:digit:]")])
  # Iv_3d <- dat %>%
  #   select(names(dat)[str_detect(names(dat), pattern = "Iv_3d[:digit:]")])
  H <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H[:digit:]")])
  Hv_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Hv_1d[:digit:]")])
  Hv_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Hv_2d[:digit:]")])
  # Hv_3d <- dat %>%
  #   select(names(dat)[str_detect(names(dat), pattern = "Hv_3d[:digit:]")])
  H_IC <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H_IC[:digit:]")])
  H_ICv_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H_ICv_1d[:digit:]")])
  H_ICv_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H_ICv_2d[:digit:]")])
  # H_ICv_3d <- dat %>%
  #   select(names(dat)[str_detect(names(dat), pattern = "H_ICv_3d[:digit:]")])
  IC <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "IC[:digit:]")])
  ICv_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "ICv_1d[:digit:]")])
  ICv_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "ICv_2d[:digit:]")])
  # ICv_3d <- dat %>%
  #   select(names(dat)[str_detect(names(dat), pattern = "ICv_3d[:digit:]")])
  D <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "D[:digit:]")])
  R <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "R[:digit:]")])
  Rv_1d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_1d[:digit:]")])
  Rv_2d <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_2d[:digit:]")])
  # Rv_3d <- dat %>%
  #   select(names(dat)[str_detect(names(dat), pattern = "Rv_3d[:digit:]")])

  rtn <- list(
    S = S,
    Shold_1d = Shold_1d,
    Sv_1d = Sv_1d,
    Shold_2d = Shold_2d,
    Sv_2d = Sv_2d,
    # Shold_3d = Shold_3d,
    # Sv_3d = Sv_3d,
    E = E,
    Ev_1d = Ev_1d,
    Ev_2d = Ev_2d,
    # Ev_3d = Ev_3d,
    I = I,
    Iv_1d = Iv_1d,
    Iv_2d = Iv_2d,
    # Iv_3d = Iv_3d,
    H = H,
    Hv_1d = Hv_1d,
    Hv_2d = Hv_2d,
    # Hv_3d = Hv_3d,
    H_IC = H_IC,
    H_ICv_1d = H_ICv_1d,
    H_ICv_2d = H_ICv_2d,
    # H_ICv_3d = H_ICv_3d,
    IC = IC[, -c(1:9)],
    ICv_1d = ICv_1d[, -c(1:9)],
    ICv_2d = ICv_2d[, -c(1:9)],
    # ICv_3d = ICv_3d[, -c(1:9)],
    D = D,
    R = R,
    Rv_1d = Rv_1d,
    Rv_2d = Rv_2d
    # Rv_3d = Rv_3d
  )

  return(rtn)
}
