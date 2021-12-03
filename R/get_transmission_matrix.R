#' Convert contact matrix to transmission matrix
#' @param x vector of relative susceptibility and infectiousnes for each age group
#' @param contact_mat contact matrix
#' @return transmission matrix
#' @keywords vacamole
#' @export
get_transmission_matrix <- function(x, contact_mat) {
  # multiply by relative susc/inf
  tmp <- sweep(contact_mat, 1, x, "*") # rows
  rtn <- sweep(tmp, 2, x, "*") # columns

  # output
  return(rtn)
}
