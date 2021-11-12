#' Convert contact matrices to transmission matrices
#' @param x data frame with contacts
#' @return list of contact matrices
#' @import dplyr
#' @import tidyr
#' @importFrom plyr aaply
#' @importFrom plyr laply
#' @keywords vacamole
#' @export
convert_contact_matrices <- function(x) {
  # create empty list to store converted matrices
  rtn <- list()
  # loop over different realisations of contact matrices
  var_names <- paste0("c_smt.", 1:200)

  for (i in 1:200) {
    tmp <- x %>%
      select(.data$part_age, .data$cnt_age, var_names[i]) %>%
      mutate(contact_type = c(
        rep("all", 81),
        rep("community", 81),
        rep("household", 81)
      )) %>%
      filter(contact_type == "all") %>%
      select(-.data$contact_type) %>%
      pivot_wider(., names_from = cnt_age, values_from = var_names[i]) %>%
      select(-.data$part_age)

    # convert to transmission matrix
    tmp1 <- as.matrix(tmp) %*% N_diag
    tmp2 <- get_transmission_matrix(rel_trans, tmp1)

    rtn[[i]] <- tmp2
  }

  # add mean to list of matrices
  tmp_mean <- aaply(laply(rtn, as.matrix), c(2, 3), mean)
  rtn$mean <- tmp_mean

  return(rtn)
}
