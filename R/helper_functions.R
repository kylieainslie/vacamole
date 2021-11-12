# ------------------------------------------------------------------
# Helper data wrangling functions ----------------------------------
# ------------------------------------------------------------------
#' Convert data from wide to long
#' @param x data frame to convert
#' @param times vector of time points. should be the same length as
#' nrow(x)
#' @return data frame in long format
#' @import tidyr
#' @import dplyr
#' @keywords vacamole
#' @export
wide_to_long <- function(x, times) {
  x %>%
    mutate(time = times) %>%
    pivot_longer(
      cols = !.data$time,
      names_to = c("state", "age_group"),
      names_sep = -1,
      values_to = "value"
    )
}

#' Post-process raw results from model run
#' @param x output from model run
#' @return data frame in long format
#' @import tidyr
#' @import dplyr
#' @importFrom stats quantile
#' @keywords vacamole
#' @export
wrangle_results <- function(x) {
  rtn_mle <- x %>%
    filter(.data$sim == 0)

  rtn_bounds <- x %>%
    filter(.data$sim != 0) %>%
    group_by(.data$time, .data$state, .data$age_group) %>%
    summarise(
      lower = stats::quantile(.data$value, probs = 0.025),
      upper = stats::quantile(.data$value, probs = 0.975)
    )

  rtn <- left_join(rtn_mle, rtn_bounds, by = c("time", "state", "age_group")) %>%
    rename(mean = .data$value) %>%
    select(-.data$sim)

  return(rtn)
}
