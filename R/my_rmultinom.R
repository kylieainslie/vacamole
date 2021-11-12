#' Postprocess output from age-structured SEIR ODE model of vaccination
#' with 2 doses and delay to protection
#' @param x vector of probabilities
#' @return List of summary results
#' @keywords vacamole
#' @importFrom stats rmultinom
#' @export
my_rmultinom <- function(x) {
  size <- x[1]
  p1 <- x[2]
  p2 <- x[3]

  if (length(x) == 4) {
    p3 <- x[4]
    prob_vec <- c(p1, p2, p3)
  } else {
    prob_vec <- c(p1, p2)
  }
  # print(size)
  # print(prob_vec)
  rmultinom(n = 1, size = size, prob = prob_vec)
}
