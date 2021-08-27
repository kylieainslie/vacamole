# --------------------------------------------------
# extra code bits
# --------------------------------------------------
# segmentation analysis ----------------------------
# my_glm <- glm(inc ~ date, data = osiris1)
# my_coef <- coef(my_glm)
# 
# p1 <- p + geom_abline(
#   intercept = my_coef[1],
#   slope = my_coef[2]
# )
# p1
# 
# # now for the actual breakpoint analysis
# my_seg <- segmented(my_glm,
#   seg.Z = ~date,
#   psi = list(date = c(
#     as.Date("2021-03-07"),
#     as.Date("2021-03-21"),
#     as.Date("2021-04-01"),
#     as.Date("2021-04-17")
#   ))
# )
# summary(my_seg)
# # the breakpoints
# my_seg$psi
# 
# # the slopes
# my_slopes <- slope(my_seg)
# # beta_mult <- (my_slopes$date[,1]/my_slopes$date[1,1])
# 
# #  get the fitted data
# my_fitted <- fitted(my_seg)
# my_model <- data.frame(Date = osiris1$date, Daily_Cases = my_fitted)
# p + geom_line(data = my_model, aes(x = Date, y = Daily_Cases), color = "blue")
# --------------------------------------------------
# Sangeeta's code from hermione
# log_likelihood <- function(t, inf_params, ip_params, fun, ...) {
#   log(fun(t, inf_params, ip_params, ...))
# }
#
# out <- mapply(
#   FUN = log_likelihood,
#   t = tvec,
#   offset = offset_vec,
#   MoreArgs = list(
#     inf_params = inf_params,
#     ip_params = ip_params,
#     fun = probability_offset
#   ),
#   SIMPLIFY = TRUE
# )
# ---------------------------------------------------
# by hand likelihood approach - grid search
# --------------------------------------------------
# my_grid <- expand.grid(
#   beta = seq(0.00001, 0.001, by = 0.00001),
#   alpha = seq(0, 1, by = 0.01)
# )
# 
# out <- apply(my_grid, 1, likelihood_func,
#              # other args for likelihood_fun
#              t = seq(0, breakpoints[1], by = 1),
#              data = osiris1,
#              params = params,
#              init = init
# )
# 
# out <- mapply(
#   FUN = likelihood_func,
#   beta = as.list(seq(0.00001, 0.001, by = 0.00001)),
#   alpha = as.list(seq(0, 1, by = 0.01)),
#   t = list(seq(0, breakpoints[1], by = 1)),
#   MoreArgs = list(
#     data = osiris1,
#     params = params,
#     init = init
#   ) # ,
#   # SIMPLIFY = TRUE
# )
