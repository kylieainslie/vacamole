library(cubature)

# determine generation interval by determining the convolution
# of the latent and infectious periods
pdf_latent <- function(x) dweibull(x, shape = 1.97, scale = 9.09)
#pdf_infectious_a <- function(y) dgamma(y, shape = 2.1, scale = 4)
pdf_infectious <- function(z) dgamma(z, shape = 2.4, scale = 4)

pdf_generation <- function(t){
  cubintegrate(function(x){pdf_latent(x = x) * pdf_infectious(z = t - x)},
               lower = c(0),
               upper = c(Inf)
               )$integral
}

vpdf_generation <- Vectorize(pdf_generation, vectorize.args = "t")
t = seq(0, 112, 1)
plot(t, vpdf_generation(t), type = "l")

# estimate Rt using incidence data and generation interval
library(EpiEstim)

results_y2o <- readRDS("results/model_results_beta03_AZ_young_to_old.rds")
results_o2y <- readRDS("results/model_results_beta03_AZ_old_to_young.rds")
results_60_plus <- readRDS("results/model_results_beta03_AZ_60_plus_first.rds")


sim_scenario <- results_o2y
vaccines_per_day <- read_xlsx("inst/extdata/data/vaccines_per_day_by_type.xlsx")
inc_mat <- sim_scenario$new_infections[,-1]
incidence <- apply(inc_mat, 2, sum)

inc_dat <- data.frame(vaccines_per_day$date[-1], 
                      I = incidence)

res_non_parametric_si <- estimate_R(inc_dat, 
                                    method="non_parametric_si",
                                    config = make_config(list(
                                      si_distr = vpdf_generation(t)))
                                   )
plot(res_non_parametric_si, legend = FALSE, type = "R")

# plot hospitalisations
hosp_o2y <- old_to_young$total_hospitalisations
hosp_y2o <- young_to_old$total_hospitalisations


