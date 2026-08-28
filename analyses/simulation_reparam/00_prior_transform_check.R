#===============================================================================
# Prior transformation check for alternative temporal parameterisations
#===============================================================================

#-------------------------------------------------------------------------------
# Load experimental package library
#-------------------------------------------------------------------------------

library(here)

reparam_lib <- normalizePath(here(".Rlib-reparam"), mustWork = TRUE)
.libPaths(c(reparam_lib, .libPaths()))

library(ETAS.inlabru, lib.loc = reparam_lib)

stopifnot(
  dirname(normalizePath(find.package("ETAS.inlabru"))) == reparam_lib)

#-------------------------------------------------------------------------------
# Sampling variables
#-------------------------------------------------------------------------------

set.seed(1)
n_draw <- 10000
bound_eps <- 1e-6

#-------------------------------------------------------------------------------
# MSE: sample from original P0 prior and transform
#-------------------------------------------------------------------------------

mse_native <- list(
  d = runif(n_draw, bound_eps, 1),
  rho = rlnorm(n_draw, meanlog = log(1.5), sdlog = 1),
  gamma = runif(n_draw, bound_eps, 1 - bound_eps))

mse_reparam <- ETAS.inlabru:::temporal_reparam_from_native(mse_native, "mse")
mse_reparam <- as.data.frame(mse_reparam)

#-------------------------------------------------------------------------------
# Rate-state: sample from original P0 prior and transform
#-------------------------------------------------------------------------------

rs_native <- list(
  B = plogis(rnorm(n_draw, mean = 7.5, sd = 1)),
  ta = rlnorm(n_draw, meanlog = log(200), sdlog = 0.5))

rs_reparam <- ETAS.inlabru:::temporal_reparam_from_native(rs_native, "rate_state")
rs_reparam <- as.data.frame(rs_reparam)

#-------------------------------------------------------------------------------
# Marginal summaries and transformed-prior dependence
#-------------------------------------------------------------------------------

summarise_parameter <- function(x) {
  c(q05 = unname(quantile(x, 0.05)),
    median = unname(median(x)),
    q95 = unname(quantile(x, 0.95)),
    meanlog = mean(log(x)),
    sdlog = sd(log(x)))
}

cat("\nMSE transformed priors\n")
print(t(sapply(mse_reparam, summarise_parameter)))

cat("\nMSE Spearman correlations\n")
print(round(cor(mse_reparam, method = "spearman"), 3))

cat("\nRate-state transformed priors\n")
print(t(sapply(rs_reparam, summarise_parameter)))

cat("\nRate-state Spearman correlations\n")
print(round(cor(rs_reparam, method = "spearman"), 3))