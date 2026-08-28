#===============================================================================
# Stage A reparameterised fitting design
#===============================================================================

#-------------------------------------------------------------------------------
# Load experimental package and baseline simulation design
#-------------------------------------------------------------------------------

library(here)

reparam_lib <- normalizePath(here(".Rlib-reparam"), mustWork = TRUE)
.libPaths(c(reparam_lib, .libPaths()))

library(ETAS.inlabru, lib.loc = reparam_lib)

stopifnot(
  dirname(normalizePath(find.package("ETAS.inlabru"))) == reparam_lib)

source(here("analyses", "simulation", "00_design.R"))

stopifnot(
  dirname(normalizePath(find.package("ETAS.inlabru"))) == reparam_lib)

#-------------------------------------------------------------------------------
# Stage-A marginal prior approximations
# Independent lognormal marginal approximations based on prior draws from 00
#-------------------------------------------------------------------------------

prior_stageA <- list(
  
  mse = list(
    decay_mse = list(meanlog = 0.9823922, sdlog = 0.9755191),
    shape_mse = list(meanlog = -0.1198424, sdlog = 1.1950273),
    gamma = prior_baseline$mse$gamma),
  
  rate_state = list(
    decay_rs = list(meanlog = 2.212921, sdlog = 1.1313568),
    shape_rs = list(meanlog = 2.006779, sdlog = 0.1373686))
  )

#-------------------------------------------------------------------------------
# Forward links: latent Gaussian scale to Stage-A parameter scale
#-------------------------------------------------------------------------------

make_links_stageA <- function(kernel) {
  
  stopifnot(kernel %in% c("mse", "rate_state"))
  
  common <- list(
    mu = \(x) gamma_t(x, prior_baseline$mu$shape, prior_baseline$mu$rate),
    K = \(x) loggaus_t(x, prior_baseline$K$meanlog, prior_baseline$K$sdlog),
    alpha = \(x) unif_t(x, prior_baseline$alpha$min, prior_baseline$alpha$max))
  
  if (kernel == "mse") {
    return(c(common, list(
      decay_mse = \(x) loggaus_t(x, prior_stageA$mse$decay_mse$meanlog, 
                                 prior_stageA$mse$decay_mse$sdlog),
      shape_mse = \(x) loggaus_t(x, prior_stageA$mse$shape_mse$meanlog, 
                                 prior_stageA$mse$shape_mse$sdlog),
      gamma = \(x) unif_t(x, prior_stageA$mse$gamma$min, 
                          prior_stageA$mse$gamma$max))))
  }
  
  if (kernel == "rate_state") {
    return(c(common, list(
      decay_rs = \(x) loggaus_t(x, prior_stageA$rate_state$decay_rs$meanlog, 
                                prior_stageA$rate_state$decay_rs$sdlog),
      shape_rs = \(x) loggaus_t(x, prior_stageA$rate_state$shape_rs$meanlog, 
                                prior_stageA$rate_state$shape_rs$sdlog))))
  }
}

#-------------------------------------------------------------------------------
# Inverse links: Stage-A parameter scale to latent Gaussian scale
#-------------------------------------------------------------------------------

make_inverse_links_stageA <- function(kernel) {
  
  stopifnot(kernel %in% c("mse", "rate_state"))
  
  common <- list(
    mu = \(x) inv_gamma_t(x, prior_baseline$mu$shape, prior_baseline$mu$rate),
    K = \(x) inv_loggaus_t(x, prior_baseline$K$meanlog, prior_baseline$K$sdlog),
    alpha = \(x) inv_unif_t(x, prior_baseline$alpha$min, prior_baseline$alpha$max))
  
  if (kernel == "mse") {
    return(c(common, list(
      decay_mse = \(x) inv_loggaus_t(x, prior_stageA$mse$decay_mse$meanlog,
                                     prior_stageA$mse$decay_mse$sdlog),
      shape_mse = \(x) inv_loggaus_t(x, prior_stageA$mse$shape_mse$meanlog, 
                                     prior_stageA$mse$shape_mse$sdlog),
      gamma = \(x) inv_unif_t(x, prior_stageA$mse$gamma$min,
                              prior_stageA$mse$gamma$max))))
  }
  
  if (kernel == "rate_state") {
    return(c(common, list(
      decay_rs = \(x) inv_loggaus_t(x, prior_stageA$rate_state$decay_rs$meanlog, 
                                    prior_stageA$rate_state$decay_rs$sdlog),
      shape_rs = \(x) inv_loggaus_t(x, prior_stageA$rate_state$shape_rs$meanlog, 
                                    prior_stageA$rate_state$shape_rs$sdlog))))
  }
}

#-------------------------------------------------------------------------------
# Original fitting starting values expressed in Stage-A coordinates
#-------------------------------------------------------------------------------

mse_start <- ETAS.inlabru:::temporal_reparam_from_native(
  as.list(initials$mse[c("d", "rho", "gamma")]), "mse")

rs_start <- ETAS.inlabru:::temporal_reparam_from_native(
  as.list(initials$rate_state[c("B", "ta")]), "rate_state")

initials_stageA <- list(
  
  mse = c(initials$mse[c("mu", "K", "alpha")], decay_mse = mse_start$decay_mse,
          shape_mse = mse_start$shape_mse,gamma = mse_start$gamma),
  
  rate_state = c(initials$rate_state[c("mu", "K", "alpha")], 
                 decay_rs = rs_start$decay_rs, shape_rs = rs_start$shape_rs)
)

#-------------------------------------------------------------------------------
# inlabru fitting options for Stage A
#-------------------------------------------------------------------------------

make_bru_options_stageA <- function(kernel, rel_tol = 0.1, max_iter = 100) {
  
  stopifnot(kernel %in% c("mse", "rate_state"))
  
  inv <- make_inverse_links_stageA(kernel)
  init <- initials_stageA[[kernel]]
  
  if (kernel == "mse") {
    th_init <- list(
      th.mu = inv$mu(init["mu"]),
      th.K = inv$K(init["K"]),
      th.alpha = inv$alpha(init["alpha"]),
      th.decay_mse = inv$decay_mse(init["decay_mse"]),
      th.shape_mse = inv$shape_mse(init["shape_mse"]),
      th.gamma = inv$gamma(init["gamma"]))
  }
  
  if (kernel == "rate_state") {
    th_init <- list(
      th.mu = inv$mu(init["mu"]),
      th.K = inv$K(init["K"]),
      th.alpha = inv$alpha(init["alpha"]),
      th.decay_rs = inv$decay_rs(init["decay_rs"]),
      th.shape_rs = inv$shape_rs(init["shape_rs"]))
  }
  
  list(bru_verbose = 0, bru_max_iter = max_iter, bru_initial = th_init,
       bru_rel_tol = rel_tol)
}

#-------------------------------------------------------------------------------
# Check forward and inverse links recover Stage-A starting values
#-------------------------------------------------------------------------------

for (kernel in c("mse", "rate_state")) {
  
  links <- make_links_stageA(kernel)
  inv <- make_inverse_links_stageA(kernel)
  init <- initials_stageA[[kernel]]
  
  latent <- sapply(names(init), \(name) inv[[name]](init[[name]]))
  
  recovered <- sapply(names(init), \(name) links[[name]](latent[[name]]))
  
  cat("\n", kernel, "\n", sep = "")
  print(rbind(initial = init, recovered = recovered))
  
  stopifnot(max(abs(init - recovered)) < 1e-8)
}