#===============================================================================
# Agreement of (D,S) parametrisation with native rate-state implementation
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
# Moderate parameter values
#-------------------------------------------------------------------------------

D_values <- c(0.1, 1, 5, 60)
S_values <- c(0.1, 1, 5, 7, 10, 12)

tt <- c(0, 0.001, 0.1, 1, 10, 100, 2371)
aa <- c(0, 0, 0.1, 1, 0)
bb <- c(0.1, 1, 1, 10, 2371)

kernel_error <- numeric()
integral_error <- numeric()

for (D in D_values) {
  for (S in S_values) {
    
    native <- ETAS.inlabru:::temporal_reparam_to_native(
      list(decay_rs = D, shape_rs = S), "rate_state")
    
    old_g <- temporal_kernel(tt, native, "rate_state")
    new_g <- ETAS.inlabru:::rate_state_reparam_kernel(tt, D, S)
    
    old_I <- temporal_kernel_integral(aa, bb, native, "rate_state")
    new_I <- ETAS.inlabru:::rate_state_reparam_integral(aa, bb, D, S)
    
    kernel_error <- c(kernel_error, abs(old_g - new_g))
    
    integral_error <- c(integral_error, abs(old_I - new_I))
  }
}

cat("max kernel difference:", max(kernel_error), "\n")

cat("max integral difference:", max(integral_error), "\n")

stopifnot(max(kernel_error) < 1e-7, max(integral_error) < 1e-7)

cat("Native-equivalence check PASSED\n")

#-------------------------------------------------------------------------------
# Density and integral calculations for extreme parameter values
#-------------------------------------------------------------------------------

D_extreme <- c(1e-300, 1e-100, 1e-20, 1e-6, 0.1, 1, 60, 1e6, 1e20, 1e100, 1e300)

S_extreme <- c(1e-300, 1e-20, 1e-8, 1e-3, 0.1, 1, 10, 50, 100, 500, 700, 1e6)

tt_extreme <- c(0, 1e-12, 1e-6, 0.001, 0.1, 1, 100, 2371)

aa_extreme <- c(0, 0, 0.001, 1, 100)
bb_extreme <- c(1e-6, 1, 1, 100, 2371)

for (D in D_extreme) {
  for (S in S_extreme) {
    
    g <- ETAS.inlabru:::rate_state_reparam_kernel(tt_extreme, D, S)
    
    I <- ETAS.inlabru:::rate_state_reparam_integral(aa_extreme, bb_extreme, D, S)
    
    width <- bb_extreme - aa_extreme
    
    stopifnot(
      all(is.finite(g)),
      all(g >= 0),
      all(g <= 1 + 1e-12),
      all(is.finite(I)),
      all(I >= 0),
      all(I <= width + 1e-10 * pmax(1, width)))
  }
}

cat("Extreme parameter check PASSED\n")

#-------------------------------------------------------------------------------
# Perturb decay_rs around the actual starting value (tests for discontinuities)
#-------------------------------------------------------------------------------

D0 <- 5.01
S0 <- 6.216606

eps <- c(-1e-5, -1e-7, 0, 1e-7, 1e-5)

for (e in eps) {
  
  D <- D0 * exp(e)
  
  g <- ETAS.inlabru:::rate_state_reparam_kernel(tt_extreme, D, S0)
  
  I <- ETAS.inlabru:::rate_state_reparam_integral(aa_extreme, bb_extreme, D, S0)
  
  stopifnot(all(is.finite(g)), all(is.finite(I)))
}

cat("Local decay_rs perturbation check PASSED\n")
cat("\nALL DIRECT RATE-STATE CHECKS PASSED\n")