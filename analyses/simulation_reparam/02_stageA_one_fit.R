#===============================================================================
# Stage A: first reparameterised convergence test
#===============================================================================

#-------------------------------------------------------------------------------
# Load Stage-A design
#-------------------------------------------------------------------------------

library(here)
source(here("analyses", "simulation_reparam", "01_stageA_design.R"))

#-------------------------------------------------------------------------------
# Use one numerical thread
#-------------------------------------------------------------------------------

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
INLA::inla.setOption(num.threads = "1:1")

#-------------------------------------------------------------------------------
# Read copied production results from reparametrisation worktree
#-------------------------------------------------------------------------------

production_main_fit_dir <- normalizePath(
  file.path(here(), "..",
    "seismology-project", "results", "simulation", "fits", "main_simulation"),
  mustWork = TRUE)

problem_fits_file <- file.path(production_main_fit_dir, "problem_fits.csv")
simulation_manifest_file <- file.path(catalogue_dir, "simulation_manifest.csv")

stopifnot(file.exists(problem_fits_file), file.exists(simulation_manifest_file))

problem_fits <- read.csv(problem_fits_file, stringsAsFactors = FALSE)

simulation_manifest <- read.csv(simulation_manifest_file, stringsAsFactors = FALSE)

#-------------------------------------------------------------------------------
# Pick a known troublesome fit
#-------------------------------------------------------------------------------

candidates <- problem_fits[problem_fits$truth_kernel == "mse" &
                             problem_fits$fitted_kernel == "rate_state" &
                             problem_fits$hit_max %in% TRUE, ]

stopifnot(nrow(candidates) > 0)

job <- candidates[1, ]

truth_i <- job$truth_kernel
fitted_i <- job$fitted_kernel
rep_i <- job$rep
catalogue_seed_i <- job$catalogue_seed
fit_seed_i <- job$fit_seed

cat("Stage A test:\n", "truth =", truth_i, "| rep =", rep_i, 
    "| fitted =", fitted_i, "| original iterations =", job$n_iter, "\n\n")

#-------------------------------------------------------------------------------
# Load exact production catalogue
#-------------------------------------------------------------------------------

cat_row <- simulation_manifest[
  simulation_manifest$kernel == truth_i & simulation_manifest$rep == rep_i, ]

stopifnot(nrow(cat_row) == 1)

obj_cat <- readRDS(
  file.path(catalogue_dir, cat_row$file))

stopifnot(obj_cat$truth_kernel == truth_i, obj_cat$rep == rep_i, 
           obj_cat$seed == catalogue_seed_i)

catalogue_i <- obj_cat$catalogue
catalogue_i <- catalogue_i[
  catalogue_i$ts >= T_fit_start &
    catalogue_i$ts <= T_fit_end, ]

catalogue_i <- catalogue_i[order(catalogue_i$ts), ]
catalogue_i$idx.p <- seq_len(nrow(catalogue_i))

#-------------------------------------------------------------------------------
# Reparameterised fit
#-------------------------------------------------------------------------------

set.seed(fit_seed_i)

link_i <- make_links_stageA(fitted_i)

bru_i <- make_bru_options_stageA(fitted_i)

bru_i$control.compute <- list(
  config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE, mlik = TRUE)

start_i <- Sys.time()

fit_i <- ETAS.inlabru::Temporal.ETAS(
  total.data = catalogue_i,
  M0 = M0,
  T1 = T_fit_start,
  T2 = T_fit_end,
  link.functions = link_i,
  coef.t. = temporal_binning$coef.t,
  delta.t. = temporal_binning$delta.t,
  N.max. = temporal_binning$N.max,
  bru.opt = bru_i,
  kernel = fitted_i)

runtime_i <- as.numeric(difftime(Sys.time(), start_i, units = "mins"))

#-------------------------------------------------------------------------------
# Convergence
#-------------------------------------------------------------------------------

log_i <- as.character(inlabru::bru_log(fit_i))

inla_failure_i <- any(grepl(
  paste(
    "Problem in inla",
    "Giving up and returning last successfully obtained result",
    "inla-program exited with an error",
    "maximum number of tries has been reached",
    "Newton-Raphson optimizer did not converge",
    sep = "|"), log_i))

converged_i <- !inla_failure_i && any(
  grepl("Convergence criterion met", log_i, fixed = TRUE))

hit_max_i <- !inla_failure_i && any(
  grepl("Maximum iterations reached", log_i, fixed = TRUE))

n_iter_i <- max(fit_i$bru_iinla$track$iteration, na.rm = TRUE)

cat("\nStage A result:\n", "converged =", converged_i, "\n",
    "hit_max =", hit_max_i, "\n", "iterations =", n_iter_i, "\n",
    "inla_failure =", inla_failure_i, "\n",
    "runtime minutes =", round(runtime_i, 2), "\n")

#-------------------------------------------------------------------------------
# Save Stage-A results
#-------------------------------------------------------------------------------

stageA_dir <- here("results", "simulation_reparam", "stageA")

dir.create(stageA_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  list(
    fit = fit_i,
    truth_kernel = truth_i,
    fitted_kernel = fitted_i,
    rep = rep_i,
    fit_seed = fit_seed_i,
    catalogue_seed = catalogue_seed_i,
    original_n_iter = job$n_iter,
    converged = converged_i,
    hit_max = hit_max_i,
    inla_failure = inla_failure_i,
    n_iter = n_iter_i,
    runtime_minutes = runtime_i),
  file.path(
    stageA_dir,sprintf( "truth_%s_rep_%04d_fit_%s.rds", truth_i, rep_i, fitted_i))
  )