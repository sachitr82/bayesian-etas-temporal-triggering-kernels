#===============================================================================
# Stage A: small convergence benchmark
#===============================================================================

#-------------------------------------------------------------------------------
# Load Stage-A design
#-------------------------------------------------------------------------------

library(dplyr)
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
# Select five difficult fits in each main cross-kernel direction
#-------------------------------------------------------------------------------

mse_to_rs <- problem_fits[
  problem_fits$truth_kernel == "mse" &
    problem_fits$fitted_kernel == "rate_state"&
    problem_fits$hit_max %in% TRUE, ]

rs_to_mse <- problem_fits[
  problem_fits$truth_kernel == "rate_state" &
    problem_fits$fitted_kernel == "mse" &
    problem_fits$hit_max %in% TRUE,]

mse_to_rs <- mse_to_rs[order(mse_to_rs$rep), ][1:5, ]
rs_to_mse <- rs_to_mse[order(rs_to_mse$rep), ][1:5, ]

benchmark_jobs <- rbind(mse_to_rs, rs_to_mse)
rownames(benchmark_jobs) <- NULL

stopifnot(nrow(benchmark_jobs) == 10, 
          all(benchmark_jobs$n_iter == fit_control$max_iter))

#-------------------------------------------------------------------------------
# Output directory
#-------------------------------------------------------------------------------

stageA_dir <- here("results", "simulation_reparam", "stageA")

dir.create(stageA_dir, recursive = TRUE, showWarnings = FALSE)

benchmark_rows <- vector("list", nrow(benchmark_jobs))

#-------------------------------------------------------------------------------
# Run benchmark
#-------------------------------------------------------------------------------

for (i in seq_len(nrow(benchmark_jobs))) {
  
  job <- benchmark_jobs[i, ]
  
  truth_i <- job$truth_kernel
  fitted_i <- job$fitted_kernel
  rep_i <- job$rep
  catalogue_seed_i <- job$catalogue_seed
  fit_seed_i <- job$fit_seed
  
  outfile_i <- file.path(stageA_dir,
    sprintf("truth_%s_rep_%04d_fit_%s.rds", truth_i, rep_i, fitted_i))
  
  cat( "\n", i, "/10 | truth =", truth_i, "| rep =", rep_i,
       "| fit =", fitted_i, "\n")
  
  # Reuse completed fits
  if (file.exists(outfile_i)) {
    
    obj_i <- readRDS(outfile_i)
    
    stopifnot(
      obj_i$truth_kernel == truth_i, obj_i$fitted_kernel == fitted_i, 
      obj_i$rep == rep_i, obj_i$fit_seed == fit_seed_i)
    
    log_i <- as.character(inlabru::bru_log(obj_i$fit))
    
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
    
    fit_status_i <- case_when(inla_failure_i ~ "inla_failure", converged_i ~ "converged",
                         hit_max_i ~ "hit_max", TRUE ~ "unknown")
    
    benchmark_rows[[i]] <- data.frame(
      truth_kernel = truth_i,
      fitted_kernel = fitted_i,
      rep = rep_i,
      original_n_iter = job$n_iter,
      reparam_n_iter = obj_i$n_iter,
      original_runtime = job$runtime_minutes,
      reparam_runtime = obj_i$runtime_minutes,
      fit_status = fit_status_i,
      converged = converged_i,
      hit_max = hit_max_i,
      inla_failure = inla_failure_i)
    
    cat("Existing Stage-A fit reused\n")
    next
  }
  
  # Load exact production catalogue
  cat_row <- simulation_manifest[
    simulation_manifest$kernel == truth_i &
      simulation_manifest$rep == rep_i, ]
  
  stopifnot(nrow(cat_row) == 1)
  
  obj_cat <- readRDS(file.path(catalogue_dir, cat_row$file))
  
  stopifnot(obj_cat$truth_kernel == truth_i, obj_cat$rep == rep_i, 
            obj_cat$seed == catalogue_seed_i)
  
  catalogue_i <- obj_cat$catalogue
  catalogue_i <- catalogue_i[
    catalogue_i$ts >= T_fit_start &
      catalogue_i$ts <= T_fit_end, ]
  
  catalogue_i <- catalogue_i[order(catalogue_i$ts), ]
  catalogue_i$idx.p <- seq_len(nrow(catalogue_i))
  
  # Reparameterised fit
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
  
  fit_status_i <- case_when(inla_failure_i ~ "inla_failure", converged_i ~ "converged",
                        hit_max_i ~ "hit_max", TRUE ~ "unknown")
  n_iter_i <- max(fit_i$bru_iinla$track$iteration, na.rm = TRUE)
  
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
      fit_status = fit_status_i,
      n_iter = n_iter_i,
      runtime_minutes = runtime_i
    ), outfile_i)
  
  benchmark_rows[[i]] <- data.frame(
    truth_kernel = truth_i,
    fitted_kernel = fitted_i,
    rep = rep_i,
    original_n_iter = job$n_iter,
    reparam_n_iter = n_iter_i,
    original_runtime = job$runtime_minutes,
    reparam_runtime = runtime_i,
    fit_status = fit_status_i,
    converged = converged_i,
    hit_max = hit_max_i,
    inla_failure = inla_failure_i)
  
  cat("fit status =", fit_status_i, "| iterations =", n_iter_i,
      "| runtime =", round(runtime_i, 2), "min\n")
}


#-------------------------------------------------------------------------------
# Benchmark results
#-------------------------------------------------------------------------------

benchmark_results <- do.call(rbind, benchmark_rows)
rownames(benchmark_results) <- NULL

cat("\nStage-A benchmark results\n\n")
print(benchmark_results, row.names = FALSE)

write.csv(benchmark_results, file.path(stageA_dir, "benchmark_10fits.csv"),
          row.names = FALSE)