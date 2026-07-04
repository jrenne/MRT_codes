# ==============================================================================
# MRT replication code: US application
# ==============================================================================
# Run this file from the root of the replication folder.
# The script uses stored intermediate objects by default, so a standard run
# reproduces figures/tables without recomputing the slow estimation steps.

rm(list = ls())

message("MRT replication: starting US workflow.")

quiet_source <- function(file, env = parent.frame(), echo = FALSE) {
  suppressPackageStartupMessages(
    invisible(capture.output(source(file, local = env, echo = echo)))
  )
}

# ------------------------------------------------------------------------------
# User switches
# ------------------------------------------------------------------------------

# Choose how much of the pipeline to run. The default reproduces tables and
# figures from saved estimates. Re-estimation is intentionally opt-in.
# - "outputs_from_saved": use stored model estimates and regenerate outputs.
# - "estimate_from_saved": re-estimate, initialized from stored estimates.
# - "estimate_from_generic": re-estimate, initialized from simple generic values.
run_mode <- "outputs_from_saved"

allowed_run_modes <- c(
  "outputs_from_saved",
  "estimate_from_saved",
  "estimate_from_generic"
)
if (!run_mode %in% allowed_run_modes) {
  stop("run_mode must be one of: ", paste(allowed_run_modes, collapse = ", "))
}

# Choose which US model specifications to process. "Process" means the full
# workflow for that specification: load or estimate the model according to
# `run_mode`, then produce its figures, tables, and workbooks.
# Use all three for the full paper replication, or a subset for quicker checks.
# Allowed names:
# - "baseline": variance and skewness observables; the paper's main model.
# - "all_moments": variance, skewness, and kurtosis observables.
# - "no_3rd_4th": variance observables only.
# Example: model_specifications_to_process <- "baseline"
model_specifications_to_process <- c("baseline",
                                     "all_moments",
                                     "no_3rd_4th")

# Model files used by the three US specifications. Re-estimation starts from
# `initial_model_files`, stores new estimates in `estimation_save_model_files`,
# and `outputs_from_saved` reads from `chart_model_files`.
canonical_model_files <- c(
  all_moments = "results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.4th.errors.no.4Q.best.final.new.RData",
  no_3rd_4th  = "results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.no.3rd.4th.no.4Q.RData",
  baseline    = "results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.only.errors.no.4Q.best.RData"
)
estimation_save_model_files <- c(
  all_moments = "results/US/estimation_runs/US.Model.all.moments.reestimated.RData",
  no_3rd_4th  = "results/US/estimation_runs/US.Model.no.3rd.4th.reestimated.RData",
  baseline    = "results/US/estimation_runs/US.Model.baseline.reestimated.RData"
)

initial_model_files <- canonical_model_files # starts from there if new estimation
chart_model_files   <- canonical_model_files # these files are the ones used to produce charts


# Set to TRUE to recompute the smoothed survey-distribution inputs from the raw
# SPF histograms before building the Kalman-filter observables. The default uses
# the precomputed workbook inputs included in data/processed/.
recompute_mixture_inputs <- FALSE
indic.compute.mixture <- recompute_mixture_inputs

# Set to TRUE only when you want the long diagnostic PDFs from the smoothing
# step. These are not needed to reproduce the paper figures and tables.
save_mixture_diagnostic_plots <- FALSE

# Show progress messages for the Monte Carlo simulation used in the
# stagflation-probability figure.
show_stagflation_progress <- TRUE

# Stagflation probabilities are based on this many Monte Carlo draws for each
# date. Use 10,000 for the final paper replication; 2,000 is faster while
# checking the workflow.
stagflation_nb_sim <- 2000
stagflation_batch_size <- 500

# Show progress messages when loading or recomputing model-implied
# distributions.
show_implied_distribution_progress <- TRUE

# Optimization controls used only when run_mode is "estimate_from_saved" or
# "estimate_from_generic". Each outer loop runs:
#   Nelder-Mead -> nlminb -> Nelder-Mead.
# These defaults reproduce the historical estimation schedule.
optimization_setup <- list(
  outer_loops = 3,
  nelder_mead_first_maxit = 2000,
  nlminb_maxit = 30,
  nelder_mead_second_maxit = 2000,
  trace = TRUE,
  compute_hessian = FALSE,
  nlminb_kkt = FALSE
)

# Set a FRED key with Sys.setenv(FRED_API_KEY = "...") before running this file
# to reproduce the term-premium and output-gap comparison overlays exactly.

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tis)      # NBER recession dates used in the figures
  library(Rcpp)     # sourceCpp for the C++ routines
})

message("MRT replication: creating output folders.")
required_dirs <- c(
  "graphs/US_2024",
  "graphs/US_2024/Baseline",
  "graphs/US_2024/All_moments",
  "graphs/US_2024/No_3rd_4th",
  "results/US/estimation_runs",
  "tables",
  "tables/US_2024/Baseline",
  "tables/US_2024/All_moments",
  "tables/US_2024/No_3rd_4th"
)
invisible(lapply(required_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Shared R and C++ routines.
message("MRT replication: loading shared R and C++ routines.")
quiet_source("RScripts/helpers/0-procedure.R")
quiet_source("RScripts/helpers/1-procedure.R")
Rcpp::sourceCpp("RScripts/helpers/function_loadings.cpp")
Rcpp::sourceCpp("RScripts/helpers/kalman_cpp.cpp")

# Raw and pre-processed US data used throughout the replication.
message("MRT replication: loading US raw and processed data.")
.data_load_plot_device <- grDevices::dev.cur()
grDevices::pdf(NULL)
.data_load_plot_device <- grDevices::dev.cur()
quiet_source("RScripts/US/load.data.US.R")
quiet_source("RScripts/US/load.data.US.G.R")
if (.data_load_plot_device %in% grDevices::dev.list()) {
  grDevices::dev.off(.data_load_plot_device)
}
unlink("Rplots.pdf")

if (indic.compute.mixture) {
  message("MRT replication: recomputing smoothed SPF distribution inputs.")
  message("MRT replication: this step fits the survey distributions and rebuilds cumulant-target observables.")
  quiet_source("RScripts/US/mixture-distribution_estimate_US.R")
  quiet_source("RScripts/US/create_data_for_KF_with_k4_US.R")
  quiet_source("RScripts/US/mixture-distribution_estimate_US_G.R")
  quiet_source("RScripts/US/create_data_for_KF_with_k4_US_G.R")
} else {
  message("MRT replication: using precomputed smoothed SPF distribution inputs from data/processed/.")
}

# ------------------------------------------------------------------------------
# Model/output runner
# ------------------------------------------------------------------------------

source_script <- function(file, env = parent.frame(), echo = FALSE, quiet = TRUE) {
  message("MRT replication: running ", file)
  if (quiet) {
    quiet_source(file, env = env, echo = echo)
  } else {
    source(file, local = env, echo = echo)
  }
}

source_outputs <- function(files, env = parent.frame()) {
  invisible(lapply(files, source_script, env = env))
}

run_model_specification <- function(path_graph,
                                    specification,
                                    use_skewness,
                                    use_kurtosis,
                                    variance_only,
                                    compute_implied_distribution,
                                    implied_distribution_file,
                                    include_full_output_set = FALSE,
                                    write_parameter_table = FALSE,
                                    write_distribution_tables = FALSE,
                                    write_stagflation_figure = FALSE,
                                    path_table = NULL) {
  message("")
  message("------------------------------------------------------------")
  message("MRT replication: starting model specification -> ", specification)
  message("------------------------------------------------------------")
  # These switches are read by the estimation and output scripts.
  indic.3rd.use <- use_skewness
  indic.4th.use <- use_kurtosis
  indic.model.var.only <- variance_only
  indic.make.implied.distri <- if (compute_implied_distribution) "TRUE" else "FALSE"
  path_implied_distribution <- implied_distribution_file
  show_stagflation_progress <- show_stagflation_progress
  stagflation_nb_sim <- stagflation_nb_sim
  stagflation_batch_size <- stagflation_batch_size
  show_implied_distribution_progress <- show_implied_distribution_progress
  indic.estimate <- run_mode != "outputs_from_saved"
  estimation_start <- if (run_mode == "estimate_from_generic") "generic" else "saved"
  model_files <- list(
    initial = unname(initial_model_files[[specification]]),
    output = unname(chart_model_files[[specification]]),
    save = unname(estimation_save_model_files[[specification]])
  )
  save_estimation_results <- indic.estimate &&
    !is.null(model_files$save) &&
    !is.na(model_files$save) &&
    nzchar(model_files$save)
  optimization_setup <- optimization_setup

  run_env <- environment()
  # The model helpers use the specification-specific horizon objects (`H`,
  # `select.inflation.types`, `freq`) through lexical scoping.
  source_script("RScripts/helpers/0-procedure.R", env = run_env)
  source_script("RScripts/helpers/1-procedure.R", env = run_env)
  source_script(
    "RScripts/US/Estimate_joint_model_4_5_quarterly.US.trend.cycle.quartely.infl.gdp.short.R",
    env = run_env,
    quiet = !indic.estimate
  )
  source_script("RScripts/US/model.implied.distribution.joint.R", env = run_env)

  core_figures <- c(
    "RScripts/make_figures/make_figure_factors.R",
    "RScripts/make_figures/make_figure_correl.R",
    "RScripts/make_figures/make_figure_fit_Surveys.R",
    "RScripts/make_figures/make_figure_model_implied_distri.R",
    "RScripts/make_figures/make_figure_TrendCycle.R"
  )
  source_outputs(core_figures, env = run_env)

  if (include_full_output_set) {
    extra_figures <- c(
      "RScripts/make_figures/make_figure_motiv.R",
      "RScripts/make_figures/make_figure_fit_Surveys_matrix1Horizon.R",
      "RScripts/make_figures/make_figure_fit_GDP_Infl.R",
      "RScripts/make_figures/make_figure_TrendCycle.R"
    )
    source_outputs(extra_figures, env = run_env)
  }

  # These scripts need echo=TRUE because some ggplot objects are printed only
  # when the sourced expressions are echoed.
  source_script("RScripts/make_figures/make_figure_DemandSupply.R", env = run_env, echo = TRUE)
  if (include_full_output_set) {
    source_script("RScripts/make_figures/make_figure_DemandSupply_var_k3rd.R", env = run_env, echo = TRUE)
  }
  if (write_stagflation_figure) {
    source_script("RScripts/make_figures/make_figure_stagflation_probability.R", env = run_env)
  }

  if (write_parameter_table) {
    if (is.null(path_table)) {
      stop("Set path_table when writing tables.", call. = FALSE)
    }
    source_script("RScripts/make_tables/make_table_param.R", env = run_env)
  }
  if (write_distribution_tables) {
    if (is.null(path_table)) {
      stop("Set path_table when writing tables.", call. = FALSE)
    }
    source_script("RScripts/make_tables/make_tables_distribution_diagnostics.R", env = run_env)
  }

  source_script("RScripts/make_workbooks/make_workbooks.R", env = run_env)
  message("MRT replication: completed model specification -> ", specification)
}

run_if_selected <- function(specification, ...) {
  if (specification %in% model_specifications_to_process) {
    run_model_specification(specification = specification, ...)
  } else {
    message("MRT replication: skipping model specification -> ", specification)
  }
}

# ------------------------------------------------------------------------------
# Replication runs
# ------------------------------------------------------------------------------

# Comparison specification: variance, skewness, and kurtosis observables.
run_if_selected(
  path_graph = "graphs/US_2024/All_moments/",
  specification = "all_moments",
  use_skewness = TRUE,
  use_kurtosis = TRUE,
  variance_only = FALSE,
  compute_implied_distribution = TRUE,
  implied_distribution_file = "results/US/US.Model.Implied.Distribution.all.moments.RData",
  include_full_output_set = TRUE,
  write_parameter_table = FALSE
)

# Robustness: variance observables only.
run_if_selected(
  path_graph = "graphs/US_2024/No_3rd_4th/",
  specification = "no_3rd_4th",
  use_skewness = FALSE,
  use_kurtosis = FALSE,
  variance_only = TRUE,
  compute_implied_distribution = TRUE,
  implied_distribution_file = "results/US/US.Model.Implied.Distribution.no.3rd.4th.RData"
)

# Baseline used in the paper: variance and skewness observables.
run_if_selected(
  path_graph = "graphs/US_2024/Baseline/",
  specification = "baseline",
  use_skewness = TRUE,
  use_kurtosis = FALSE,
  variance_only = FALSE,
  compute_implied_distribution = TRUE,
  implied_distribution_file = "results/US/US.Model.Implied.Distribution.with.3rd.only.RData",
  write_parameter_table = TRUE,
  write_distribution_tables = TRUE,
  write_stagflation_figure = TRUE,
  path_table = "tables/US_2024/Baseline/"
)

# RScripts/check_paper_outputs.R checks that all paper inputs are regenerated
# by the replication workflow.

message("MRT replication: US workflow complete.")
