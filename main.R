# ==============================================================================
# MRT replication code: US application
# ==============================================================================
# Run this file from the root of the replication folder.
# The script uses stored intermediate objects by default, so a standard run
# reproduces figures/tables without recomputing the slow estimation steps.

rm(list = ls())

message("MRT replication: starting US workflow.")

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

# Set to TRUE to re-estimate the Gaussian-mixture approximations to the SPF
# histograms and rebuild the Kalman-filter observables. The default uses the
# stored workbook inputs shipped with the replication package.
indic.compute.mixture <- FALSE

# Set a FRED key with Sys.setenv(FRED_API_KEY = "...") before running this file
# to reproduce the term-premium and output-gap comparison overlays exactly.

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

library(tis)      # NBER recession dates used in the figures
library(Rcpp)     # sourceCpp for the C++ routines

message("MRT replication: creating output folders.")
required_dirs <- c(
  "graphs/US_2024",
  "graphs/US_2024/Baseline",
  "graphs/US_2024/All_moments",
  "graphs/US_2024/No_3rd_4th",
  "tables",
  "tables/No_3rd_4th",
  "tables/No_4th"
)
invisible(lapply(required_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Shared R and C++ routines.
message("MRT replication: loading shared R and C++ routines.")
source("RScripts/0-procedure.R")
source("RScripts/1-procedure.R")
Rcpp::sourceCpp("RScripts/function_loadings.cpp")
Rcpp::sourceCpp("RScripts/kalman_cpp.cpp")

# Raw and pre-processed US data used throughout the replication.
message("MRT replication: loading US raw and processed data.")
source("RScripts/US/load.data.US.R")
source("RScripts/US/load.data.US.G.R")

if (indic.compute.mixture) {
  message("MRT replication: recomputing SPF mixture inputs.")
  # Rebuild the survey-distribution inputs from the raw SPF histograms.
  source("RScripts/US/mixture-distribution_estimate_US.R")
  source("RScripts/US/create_data_for_KF_with_k4_US.R")
  source("RScripts/US/mixture-distribution_estimate_US_G.R")
  source("RScripts/US/create_data_for_KF_with_k4_US_G.R")
}

# ------------------------------------------------------------------------------
# Model/output runner
# ------------------------------------------------------------------------------

source_script <- function(file, env = parent.frame(), echo = FALSE) {
  message("MRT replication: running ", file)
  source(file, local = env, echo = echo)
}

source_outputs <- function(files, env = parent.frame()) {
  invisible(lapply(files, source_script, env = env))
}

run_model_specification <- function(path_graph,
                                    use_skewness,
                                    use_kurtosis,
                                    variance_only,
                                    compute_implied_distribution,
                                    include_full_output_set = FALSE,
                                    write_parameter_table = FALSE) {
  message("MRT replication: starting model specification -> ", path_graph)
  # These switches are read by the estimation and output scripts.
  indic.3rd.use <- use_skewness
  indic.4th.use <- use_kurtosis
  indic.model.var.only <- variance_only
  indic.make.implied.distri <- if (compute_implied_distribution) "TRUE" else "FALSE"
  indic.estimate <- run_mode != "outputs_from_saved"
  estimation_start <- if (run_mode == "estimate_from_generic") "generic" else "saved"
  save_estimation_results <- indic.estimate

  run_env <- environment()
  # The model helpers use the specification-specific horizon objects (`H`,
  # `select.inflation.types`, `freq`) through lexical scoping.
  source_script("RScripts/0-procedure.R", env = run_env)
  source_script("RScripts/1-procedure.R", env = run_env)
  source_script("RScripts/US/Estimate_joint_model_4_5_quarterly.US.trend.cycle.quartely.infl.gdp.short.R", env = run_env)
  source_script("RScripts/US/model.implied.distribution.joint.R", env = run_env)

  core_figures <- c(
    "RScripts/make_outputs/make_figure_factors.R",
    "RScripts/make_outputs/make_figure_correl.R",
    "RScripts/make_outputs/make_figure_fit_Surveys.R",
    "RScripts/make_outputs/make_figure_model_implied_distri.R",
    "RScripts/make_outputs/make_figure_TrendCycle.R"
  )
  source_outputs(core_figures, env = run_env)

  if (include_full_output_set) {
    extra_figures <- c(
      "RScripts/make_outputs/make_figure_motiv.R",
      "RScripts/make_outputs/make_figure_fit_Surveys_matrix1Horizon.R",
      "RScripts/make_outputs/make_figure_fit_GDP_Infl.R",
      "RScripts/make_outputs/make_figure_TrendCycle.R"
    )
    source_outputs(extra_figures, env = run_env)
  }

  # These scripts need echo=TRUE because some ggplot objects are printed only
  # when the sourced expressions are echoed.
  source_script("RScripts/make_outputs/make_figure_DemandSupply.R", env = run_env, echo = TRUE)
  if (include_full_output_set) {
    source_script("RScripts/make_outputs/make_figure_DemandSupply_var_k3rd.R", env = run_env, echo = TRUE)
  }

  if (write_parameter_table) {
    source_script("RScripts/make_outputs/make_table_param.R", env = run_env)
  }

  source_script("RScripts/make_outputs/make_workbooks.R", env = run_env)
  message("MRT replication: completed model specification -> ", path_graph)
}

# ------------------------------------------------------------------------------
# Replication runs
# ------------------------------------------------------------------------------

# Comparison specification: variance, skewness, and kurtosis observables.
run_model_specification(
  path_graph = "graphs/US_2024/All_moments/",
  use_skewness = TRUE,
  use_kurtosis = TRUE,
  variance_only = FALSE,
  compute_implied_distribution = FALSE,
  include_full_output_set = TRUE,
  write_parameter_table = FALSE
)

# Robustness: variance observables only.
run_model_specification(
  path_graph = "graphs/US_2024/No_3rd_4th/",
  use_skewness = FALSE,
  use_kurtosis = FALSE,
  variance_only = TRUE,
  compute_implied_distribution = TRUE
)

# Baseline used in the paper: variance and skewness observables.
run_model_specification(
  path_graph = "graphs/US_2024/Baseline/",
  use_skewness = TRUE,
  use_kurtosis = FALSE,
  variance_only = FALSE,
  compute_implied_distribution = TRUE,
  write_parameter_table = TRUE
)

# Some paper outputs are shipped as static files in graphs/US_2024/Baseline/
# and tables/ when the active code archive does not regenerate them directly.
# RScripts/check_paper_outputs.R checks that all paper inputs are covered.

message("MRT replication: US workflow complete.")
