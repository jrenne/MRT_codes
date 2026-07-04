# ==============================================================================
# MRT auxiliary code: euro-area application
# ==============================================================================
# This script reproduces the EA workflow present in the original project folder.
# It is included as auxiliary material; main.R is the US replication used by the
# current paper.

rm(list = ls())

message("MRT replication: starting auxiliary EA workflow.")

# ------------------------------------------------------------------------------
# User switches
# ------------------------------------------------------------------------------

# Choose how much of the pipeline to run:
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

# Set to TRUE to rebuild Gaussian-mixture survey inputs from the raw ECB SPF
# files. The default uses the stored workbook inputs shipped with this folder.
indic.compute.mixture <- FALSE

# The EA data-loading scripts use Eurostat for HICP data. Internet access may
# therefore be needed even when the cached model inputs are used.

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

library(tis)
library(Rcpp)

message("MRT replication: creating EA output folders.")
required_dirs <- c(
  "graphs/EA_2024/Baseline",
  "tables/EA/Baseline"
)
invisible(lapply(required_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

message("MRT replication: loading shared R and C++ routines.")
source("RScripts/0-procedure.R")
source("RScripts/1-procedure.R")
Rcpp::sourceCpp("RScripts/function_loadings.cpp")
Rcpp::sourceCpp("RScripts/kalman_cpp.cpp")

if (indic.compute.mixture) {
  message("MRT replication: loading raw EA data and recomputing mixture inputs.")
  source("RScripts/EA/load.data.EA_new_new.R")
  source("RScripts/EA/load.data.GDP.EA_new_new.R")
  source("RScripts/EA/mixture-distribution_estimate_EA_calendar.R")
  source("RScripts/EA/create_data_for_KF_with_k4_EA_calendar.R")
  source("RScripts/EA/mixture-distribution_GDP_estimate_EA_calendar.R")
  source("RScripts/EA/create_data_for_KF_with_k4_GDP_EA_calendar.R")
} else {
  message("MRT replication: loading local EA macro series needed for output figures.")
  hcpi_raw <- read.table("data/EA/raw/prc_hicp_midx_page_linear.csv.gz",
                         sep = ",", quote = "\"", header = TRUE)
  hcpi <- hcpi_raw %>%
    filter(geo == "EA", coicop == "CP00", unit == "I15") %>%
    mutate(date = as.Date(paste0(TIME_PERIOD, "-01")) + 14) %>%
    dplyr::select(date, OBS_VALUE) %>%
    rename(EA.hcpi = OBS_VALUE)
  ts_hcpi <- ts(hcpi$EA.hcpi, frequency = 12, start = c(1996, 1))
  stl_hcpi <- stl(ts_hcpi, s.window = "periodic")
  hcpi$EA.hcpi.deseasonalized <- rowSums(stl_hcpi$time.series[, c("trend", "remainder")])
  DATA <- as.data.frame(hcpi)

  gdp_raw <- read.table("data/EA/raw/namq_10_gdp_page_linear.csv.gz",
                        sep = ",", quote = "\"", header = TRUE)
  gdp <- gdp_raw %>%
    filter(geo == "EA") %>%
    mutate(date = as.Date(as.yearqtr(paste(substr(TIME_PERIOD, 1, 4), substr(TIME_PERIOD, 6, 7)))) + 14) %>%
    dplyr::select(date, OBS_VALUE) %>%
    rename(EA.gdp = OBS_VALUE)
  ts_gdp <- ts(gdp$EA.gdp, frequency = 4, start = c(1995, 1))
  stl_gdp <- stl(ts_gdp, s.window = "periodic")
  gdp$EA.gdp.deseasonalized <- rowSums(stl_gdp$time.series[, c("trend", "remainder")])
  DATA.G <- as.data.frame(gdp)
}

source_script <- function(file, env = parent.frame(), echo = FALSE) {
  message("MRT replication: running ", file)
  source(file, local = env, echo = echo)
}

# ------------------------------------------------------------------------------
# Baseline EA specification
# ------------------------------------------------------------------------------

indic.4th.use <- FALSE
indic.3rd.use <- TRUE
indic.model.var.only <- FALSE
indic.make.implied.distri <- "FALSE"
indic.estimate <- run_mode != "outputs_from_saved"
estimation_start <- if (run_mode == "estimate_from_generic") "generic" else "saved"
save_estimation_results <- indic.estimate

run_env <- environment()
source_script("RScripts/0-procedure.R", env = run_env)
source_script("RScripts/1-procedure.R", env = run_env)
source_script("RScripts/EA/Estimate_joint_EA_model_no_higher_order_calendar_short.R", env = run_env)
source_script("RScripts/EA/model.implied.distribution.joint.EA.R", env = run_env)

path_graph <- "graphs/EA_2024/Baseline/"
source_script("RScripts/make_outputs/make_figure_factors.R", env = run_env)
source_script("RScripts/make_outputs/make_figure_correl.R", env = run_env)
source_script("RScripts/make_outputs/make_figure_fit_Surveys.R", env = run_env)
if (indic.compute.mixture) {
  source_script("RScripts/make_outputs/make_figure_model_implied_distri.R", env = run_env)
} else {
  message("MRT replication: skipping raw SPF distribution plots in cached EA mode.")
}
source_script("RScripts/make_outputs/make_figure_TrendCycle.R", env = run_env)

# This script needs echo=TRUE because some ggplot objects are printed only when
# the sourced expressions are echoed.
source_script("RScripts/make_outputs/make_figure_DemandSupply.R", env = run_env, echo = TRUE)

path_table <- "tables/EA/Baseline/"
source_script("RScripts/make_outputs/make_table_param.R", env = run_env)
source_script("RScripts/make_outputs/make_workbooks.R", env = run_env)

message("MRT replication: auxiliary EA workflow complete.")
