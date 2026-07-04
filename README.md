# Codes for "Inflation and Growth Risk: Balancing the Scales with Surveys"

This repository contains the codes and data needed to reproduce the empirical results of the paper:

**Inflation and Growth Risk: Balancing the Scales with Surveys**  
Sarah Mouabbi, Jean-Paul Renne, and Adrien Tschopp  
Paper link: [SSRN working paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6336238)

The current paper uses the US application. Euro-area data and scripts are also included as auxiliary material.

## How To Run

Open `MRT_codes.Rproj` in RStudio, or start R from this folder. Then run:

```r
source("install_packages.R")
source("main.R")
```

The script uses relative paths and should be run from the root of this folder. It prints progress messages before each major data, model, figure, table, and workbook step.

By default, `main.R` sets:

```r
run_mode <- "outputs_from_saved"
```

This reproduces the paper figures and tables from the stored estimation results. Users who want to re-estimate the model can change `run_mode` at the top of `main.R`:

- `"outputs_from_saved"`: regenerate figures and tables from saved estimates.
- `"estimate_from_saved"`: re-estimate, initialized from saved estimates.
- `"estimate_from_generic"`: re-estimate, initialized from generic values.

To recompute the Gaussian-mixture smoothing inputs from the raw SPF histograms, set:

```r
indic.compute.mixture <- TRUE
```

The euro-area auxiliary workflow can be run separately with:

```r
source("main_EA.R")
```

## Software Requirements

The code was prepared for R and uses standard CRAN packages. Missing packages can be installed by running `install_packages.R`.

Main required packages include:

`Matrix`, `doParallel`, `foreach`, `doSNOW`, `lubridate`, `readxl`, `dplyr`, `tidyr`, `zoo`, `matrixcalc`, `tidyverse`, `numDeriv`, `MASS`, `optimx`, `sandwich`, `gdata`, `ggplot2`, `moments`, `mFilter`, `roll`, `Rcpp`, `RcppEigen`, `tis`, `fredr`, `openxlsx`, `extraDistr`, `writexl`, `eurostat`, and `readODS`.

The C++ routines used by the state-space calculations are compiled through `Rcpp::sourceCpp`.

## Data Sources

The replication package includes the US inputs used by the paper under `data/US/` and cached processed inputs under `data/processed/`.

Main US data sources include:

- Survey of Professional Forecasters density and point forecasts.
- FRED/BEA macroeconomic series used in the US application.
- PDS data used by the original project code.
- Cached processed workbooks and estimation inputs used by the default replication run.

The default run uses:

- `data/processed/Output.US.xlsx`
- `data/processed/survey.DATA.US.with.param.xlsx`
- `data/processed/all_mat_sparse_4_5_new.RData`
- `results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.only.errors.no.4Q.best.RData`
- `results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.4th.errors.no.4Q.best.final.new.RData`
- `results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.no.3rd.4th.no.4Q.RData`
- `results/US/US.Model.Implied.Distribution.RData`

The file `results/US/US.Model.gdp.trend.cycle.4.5.quarterly.best.corr.gdp.with.3rd.4th.RData` is used only as an initialization point when `run_mode <- "estimate_from_saved"`.

The auxiliary euro-area workflow uses:

- `data/EA/`
- `data/processed/Output.EA.xlsx`
- `data/processed/survey.DATA.EA.with.param.xlsx`
- `results/EA/EA.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.only.no.4Q.RData`
- `results/EA/EA.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.only.no.4Q.mean.higher.errors.RData`
- `results/EA/EA.Model.Implied.Distribution.RData`

Two US comparison figures can use FRED through the `fredr` package. To reproduce those live FRED calls, set a FRED API key before running `main.R`:

```r
Sys.setenv(FRED_API_KEY = "your_key_here")
```

If `FRED_API_KEY` is not set, the code skips the live FRED term-premium call and copies the paper-ready cached PDF instead. It also falls back to the local GDP series for the trend/cycle comparison.

### Data Redistribution

The MIT License applies to the code in this repository. Some raw data files are derived from third-party providers, including the Survey of Professional Forecasters, FRED/BEA, and euro-area statistical sources. These data remain subject to the terms and redistribution policies of the original providers.

The repository includes the input files and cached processed objects needed for academic replication of the paper. Users who redistribute, fork, or publish a modified copy of the repository should verify that any included third-party data files can be redistributed under the relevant provider terms. When in doubt, remove the raw data files and provide instructions for downloading them from the original sources.

## Main Outputs

The US workflow writes figures, workbooks, and tables to:

- `graphs/US_2024/Baseline/`
- `graphs/US_2024/All_moments/`
- `graphs/US_2024/No_3rd_4th/`
- `tables/`

The main-paper tables are written to `tables/`. `table_param.txt` is generated by the default US workflow; `table_distribution_divergence_VaR.txt` is shipped as a static paper output because the active code archive does not regenerate its final formatting.

The paper baseline is `graphs/US_2024/Baseline/`, which includes first-, second-, and third-order SPF moments and excludes fourth-order moments. `graphs/US_2024/All_moments/` keeps the specification that also includes fourth-order moments, and `graphs/US_2024/No_3rd_4th/` keeps the variance-only robustness specification.

Paper figures that are not generated by the active scripts are shipped directly in the relevant output folder, mainly `graphs/US_2024/Baseline/`. Additional figures produced by the original scripts are kept in the corresponding specification folders under `graphs/US_2024/`, since they may be useful for further analysis.

To check that all current paper outputs are covered by the replication package, run:

```r
source("RScripts/check_paper_outputs.R")
```

The EA workflow writes outputs to:

- `graphs/EA_2024/Baseline/`
- `tables/EA/Baseline/`

## Code Organization

- `main.R`: main US replication script.
- `main_EA.R`: auxiliary euro-area workflow.
- `Main.US.short.R`: backward-compatible wrapper that sources `main.R`.
- `install_packages.R`: package installation helper.
- `RScripts/0-procedure.R`, `RScripts/1-procedure.R`: shared R functions.
- `RScripts/function_loadings.cpp`, `RScripts/kalman_cpp.cpp`: C++ routines.
- `RScripts/US/`: US data, smoothing, estimation, and model-implied-distribution scripts.
- `RScripts/EA/`: euro-area data, smoothing, estimation, and model-implied-distribution scripts.
- `RScripts/make_outputs/`: figure, table, and workbook scripts.
- `RScripts/check_paper_outputs.R`: static audit of paper figures and tables.
- `data/US/`: US raw inputs.
- `data/EA/`: euro-area raw inputs.
- `data/processed/`: cached processed inputs and paper-ready cached artifacts.
- `results/US/`: stored US model and distribution objects used by `main.R`.
- `results/EA/`: stored EA model and distribution objects used by `main_EA.R`.

## Reproducibility Notes

The default run is designed to reproduce the tables and figures used in the paper without rerunning the slowest smoothing and estimation steps.

Re-estimation is available but can be time-consuming. It may also overwrite stored results if `save_estimation_results` is enabled in the relevant scripts.

Freshly generated PDFs may not be byte-identical to the Overleaf PDFs because PDF metadata can differ across systems and runs. Cached paper artifacts are stored exactly when the active code archive does not regenerate a paper figure directly.

## License

This code is distributed under the MIT License. See `LICENSE`. Third-party data files are not covered by the MIT License and remain subject to their original providers' terms.
