# ==============================================================================
# Static audit of paper outputs
# ==============================================================================
# This checks whether each US figure/table input by the current paper is covered
# by the replication package, either as an active generated output or as a
# static paper output shipped in the natural output folder.

paper_figures <- c(
  "Figures_US/figure_distributions_noModeled_gdp_2015.pdf",
  "Figures_US/figure_distributions_noModeled_gdp_2020.pdf",
  "Figures_US/figure_distributions_noModeled_infl_2015.pdf",
  "Figures_US/figure_distributions_noModeled_infl_2020.pdf",
  "Figures_US_May2025/No_4th/6.Fitted.SPF.US.8Q.beta.pe.pdf",
  "Figures_US_May2025/No_4th/11.Fitted.SPF.US.8Q.beta.var.pdf",
  "Figures_US_May2025/No_4th/16.Fitted.SPF.US.8Q.beta.k3rd.pdf",
  "Figures_US_May2025/No_4th/27.Fitted.SPF.US.G.8Q.beta.pe.pdf",
  "Figures_US_May2025/No_4th/32.Fitted.SPF.US.G.8Q.beta.var.pdf",
  "Figures_US_May2025/No_4th/37.Fitted.SPF.US.G.8Q.beta.k3rd.pdf",
  "Figures_US_May2025/No_4th/Probability.stagflation.y.o.y.4.6.8.Q.pdf",
  "Figures_US_May2025/No_4th/Y.t.pdf",
  "Figures_US_May2025/No_4th/z.t.pdf",
  "Figures_US_May2025/No_4th/figure_distributions_gdp_1985.pdf",
  "Figures_US_May2025/No_4th/figure_distributions_infl_1985.pdf",
  "Figures_US_May2025/No_4th/log.hcpi.gdp.decomposition.pdf",
  "Figures_US_May2025/No_4th/log.gdp.cyclical.decomposition.pdf",
  "Figures_US_May2025/No_4th/annual.gdp.growth.decomposition.pdf",
  "Figures_US_May2025/No_4th/log.hcpi.cyclical.decomposition.pdf",
  "Figures_US_May2025/No_4th/y.o.y.inflation.decomposition.pdf",
  "Figures_US_May2025/No_4th/correlation_singlePlot.pdf",
  "Figures_US_May2025/No_4th/correlation_compareTP.pdf",
  "Figures_US_May2025/No_3rd_4th/annual.gdp.growth.decomposition.pdf",
  "Figures_US_May2025/No_3rd_4th/y.o.y.inflation.decomposition.pdf",
  "Figures_US_May2025/All_moments/annual.gdp.growth.decomposition.pdf",
  "Figures_US_May2025/All_moments/y.o.y.inflation.decomposition.pdf"
)

static_baseline_figures <- c(
  "Figures_US/figure_distributions_noModeled_gdp_2015.pdf",
  "Figures_US/figure_distributions_noModeled_gdp_2020.pdf",
  "Figures_US/figure_distributions_noModeled_infl_2015.pdf",
  "Figures_US/figure_distributions_noModeled_infl_2020.pdf",
  "Figures_US_May2025/No_4th/Probability.stagflation.y.o.y.4.6.8.Q.pdf",
  "Figures_US_May2025/No_4th/correlation_compareTP.pdf"
)

map_generated_figure <- function(path) {
  file_name <- basename(path)
  if (path %in% static_baseline_figures) {
    return(file.path("graphs/US_2024/Baseline", file_name))
  }
  if (grepl("^Figures_US_May2025/No_4th/", path)) {
    return(file.path("graphs/US_2024/Baseline", file_name))
  }
  if (grepl("^Figures_US_May2025/No_3rd_4th/", path)) {
    return(file.path("graphs/US_2024/No_3rd_4th", file_name))
  }
  if (grepl("^Figures_US_May2025/All_moments/", path)) {
    return(file.path("graphs/US_2024/All_moments", file_name))
  }
  NA_character_
}

audit_one <- function(path) {
  target <- map_generated_figure(path)
  if (!is.na(target) && file.exists(target)) {
    status <- if (path %in% static_baseline_figures) "static_in_output_folder" else "generated_by_main"
    return(data.frame(path = path, status = status, package_path = target))
  }
  data.frame(path = path, status = "missing", package_path = NA_character_)
}

figure_audit <- do.call(rbind, lapply(paper_figures, audit_one))

generated_tables <- c("table_param.txt")
static_tables <- c("table_distribution_divergence_VaR.txt")

table_audit <- rbind(
  do.call(rbind, lapply(generated_tables, function(file_name) {
    target <- file.path("tables", file_name)
    data.frame(
      path = file.path("Tables", file_name),
      status = if (file.exists(target)) "generated_by_main" else "missing",
      package_path = target
    )
  })),
  do.call(rbind, lapply(static_tables, function(file_name) {
    target <- file.path("tables", file_name)
    data.frame(
      path = file.path("Tables", file_name),
      status = if (file.exists(target)) "static_in_output_folder" else "missing",
      package_path = target
    )
  }))
)

print(figure_audit, row.names = FALSE)
print(table_audit, row.names = FALSE)

if (any(figure_audit$status == "missing") || any(table_audit$status == "missing")) {
  stop("Some paper outputs are not covered by the package.")
}
