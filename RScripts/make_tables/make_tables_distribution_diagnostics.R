# ==============================================================================
# Distribution-fit diagnostic tables
# ==============================================================================

decimal <- 3

if (!exists("path_table")) {
  path_table <- "tables/US_2024/Baseline/"
}
dir.create(path_table, recursive = TRUE, showWarnings = FALSE)

specs <- list(
  baseline = list(
    label = "Baseline ($1^{st}$,$2^{nd}$, $3^{rd}$)",
    file = "results/US/US.Model.Implied.Distribution.with.3rd.only.RData"
  ),
  all_moments = list(
    label = "All Moments ($1^{st}$,$2^{nd}$,$3^{rd}$,$4^{th}$)",
    file = "results/US/US.Model.Implied.Distribution.all.moments.RData"
  ),
  no_higher = list(
    label = "No Higher-Order Mom. ($1^{st}$,$2^{nd}$)",
    file = "results/US/US.Model.Implied.Distribution.no.3rd.4th.RData"
  )
)

missing_files <- vapply(specs, function(spec) !file.exists(spec$file), logical(1))
if (any(missing_files)) {
  stop(
    "Missing implied-distribution file(s): ",
    paste(vapply(specs[missing_files], `[[`, character(1), "file"), collapse = ", "),
    call. = FALSE
  )
}

implied_distributions <- lapply(specs, function(spec) readRDS(spec$file))

normalize_pdf <- function(x) {
  x[x < 0] <- 0
  x / sum(x, na.rm = TRUE)
}

cdf_quantile <- function(x, cdf, prob) {
  cdf <- cummax(cdf)
  cdf <- pmin(pmax(cdf, 0), 1)
  approx(cdf, x, xout = prob, ties = "ordered", rule = 2)$y
}

summarise_vector <- function(x) {
  c(
    min = min(x, na.rm = TRUE),
    q25 = unname(quantile(x, 0.25, na.rm = TRUE)),
    median = median(x, na.rm = TRUE),
    q75 = unname(quantile(x, 0.75, na.rm = TRUE)),
    max = max(x, na.rm = TRUE)
  )
}

empty_metrics <- function() {
  list(tvd = numeric(), kld = numeric(), p05 = numeric(), p95 = numeric())
}

diagnostics <- lapply(names(specs), function(.) {
  setNames(
    replicate(2, setNames(replicate(4, empty_metrics(), simplify = FALSE), paste0(5:8, "Q")), simplify = FALSE),
    c("infl", "gdp")
  )
})
names(diagnostics) <- names(specs)

for (type_var in c("infl", "gdp")) {
  area_var <- if (type_var == "infl") 1 else 2
  for (i in seq_len(nrow(observables))) {
    current_date <- as.Date(vec.dates$date)[i]
    if (current_date %in% as.Date(c("1985-01-15", "1986-01-15", "1990-01-15"))) {
      next
    }
    horizon <- switch(format(current_date, "%m"),
                      "01" = 8,
                      "04" = 7,
                      "07" = 6,
                      "10" = 5,
                      NA)
    if (is.na(horizon)) {
      next
    }

    for (spec_name in names(specs)) {
      distri <- implied_distributions[[spec_name]]
      all_distribution <- compute.model.implied(
        current_date,
        type_var,
        horizon,
        distri$PDF.all,
        distri$PDF.x.all
      )

      pdf_implied <- normalize_pdf(all_distribution$pdf.implied)
      pdf_mixture <- normalize_pdf(all_distribution$pdf.mixture)
      cdf_implied <- cumsum(pdf_implied)
      cdf_mixture <- cumsum(pdf_mixture)

      horizon_name <- paste0(horizon, "Q")
      diagnostics[[spec_name]][[type_var]][[horizon_name]]$tvd <-
        c(diagnostics[[spec_name]][[type_var]][[horizon_name]]$tvd,
          total_variation_distance(pdf_implied, pdf_mixture))
      diagnostics[[spec_name]][[type_var]][[horizon_name]]$kld <-
        c(diagnostics[[spec_name]][[type_var]][[horizon_name]]$kld,
          KL_divergence(pdf_implied, pdf_mixture))

      implied_p05 <- cdf_quantile(all_distribution$x.implied, cdf_implied, 0.05)
      mixture_p05 <- cdf_quantile(all_distribution$x.mixture, cdf_mixture, 0.05)
      implied_p95 <- cdf_quantile(all_distribution$x.implied, cdf_implied, 0.95)
      mixture_p95 <- cdf_quantile(all_distribution$x.mixture, cdf_mixture, 0.95)
      diagnostics[[spec_name]][[type_var]][[horizon_name]]$p05 <-
        c(diagnostics[[spec_name]][[type_var]][[horizon_name]]$p05, abs(implied_p05 - mixture_p05))
      diagnostics[[spec_name]][[type_var]][[horizon_name]]$p95 <-
        c(diagnostics[[spec_name]][[type_var]][[horizon_name]]$p95, abs(implied_p95 - mixture_p95))
    }
  }
}

format_min_max <- function(x) {
  paste0("(", make.entry(x["min"], decimal), ", ", make.entry(x["max"], decimal), ")")
}

format_iqr <- function(x) {
  paste0("(", make.entry(x["q25"], decimal), ", ", make.entry(x["q75"], decimal), ")")
}

metric_values <- function(type_var, horizon, metric) {
  lapply(diagnostics, function(spec) spec[[type_var]][[paste0(horizon, "Q")]][[metric]])
}

metric_values_overall <- function(type_var, metric) {
  lapply(diagnostics, function(spec) {
    unlist(lapply(paste0(5:8, "Q"), function(h) spec[[type_var]][[h]][[metric]]), use.names = FALSE)
  })
}

add_metric_rows <- function(latex_table, values) {
  stats <- lapply(values, summarise_vector)
  latex_table <- rbind(
    latex_table,
    paste0("Min, Max&", paste(vapply(stats, format_min_max, character(1)), collapse = "&"), "\\\\"),
    paste0("25th, 75th&", paste(vapply(stats, format_iqr, character(1)), collapse = "&"), "\\\\"),
    paste0("Median&", paste(vapply(stats, function(x) make.entry(x["median"], decimal), character(1)), collapse = "&"), "\\\\"),
    "\\\\"
  )
  latex_table
}

make_divergence_table <- function() {
  latex_table <- rbind(
    "\\begin{table}[ph!]",
    "\\begingroup \\tiny",
    "\\caption{Distribution Divergence Comparison Table}",
    "\\label{tab:distribution.divergence}",
    "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}rrrrrr}",
    "\\hline",
    "\\hline",
    " \\multicolumn{2}{c}{a) Inflation} \\\\",
    paste0(" & \\multicolumn{2}{c}{", specs$baseline$label, "} & \\multicolumn{2}{c}{",
           specs$all_moments$label, "} & \\multicolumn{2}{c}{", specs$no_higher$label, "} \\\\"),
    "\\cmidrule(lr){2-3}",
    "\\cmidrule(lr){4-5}",
    "\\cmidrule(lr){6-7}",
    " & $d_{TV}$ & $d_{KL}$ & $d_{TV}$ & $d_{KL}$ & $d_{TV}$ & $d_{KL}$ \\\\",
    "\\hline"
  )

  for (type_var in c("infl", "gdp")) {
    if (type_var == "gdp") {
      latex_table <- rbind(latex_table, "\\hline", " \\multicolumn{2}{c}{b) Real GDP growth} \\\\", "\\\\")
    }
    for (horizon in 5:8) {
      latex_table <- rbind(
        latex_table,
        paste0(" \\multicolumn{1}{c}{horizon = ", horizon, "Q} \\\\"),
        "\\cmidrule(lr){1-1}"
      )
      values <- unlist(
        Map(function(tvd, kld) list(tvd, kld),
            metric_values(type_var, horizon, "tvd"),
            metric_values(type_var, horizon, "kld")),
        recursive = FALSE
      )
      latex_table <- add_metric_rows(latex_table, values)
    }
    latex_table <- rbind(
      latex_table,
      " \\multicolumn{1}{c}{Overall} \\\\",
      "\\cmidrule(lr){1-1}"
    )
    values <- unlist(
      Map(function(tvd, kld) list(tvd, kld),
          metric_values_overall(type_var, "tvd"),
          metric_values_overall(type_var, "kld")),
      recursive = FALSE
    )
    latex_table <- add_metric_rows(latex_table, values)
  }

  rbind(
    latex_table,
    "\\hline",
    "\\end{tabular*}",
    "\\begin{footnotesize}",
    "\\parbox{\\linewidth}{\\textit{Notes}: The table reports distribution-fit statistics comparing model-implied predictive distributions with the smooth distributions obtained from Gaussian mixtures fitted to SPF histograms.}",
    "\\end{footnotesize}",
    "\\endgroup",
    "\\end{table}"
  )
}

make_percentile_table <- function() {
  latex_table <- rbind(
    "\\begin{table}[ph!]",
    "\\begingroup \\tiny",
    "\\caption{Distance of model-implied percentiles from observed ones}",
    "\\label{tab:distribution.divergence.VaR}",
    "\\renewcommand{\\arraystretch}{0.6}",
    "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}rrrrrr}",
    "\\hline",
    "\\hline",
    paste0(" & \\multicolumn{2}{c}{", specs$baseline$label, "} & \\multicolumn{2}{c}{",
           specs$all_moments$label, "} & \\multicolumn{2}{c}{", specs$no_higher$label, "} \\\\"),
    "\\cmidrule(lr){2-3}",
    "\\cmidrule(lr){4-5}",
    "\\cmidrule(lr){6-7}",
    " & $5^{th}$ & $95^{th}$ & $5^{th}$ & $95^{th}$ & $5^{th}$ & $95^{th}$ \\\\",
    "\\hline",
    "\\\\",
    "\\multicolumn{2}{c}{(a) Inflation} \\\\",
    "\\\\"
  )

  for (type_var in c("infl", "gdp")) {
    if (type_var == "gdp") {
      latex_table <- rbind(latex_table, "\\\\", " \\multicolumn{2}{c}{(b) Real GDP growth} \\\\", "\\\\")
    }
    for (horizon in 5:8) {
      latex_table <- rbind(
        latex_table,
        paste0(" \\multicolumn{1}{c}{horizon = ", horizon, "Q} \\\\"),
        "\\cmidrule(lr){1-1}"
      )
      values <- unlist(
        Map(function(p05, p95) list(p05, p95),
            metric_values(type_var, horizon, "p05"),
            metric_values(type_var, horizon, "p95")),
        recursive = FALSE
      )
      latex_table <- add_metric_rows(latex_table, values)
    }
    latex_table <- rbind(
      latex_table,
      " \\multicolumn{1}{c}{Overall} \\\\",
      "\\cmidrule(lr){1-1}"
    )
    values <- unlist(
      Map(function(p05, p95) list(p05, p95),
          metric_values_overall(type_var, "p05"),
          metric_values_overall(type_var, "p95")),
      recursive = FALSE
    )
    latex_table <- add_metric_rows(latex_table, values)
  }

  rbind(
    latex_table,
    "\\hline",
    "\\end{tabular*}",
    "\\begin{footnotesize}",
    "\\parbox{\\linewidth}{\\textit{Notes}: The table reports the absolute difference between model-implied percentiles and those from the smooth distributions obtained from Gaussian mixtures fitted to SPF histograms.}",
    "\\end{footnotesize}",
    "\\endgroup",
    "\\end{table}"
  )
}

write(make_divergence_table(), file.path(path_table, "table_distribution_divergence.txt"))
write(make_percentile_table(), file.path(path_table, "table_distribution_divergence_VaR.txt"))
