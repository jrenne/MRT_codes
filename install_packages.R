packages <- c(
  "Matrix",
  "doParallel",
  "foreach",
  "doSNOW",
  "lubridate",
  "readxl",
  "dplyr",
  "tidyr",
  "zoo",
  "matrixcalc",
  "tidyverse",
  "numDeriv",
  "MASS",
  "bitops",
  "RCurl",
  "optimx",
  "sandwich",
  "gdata",
  "stringr",
  "ggplot2",
  "moments",
  "mFilter",
  "roll",
  "Rcpp",
  "RcppEigen",
  "rlang",
  "tis",
  "fredr",
  "openxlsx",
  "extraDistr",
  "writexl",
  "eurostat",
  "readODS"
)

missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
