# ==============================================================================
# Remove generated outputs from a local run
# ==============================================================================
# This keeps the repository lightweight after running the replication scripts.

generated_graphs <- list.files(
  "graphs",
  pattern = "[.](pdf|xlsx)$",
  recursive = TRUE,
  full.names = TRUE
)

generated_tables <- list.files(
  "tables",
  pattern = "[.]txt$",
  recursive = TRUE,
  full.names = TRUE
)

files_to_remove <- c(generated_graphs, generated_tables, "Rplots.pdf")
files_to_remove <- files_to_remove[file.exists(files_to_remove)]

if (length(files_to_remove) > 0) {
  unlink(files_to_remove)
}

message("Removed ", length(files_to_remove), " generated output files.")
