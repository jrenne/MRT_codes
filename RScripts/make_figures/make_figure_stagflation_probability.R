# ==============================================================================
# Model-implied stagflation probabilities
# ==============================================================================

if (area != "US") {
  stop("The stagflation probability figure is currently defined for the US paper sample.")
}

set.seed(20260704)

if (!exists("show_stagflation_progress")) show_stagflation_progress <- TRUE
if (!exists("stagflation_nb_sim")) stagflation_nb_sim <- 2000
if (!exists("stagflation_batch_size")) stagflation_batch_size <- 500

nb_sim <- stagflation_nb_sim
batch_size <- min(stagflation_batch_size, nb_sim)
target_horizons <- c(4, 6, 8)
max_horizon <- max(target_horizons)

X0 <- KF.result4$xi.tT
T_sample <- nrow(X0)

if (show_stagflation_progress) {
  message(
    "MRT replication: simulating stagflation probabilities: ",
    format(nb_sim, big.mark = ","),
    " draws per date (",
    format(T_sample * nb_sim, big.mark = ","),
    " total paths), horizons ",
    paste(target_horizons, collapse = "/"),
    " quarters."
  )
}
simulation_start_time <- Sys.time()

delta_infl <- matrix(c(Model.final$delta[, 1], rep(0, Model.final$q)), ncol = 1)
delta_gdp <- matrix(c(Model.final$delta[, 2], rep(0, Model.final$q)), ncol = 1)

probability_data <- data.frame(date = as.Date(vec.dates$date))
event_counts <- matrix(
  0,
  nrow = T_sample,
  ncol = length(target_horizons),
  dimnames = list(NULL, paste0(target_horizons, "Q"))
)

batch_starts <- seq(1, nb_sim, by = batch_size)
for (batch_id in seq_along(batch_starts)) {
  batch_start <- batch_starts[batch_id]
  batch_end <- min(batch_start + batch_size - 1, nb_sim)
  current_batch_size <- batch_end - batch_start + 1
  sim_index <- rep(seq_len(T_sample), each = current_batch_size)
  Y0_sim <- X0[sim_index, seq_len(Model.final$n), drop = FALSE]
  z0_sim <- X0[sim_index, (Model.final$n + 1):(Model.final$n + Model.final$q), drop = FALSE]
  
  if (show_stagflation_progress) {
    message(
      "MRT replication: stagflation simulation batch ",
      batch_id, "/", length(batch_starts),
      " (draws ", format(batch_start, big.mark = ","),
      "-", format(batch_end, big.mark = ","),
      " per date)."
    )
  }
  simulated_X <- simul.model(Model.final, Y0_sim, z0_sim, max_horizon)
  
  for (horizon_index in seq_along(target_horizons)) {
    horizon <- target_horizons[horizon_index]
    X_h <- simulated_X[, , horizon + 1, drop = FALSE][, , 1]
    infl_yoy <- as.numeric(4 * Model.final$pi.bar[1] + X_h %*% delta_infl)
    gdp_yoy <- as.numeric(4 * Model.final$pi.bar[2] + X_h %*% delta_gdp)
    event <- infl_yoy > 4 & gdp_yoy < 0
    event_counts[, horizon_index] <- event_counts[, horizon_index] +
      as.numeric(tapply(event, sim_index, sum, na.rm = TRUE))
  }
}
probability_data[, paste0(target_horizons, "Q")] <- event_counts / nb_sim

if (show_stagflation_progress) {
  message(
    "MRT replication: stagflation simulations completed in ",
    round(as.numeric(difftime(Sys.time(), simulation_start_time, units = "secs")), 1),
    " seconds."
  )
}

path <- paste0(path_graph, "Probability.stagflation.y.o.y.4.6.8.Q.pdf")
if (show_stagflation_progress) {
  message("MRT replication: writing stagflation probability figure to ", path)
}
pdf(path, width = 7, height = 4, pointsize = 11)
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)
par(mar = c(3.0, 3.0, 0.6, 0.4), mgp = c(1.8, 0.5, 0), tcl = -0.25)

matplot(
  probability_data$date,
  probability_data[, paste0(target_horizons, "Q")],
  type = "n",
  lty = c(1, 2, 3),
  lwd = c(2.6, 3.2, 4.8),
  col = c("black", "dark grey", "grey45"),
  xlab = "",
  ylab = "",
  las = 1,
  ylim = c(0, max(probability_data[, paste0(target_horizons, "Q")], na.rm = TRUE) * 1.1)
)
grid()
make_recessions()
matlines(
  probability_data$date,
  probability_data[, paste0(target_horizons, "Q")],
  lty = c(1, 2, 3),
  lwd = c(2.6, 3.2, 4.8),
  col = c("black", "dark grey", "grey45")
)
legend(
  "topright",
  legend = paste(target_horizons, "quarters ahead"),
  col = c("black", "dark grey", "grey45"),
  lty = c(1, 2, 3),
  lwd = c(2.6, 3.2, 4.8),
  bg = "white",
  cex = 0.9
)

dev.off()
