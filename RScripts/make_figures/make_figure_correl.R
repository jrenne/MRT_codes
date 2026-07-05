# ==============================================================================
#       PLOTS CORRELATION
# ==============================================================================

# Correlations
X.tT <- KF.result4$xi.tT
HH <- c(1,4,8)
res.test     <- compute.var.cov.X.tT(Model.final, X.tT, HH)

first.part   <- matrix(c(Model.final$delta[,2],rep(0,q)),1,Model.final$n+Model.final$q) %x% 
  matrix(c(Model.final$delta[,1],rep(0,q)),1,Model.final$n+Model.final$q)

second.part  <- matrix(c(Model.final$delta[,1],rep(0,q)),1,Model.final$n+Model.final$q) %x% 
  matrix(c(Model.final$delta[,1],rep(0,q)),1,Model.final$n+Model.final$q)

third.part   <- matrix(c(Model.final$delta[,2],rep(0,q)),1,Model.final$n+Model.final$q) %x% 
  matrix(c(Model.final$delta[,2],rep(0,q)),1,Model.final$n+Model.final$q)

all.corr     <- matrix(NA,dim(KF.result4$xi.tT)[1],length(HH)+1)
all.corr[,1] <- vec.dates$date

make_percent_axis <- function(side, limits, by = 0.1, col = "black", col.axis = col) {
  ticks <- seq(limits[1], limits[2], by = by)
  axis(
    side,
    at = ticks,
    labels = paste0(round(100 * ticks), "%"),
    las = 1,
    col = col,
    col.axis = col.axis
  )
}

standardize_for_plot <- function(x) {
  center <- mean(x, na.rm = TRUE)
  scale <- sd(x, na.rm = TRUE)
  if (!is.finite(scale) || scale == 0) {
    scale <- 1
  }
  list(
    z = (x - center) / scale,
    center = center,
    scale = scale
  )
}

axis_positions <- function(ticks, standardization) {
  (ticks - standardization$center) / standardization$scale
}

draw_standardized_correlation_tp <- function(correlation_dates,
                                             correlation_values,
                                             term_premium_dates,
                                             term_premium_bps,
                                             term_premium_label) {
  corr_std <- standardize_for_plot(correlation_values)
  tp_std <- standardize_for_plot(term_premium_bps)
  y.lim <- range(c(corr_std$z, tp_std$z), finite = TRUE)
  y.pad <- 0.08 * diff(y.lim)
  y.lim <- y.lim + c(-y.pad, y.pad)

  correlation_ticks <- pretty(correlation_values, n = 6)
  correlation_tick_locations <- axis_positions(correlation_ticks, corr_std)
  term_premium_ticks <- pretty(term_premium_bps, n = 6)
  term_premium_tick_locations <- axis_positions(term_premium_ticks, tp_std)

  plot(as.Date(correlation_dates), corr_std$z, type = "n",
       ylim = y.lim,
       xlab = "", ylab = "",
       las = 1,
       yaxt = "n",
       lwd = length(HH) + 1 - u)

  axis(2,
       at = correlation_tick_locations,
       labels = paste0(round(100 * correlation_ticks), "%"),
       las = 1)
  axis(4,
       at = term_premium_tick_locations,
       labels = paste0(round(term_premium_ticks), " bps"),
       col = "dark grey",
       col.axis = "dark grey",
       las = 1)
  mtext("Correlation (%)", side = 2, line = 3.0, cex = 0.85)
  mtext("Term premium (bps)", side = 4, line = 3.8, cex = 0.85, col = "dark grey")

  grid()
  make_recessions()
  abline(h = 0, lty = 2)
  lines(as.Date(correlation_dates), corr_std$z,
        lwd = length(HH) + 1 - u)
  lines(as.Date(term_premium_dates), tp_std$z,
        lwd = 2, col = "dark grey")

  legend("topleft",
         c(paste(horiz_in_Q, "-quarter correlation (lhs)", sep = ""),
           term_premium_label),
         bg = "white",
         col = c("black", "dark grey"),
         lwd = 2,
         cex = 0.9)
}


path <- paste0(path_graph,"correlation.pdf")
pdf(path, width = 7, height = 7)
par(plt=c(.1,.95,.15,.8))
par(mfrow=c(length(HH),1))

for (u in 1:length(HH)) {
  
  res.cov.test <- first.part%*%res.test[,,u] 
  var.1 <- second.part%*%res.test[,,u] 
  var.2 <- third.part%*%res.test[,,u] 
  res.corr.test <- res.cov.test/(var.1^0.5*var.2^0.5)
  plot(vec.dates$date,c(res.corr.test), type="l",
       ylim=c(-.3,0.3),
       ylab="",
       main=paste0("Correlation Inflation GDP - h=", HH[u], " quarter(s)"),
       las=1,lwd=2)
  
  make_recessions()
  abline(h=0, lty=2)
  grid()
  
  all.corr[,u+1] <- c(res.corr.test)
  
}

dev.off()



# Figure with a single plot and different horizons: ----------------------------

path <- paste0(path_graph,"correlation_singlePlot.pdf")
pdf(path, width = 7, height = 4,pointsize = 11)
par(plt=c(.10,.95,.15,.95))
par(mfrow=c(1,1))
y.lim.single <- c(-0.2, 0.5)

for (u in 1:length(HH)) {
  
  res.cov.test <- first.part%*%res.test[,,u] 
  var.1 <- second.part%*%res.test[,,u] 
  var.2 <- third.part%*%res.test[,,u]
  res.corr.test <- res.cov.test/(var.1^0.5*var.2^0.5)
  if(u == 1){
    plot(as.Date(vec.dates$date),c(res.corr.test), type="l",
         #ylim=c(-.3,0.3),
         ylim = y.lim.single,
         xlab="",ylab="",
         las=1,
         yaxt = "n",
         lwd=length(HH)+1-u)
    make_percent_axis(2, y.lim.single)
    grid()
    
    make_recessions()
    abline(h=0, lty=2)
    grid()
    
    legend("topleft",
           paste(HH," quarters",sep=""),
           bg="white",
           lwd=length(HH):1,
           cex=0.9,
           title="Horizon:")
  }else{
    lines(as.Date(vec.dates$date),c(res.corr.test),
          lwd=length(HH)+1-u)
  }
}

dev.off()


# Figure comparing correlations and term premiums: -----------------------------

horiz_in_Q <- 4
u <- which(HH==horiz_in_Q)

if(area=="US"){
start.date <- "1968-10-01"
end.date   <- "2024-01-01"
freq       <- "q"
fred_key <- Sys.getenv("FRED_API_KEY")
if (nzchar(fred_key) && requireNamespace("fredr", quietly = TRUE)) {
  library(fredr)
  fredr_set_key(fred_key)
  start_date <- as.Date(start.date)
  end_date   <- as.Date(end.date)
  data_TP <- fredr(series_id = "THREEFYTP10",
                   observation_start = start_date,
                   observation_end   = end_date,
                   frequency = freq,
                   aggregation_method = "avg")
} else {
  fred_url <- "https://fred.stlouisfed.org/graph/fredgraph.csv?id=THREEFYTP10"
  data_TP <- read.csv(fred_url)
  colnames(data_TP) <- c("date", "value")
  data_TP$date <- as.Date(data_TP$date)
  data_TP$value <- as.numeric(data_TP$value)
  data_TP <- data_TP[data_TP$date >= as.Date(start.date) & data_TP$date <= as.Date(end.date), ]
  data_TP <- data_TP %>%
    mutate(date = as.Date(as.yearqtr(date))) %>%
    group_by(date) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
}

path <- paste0(path_graph,"correlation_compareTP.pdf")
pdf(path, width = 9, height = 4,pointsize = 11)
par(plt=c(.11,.86,.15,.95))
par(mfrow=c(1,1))

res.cov.test <- first.part%*%res.test[,,u] 
var.1 <- second.part%*%res.test[,,u] 
var.2 <- third.part%*%res.test[,,u] 
res.corr.test <- res.cov.test/(var.1^0.5*var.2^0.5)
draw_standardized_correlation_tp(
  vec.dates$date,
  c(res.corr.test),
  data_TP$date,
  100 * data_TP$value,
  "Kim-Wright 10-year term premium (rhs)"
)

dev.off()

}

if(area=="EA"){
  library(readODS)
  # download EUTERPE 10Y term premium (Euro Area)
  url <- "https://www.unive.it/pag/fileadmin/user_upload/progetti_ricerca/euterpe/documenti/estimates_2024_02.ods"
  tmp_file <- tempfile(fileext = ".ods")
  download.file(url, destfile = tmp_file, mode = "wb")
  euro_tp <- read_ods(tmp_file, sheet = 1, skip = 2)
  euro_tp <- euro_tp[, 1:11]
  colnames(euro_tp) <- c("date", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9", "Y10")
  
  data_TP <- euro_tp %>%
    mutate(date = as.yearmon(date, "%b-%y"),
           date= as.Date(date)) %>%
    dplyr::select(date,Y10) %>%
    rename(value=Y10)
  
  path <- paste0(path_graph,"correlation_compareTP.pdf")
  pdf(path, width = 9, height = 4,pointsize = 11)
  par(plt=c(.11,.86,.15,.95))
  par(mfrow=c(1,1))
  
  res.cov.test <- first.part%*%res.test[,,u] 
  var.1 <- second.part%*%res.test[,,u] 
  var.2 <- third.part%*%res.test[,,u] 
  res.corr.test <- res.cov.test/(var.1^0.5*var.2^0.5)
  draw_standardized_correlation_tp(
    vec.dates$date,
    c(res.corr.test),
    data_TP$date,
    10000 * data_TP$value,
    "EUTERPE 10-year term premium (rhs)"
  )
  
  dev.off()
  
}






# 
# 
# u <- 1
# res.cov.test <- first.part%*%res.test[,,u] 
# var.1 <- second.part%*%res.test[,,u] 
# var.2 <- third.part%*%res.test[,,u] 
# res.corr.test.1 <- res.cov.test/(var.1^0.5*var.2^0.5)
# max(res.corr.test.1)
# 
# u <- 2
# res.cov.test <- first.part%*%res.test[,,u] 
# var.1 <- second.part%*%res.test[,,u] 
# var.2 <- third.part%*%res.test[,,u] 
# res.corr.test.4 <- res.cov.test/(var.1^0.5*var.2^0.5)
# max(res.corr.test.4)
# 
# 
# path <- paste0("graphs/US_2024/correlation.pdf")
# pdf(path, width = 5, height = 4.5)
# b <- 1
# m.1 <- matrix(seq(1,b+1),nrow = (b+1),ncol = 1,byrow = TRUE)
# layout(mat = m.1, heights = c(rep(0.925/b,b),0.075))
# par(mar=c(2, 3, 1, 1))
# 
# 
# plot(vec.dates,c(res.corr.test.1), type="l", ylim=c(min(0,c(res.corr.test.1)),max(c(res.corr.test.1),0)), ylab="",
#      main="", lwd=2, col="darkgrey")
# abline(h=0, lty=2)
# 
# dev.off()
# 
# 
