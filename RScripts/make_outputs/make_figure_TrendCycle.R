# ==============================================================================
#       PLOTS DECOMPOSITION TREND CYCLE
# ==============================================================================

res_TrendCycle  <- compute_TrendCycle(Model.final,
                                      KF.result4)
Trend.inflation <- res_TrendCycle$Trend.inflation
Trend.gdp       <- res_TrendCycle$Trend.gdp
Cycle.inflation <- res_TrendCycle$Cycle.inflation
Cycle.gdp       <- res_TrendCycle$Cycle.gdp
model.p.t       <- res_TrendCycle$model.p.t
model.gdp.t     <- res_TrendCycle$model.gdp.t

if(area=="US"){
  # Extract Output gap from FRED:
  library(fredr)
  start.date <- "1968-10-01"
  end.date   <- "2024-01-01"
  freq       <- "q"
  fred_key <- Sys.getenv("FRED_API_KEY")
  if (nzchar(fred_key)) {
    fredr_set_key(fred_key)
  } else {
    warning("FRED_API_KEY is not set; skipping the FRED output-gap comparison line.")
  }
  has_fred_output_gap <- nzchar(fred_key)
  start_date <- as.Date(start.date)
  end_date   <- as.Date(end.date)
  if (has_fred_output_gap) {
    data_GDP <- fredr(series_id = "GDPC1",
                      observation_start = start_date,
                      observation_end   = end_date,
                      frequency = freq,
                      aggregation_method = "avg")
    data_GDPPOT <- fredr(series_id = "GDPPOT",
                         observation_start = start_date,
                         observation_end   = end_date,
                         frequency = freq,
                         aggregation_method = "avg")
    z <- log(data_GDP$value/data_GDPPOT$value)
  } else {
    data_GDP <- read.csv("./data/US/raw/GDPC1.csv")
    date_column <- intersect(c("observation_date", "DATE", "date"), names(data_GDP))[1]
    if (is.na(date_column)) {
      stop("Could not identify the date column in data/US/raw/GDPC1.csv.")
    }
    data_GDP$date <- as.Date(data_GDP[[date_column]])
    data_GDP <- data_GDP[data_GDP$date >= start_date & data_GDP$date <= end_date, ]
    z <- rep(NA_real_, nrow(data_GDP))
  }
}

if(area=="EA"){
  # Import output_gap data
  output_gap <- read.csv("./data/EA/OECD.ECO.MAD,DSD_EO_114@DF_EO_114,1.0+EA17.GAP.csv", skip = 0) %>%
    rename(date = TIME_PERIOD,
           EA.output.gap = OBS_VALUE) %>%
    dplyr::select(date, EA.output.gap)
  
  # Keep only rows where the date includes a "Q" (quarterly data)
  output_gap_q <- output_gap[grepl("Q", output_gap$date), ]
  
  # Convert "1996-Q1" style to yearqtr class and sort
  output_gap_q$date <- as.Date(as.yearqtr(output_gap_q$date, format = "%Y-Q%q"))
  
  # Order the dataset by the converted date
  data_GDP <- output_gap_q[order(output_gap_q$date), ]
  z <- data_GDP$EA.output.gap/100

}

### INFLATION
# #cd534c red
# #cd534c yellow
# 0173c2 blue
# #868686 grey

path <- paste0(path_graph,"log.hcpi.gdp.decomposition.pdf")
pdf(path, width = 9, height = 5.5,pointsize = 12)

par(plt=c(.1, .95, .2, .85))
par(mfrow=c(2,2))

plot(as.Date(vec.dates$date), model.p.t,type='l', lwd=2,
     ylim=c(min(model.p.t, Trend.inflation/100,na.rm=TRUE),
            max(model.p.t, Trend.inflation/100,na.rm=TRUE)),
     las=1,
     main="(a.1) Price level and trend",xlab="",ylab="")
lines(as.Date(vec.dates$date), Trend.inflation/100, col="#868686", lwd=2)
grid()
make_recessions()

legend("bottomright",
       c("log price index","Trend"),
       lwd=2,lty=1,bg="white",
       col=c("black","#868686"),
       cex=0.86, horiz = T)

plot(as.Date(vec.dates$date), Cycle.inflation/100,type='l',col='#868686',
     lwd=2, xlab="", ylab="",main="(a.2) Cycle component of the price level",
     las=1,
     ylim=c(min(0, Cycle.inflation/100,na.rm=TRUE)-0.005,
            max(0, Cycle.inflation/100,na.rm=TRUE)+0.005))
grid()
make_recessions()
abline(h=0, lty=1,col="dark grey")


plot(as.Date(vec.dates$date), model.gdp.t,type='l', lwd=2,
     ylim=c(min(model.gdp.t, Trend.gdp/100,na.rm=TRUE),
            max(model.gdp.t, Trend.gdp/100,na.rm=TRUE)),
     las=1,
     main="(b.1) GDP and trend",xlab="",ylab="")
lines(as.Date(vec.dates$date), Trend.gdp/100, col="#868686", lwd=2)
grid()
make_recessions()

legend("bottomright",
       c("log GDP","Trend"),
       lwd=2,lty=1,bg="white",
       col=c("black","#868686"),
       cex=0.86, horiz = T)


if(area=="US"){
  plot(as.Date(vec.dates$date), Cycle.gdp/100,type='l',col='#868686',
       lwd=2, xlab="", ylab="",main="(b.2) Cycle component of GDP",
       las=1,
       ylim=c(min(0, Cycle.gdp/100,na.rm=TRUE)-0.005,
              max(0, Cycle.gdp/100,na.rm=TRUE)+0.005))
  grid()
  make_recessions()
  abline(h=0, lty=1,col="dark grey")
  if (has_fred_output_gap) {
    lines(data_GDP$date,z,lty=3,lwd=2)
    legend("topright",
           c("Model","FRED"),
           lwd=2,
           lty=c(1,3),
           bg="white",
           col=c('#868686',"black"),
           cex=0.86, horiz = T)
  } else {
    legend("topright",
           "Model",
           lwd=2,
           lty=1,
           bg="white",
           col='#868686',
           cex=0.86)
  }
}

if(area=="EA"){
  plot(as.Date(vec.dates$date), Cycle.gdp/100,type='l',col='#868686',
       lwd=2, xlab="", ylab="",main="(b.2) Cycle component of GDP",
       las=1,
       ylim=c(min(0, z,na.rm=TRUE)-0.005,
              max(0, z,na.rm=TRUE)+0.005))
  grid()
  make_recessions()
  abline(h=0, lty=1,col="dark grey")
  lines(data_GDP$date,z,lty=3,lwd=2)
  legend("bottomleft",
         c("Model","OECD"),
         lwd=2,
         lty=c(1,3),
         bg="white",
         col=c('#868686',"black"),
         cex=0.86, horiz = T)
}


dev.off()
