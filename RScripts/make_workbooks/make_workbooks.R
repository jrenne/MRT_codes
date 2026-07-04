library(openxlsx)

# fitted obs
data.obs.smoothed <- as.data.frame(KF.result4$fitted.obs.smoothed)
data.obs.smoothed <- cbind(observables.with.dates$date,data.obs.smoothed)
colnames(data.obs.smoothed) <- colnames(observables.with.dates)

# latent
data.latent <- as.data.frame(KF.result4$xi.tT)
data.latent <- cbind(observables.with.dates$date,data.latent)
colnames(data.latent) <- c("date","y1","y2","y3","y4", "y1_lag", "y2_lag", "y3_lag", "y4_lag",
                           "pi1", "pi2",  "pi3", "gdp1", "gdp2", "gdp3",
                           "z1", "z2", "z3", "z4", "z5")

# colnames(data.latent) <- c("date","y1","y2","y3","y4", "y1_lag", "y2_lag", "y3_lag", "y4_lag",
#                            "pi1", "pi2",  "pi3", "gdp1", "gdp2", "gdp3", "gdp4", "gdp5", "gdp6",
#                            "z1", "z2", "z3", "z4", "z5")

# Trend cycle
data.trend.cycle <- cbind(vec.dates,model.p.t,Trend.inflation,Cycle.inflation,model.gdp.t,Trend.gdp,Cycle.gdp)

#Supply demand
# Create data base for cycle decomposition (Inflation)
cycle.inflation.decomposition.data <- data.frame(vec.dates, Cycle.inflation=Cycle.inflation/100, Cycle.inflation.supply=Cycle.inflation.supply/100, Cycle.inflation.demand=Cycle.inflation.demand/100)

# Create data base for cycle decomposition (GDP)
cycle.gdp.decomposition.data <- data.frame(vec.dates,Cycle.gdp=Cycle.gdp/100, Cycle.gdp.supply=Cycle.gdp.supply/100, Cycle.gdp.demand=Cycle.gdp.demand/100) 

# Create data base for inflation decomposition 
inflation.y.o.y.decomposition.data <- data.frame(vec.dates, Inflation.y.o.y=model.inflation.annual, Inflation.y.o.y.supply=Inflation.y.o.y.supply, Inflation.y.o.y.demand=Inflation.y.o.y.demand)

# Create data base for gdp decomposition 
gdp.annual.decomposition.data <- data.frame(vec.dates, gdp.annual=model.gdp.annual, gdp.annual.supply=gdp.annual.supply, gdp.annual.demand=gdp.annual.demand)

# Database correlation
all.corr <- as.data.frame(all.corr)
all.corr$V1 <- as.Date(all.corr$V1)
colnames(all.corr) <- c("date", "correlation.1Q", "correlation.4Q", "correlation.8Q")

# Create a new workbook
wb <- createWorkbook()

# Add sheets to the workbook
addWorksheet(wb, "Observables")
addWorksheet(wb, "Fitted")
addWorksheet(wb, "Latent")
addWorksheet(wb, "TrendCycle")
addWorksheet(wb, "CycleDecompositionInflation")
addWorksheet(wb, "CycleDecompositionGDP")
addWorksheet(wb, "InflationYoYDecompostion")
addWorksheet(wb, "YoYGDPDecomposition")
addWorksheet(wb, "Correlation")

# Write data frames to sheets
writeData(wb, sheet = "Observables", observables.with.dates)
writeData(wb, sheet = "Fitted", data.obs.smoothed)
writeData(wb, sheet = "Latent", data.latent)
writeData(wb, sheet = "TrendCycle", data.trend.cycle)
writeData(wb, sheet = "CycleDecompositionInflation", cycle.inflation.decomposition.data)
writeData(wb, sheet = "CycleDecompositionGDP", cycle.gdp.decomposition.data)
writeData(wb, sheet = "InflationYoYDecompostion", inflation.y.o.y.decomposition.data)
writeData(wb, sheet = "YoYGDPDecomposition", gdp.annual.decomposition.data)
writeData(wb, sheet = "Correlation", all.corr)


# Save the workbook to a file
if(area=="US"){
  saveWorkbook(wb, paste0(path_graph, "Output.US.xlsx"), overwrite = TRUE)
}

if(area=="EA"){
  saveWorkbook(wb, paste0(path_graph, "Output.EA.xlsx"), overwrite = TRUE)
}

# Confirmation message
cat("Data has been written to Output.", area, ".xlsx with multiple sheets\n", sep = "")
