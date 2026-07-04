# make_recessions <- function(col="#AA55AA44"){
#   nber_dates <- nberDates()
#   start_recession_dates <- as.Date(as.character(nber_dates[,1]),"%Y%m%d")
#   end_recession_dates   <- as.Date(as.character(nber_dates[,2]),"%Y%m%d")
#   
#   for(i in 1:dim(nber_dates)[1]){
#     polygon(c(start_recession_dates[i],start_recession_dates[i],
#               end_recession_dates[i],end_recession_dates[i]),
#             c(-10000,+10000,+10000,-10000),border=NaN,col=col)
#   }
#   return(1)
# }

get_fred_recessions <- function(series_id) {
  fred_key <- Sys.getenv("FRED_API_KEY")
  if (!requireNamespace("fredr", quietly = TRUE) || !nzchar(fred_key)) {
    return(data.frame(start_date = as.Date(character()), end_date = as.Date(character())))
  }
  fredr::fredr_set_key(fred_key)
  rec_data <- fredr::fredr(series_id)
  rec_data %>%
    arrange(date) %>%
    mutate(
      start = ifelse(value == 1 & lag(value, 1) == 0, 1, 0),
      end = ifelse(value == 0 & lag(value, 1) == 1, 1, 0)
    ) %>%
    filter(start == 1 | end == 1) %>%
    mutate(start_date = date, end_date = date) %>%
    mutate(end_date = dplyr::lead(end_date, 1)) %>%
    filter(value == 1) %>%
    dplyr::select(start_date, end_date)
}

make_recessions <- function(col="#D9D9D966"){
  
  if (area == "EA"){
    # OECD based Recession Indicators for Euro Area from the Period following the 
    # Peak through the Trough 
    rec_data_plot <- get_fred_recessions("EUROREC")
  }
  
  if (area == "US"){
    # NBER based Recession Indicators for the United States from the Period 
    # following the Peak through the Trough
    nber_dates <- nberDates()
    rec_data_plot <- data.frame(
      start_date = as.Date(as.character(nber_dates[,1]), "%Y%m%d"),
      end_date   = as.Date(as.character(nber_dates[,2]), "%Y%m%d")
    )
    }
  
  if (area == "CH"){
    # OECD based Recession Indicators for Switzerland from the Period following the 
    # Peak through the Trough 
    rec_data_plot <- get_fred_recessions("CHERECD")
  }
  
  start_recession_dates <- rec_data_plot$start_date
  end_recession_dates   <- rec_data_plot$end_date
  
  if (nrow(rec_data_plot) == 0) {
    return(1)
  }
  for(i in 1:dim(rec_data_plot)[1]){
    polygon(c(start_recession_dates[i],start_recession_dates[i],
              end_recession_dates[i],end_recession_dates[i]),
            c(-10000,+10000,+10000,-10000),border=NaN,col=col)
  }
  return(1)
}


get_recessions_date <- function(area="EA"){
  
  if (area == "EA"){
    # OECD based Recession Indicators for Euro Area from the Period following the 
    # Peak through the Trough 
    rec_data_plot <- get_fred_recessions("EUROREC")
  }
  
  if (area == "US"){
    # NBER based Recession Indicators for the United States from the Period 
    # following the Peak through the Trough
    nber_dates <- nberDates()
    rec_data_plot <- data.frame(
      start_date = as.Date(as.character(nber_dates[,1]), "%Y%m%d"),
      end_date   = as.Date(as.character(nber_dates[,2]), "%Y%m%d")
    )
      }
  
  if (area == "CH"){
    # OECD based Recession Indicators for Switzerland from the Period following the 
    # Peak through the Trough 
    rec_data_plot <- get_fred_recessions("CHERECD")
  }
  
  start_date <- rec_data_plot$start_date
  end_date   <- rec_data_plot$end_date
  
  plot_data_rec <- data.frame(start_date=start_date,end_date=end_date)
 
  
  return(plot_data_rec)
}


make.entry <- function(x,decimal,dollar=1){
  format.nb <- paste("%.",decimal,"f",sep="")
  if(dollar==1){
    output <- paste("$",sprintf(format.nb,x),"$",sep="")
  }else{
    output <- sprintf(format.nb,x)
  }
  return(output)
}


compute_TrendCycle <- function(Model,KF.result){
  
  # Trend and Cycles charts
  indic.fact <- c(1:Model$m)
  Delta.Trend.inflation <- Model$pi.bar[1]+ KF.result$xi.tT[,indic.fact]%*%Model$delta.t[indic.fact,1]
  Cycle.inflation <- KF.result$xi.tT[,indic.fact]%*%Model$delta.c[indic.fact,1]
  Delta.Trend.gdp <- Model$pi.bar[2] + KF.result$xi.tT[,indic.fact]%*%Model$delta.t[indic.fact,2]
  Cycle.gdp <- KF.result$xi.tT[,indic.fact] %*% Model$delta.c[indic.fact,2]
  
  # Create the trend. To do so, we make the assumption that p.tm1 is perfectly modeled.
  # We use the real data for p.tm1 and we use c.t-1 computed. Other Methods below;
  # But suppose all p.t are perfectly modelized (which is false).
  # Trend.inflation <- log(Data.Trend.Cycle$EA.hcpi.deseasonalized) - Cycle.inflation
  # Trend.gdp <- log(Data.Trend.Cycle$EA.gdp.deseasonalized) - Cycle.gdp
  
  ## Use p.tm1 (data) and compute T.tm1 perfectly modeled
  first.date <- head(observables.with.dates$date,1) - months(3)
  
  if(area=="US"){
    p.tm1 <- DATA.US.all %>% filter(date == first.date) %>% inner_join(DATA.G.US, by="date")
    T.tm1.hcpi.initial <- p.tm1$US.log.cpi - KF.result$xi.tT[1,(Model$m+1):(Model$m+Model$m)]%*%Model$delta.c[1:Model$m,1]
    T.tm1.gdp.initial  <- p.tm1$US.log.gdp - KF.result$xi.tT[1,(Model$m+1):(Model$m+Model$m)]%*%Model$delta.c[1:Model$m,2]
  } 
  
  if(area=="EA"){
    p.tm1 <- DATA %>% filter(date == "1998-10-15") %>% inner_join(DATA.G, by="date")
    T.tm1.hcpi.initial <- log(p.tm1$EA.hcpi.deseasonalized)*100 - KF.result4$xi.tT[1,(Model.final$m+1):(Model.final$m+Model.final$m)]%*%Model.final$delta.c[1:Model.final$m,1]
    T.tm1.gdp.initial <- log(p.tm1$EA.gdp.deseasonalized)*100 - KF.result4$xi.tT[1,(Model.final$m+1):(Model.final$m+Model.final$m)]%*%Model.final$delta.c[1:Model.final$m,2]
  }
  
  # Compute T.t for the first period
  Trend.inflation <- T.tm1.hcpi.initial + Delta.Trend.inflation[1]
  Trend.gdp       <- T.tm1.gdp.initial  + Delta.Trend.gdp[1]
  
  # Loop to compute the trend at all periods
  for(i in 2:length(Delta.Trend.inflation)){
    Trend.inflation[i] <- Trend.inflation[i-1] + Delta.Trend.inflation[i]
    Trend.gdp[i] <- Trend.gdp[i-1] + Delta.Trend.gdp[i]
  }
  
  # Compute modeled p.t and gdp.t
  model.p.t <- (Trend.inflation + Cycle.inflation)/100
  model.gdp.t <- (Trend.gdp + Cycle.gdp)/100
  
  return(list(Trend.inflation = Trend.inflation,
              Trend.gdp       = Trend.gdp,
              Cycle.inflation = Cycle.inflation,
              Cycle.gdp       = Cycle.gdp,
              model.p.t       = model.p.t,
              model.gdp.t     = model.gdp.t))
}


compute_growth_rates <- function(Model,KF.result){
  
  # Compute quarterly inflation and gdp
  model.base.inflation.quarter <- Model$pi.bar[1] + KF.result$xi.tT[,1:Model$m]%*%Model$delta.t[1:Model$m,1] +
    (KF.result$xi.tT[,1:Model$m]- KF.result$xi.tT[,(Model$m+1):(Model$m+Model$m)])%*%Model$delta.c[1:Model$m,1]
  model.base.gdp.quarter <- Model$pi.bar[2] + KF.result$xi.tT[,1:Model$m]%*%Model$delta.t[1:Model$m,2] +
    (KF.result$xi.tT[,1:Model$m]- KF.result$xi.tT[,(Model$m+1):(Model$m+Model$m)])%*%Model$delta.c[1:Model$m,2]
  # 
  # Compute y.o.y inflation
  model.inflation.annual.1 <- rollsum(model.base.inflation.quarter,4)
  model.inflation.annual <- Model$pi.bar[1]*4 + KF.result$xi.tT%*% matrix(c(Model$delta[,1],rep(0,Model$q)),Model$n+Model$q,1)
  
  # Compute y.o.y gdp and annual gdp (y.o.y)
  model.gdp.y.o.y <- rollsum(model.base.gdp.quarter,4)
  model.gdp.annual <- Model$pi.bar[2]*4 + KF.result$xi.tT%*% matrix(c(Model$delta[,2],rep(0,Model$q)),Model$n+Model$q,1)
  
  return(list(model.inflation.annual.1     = model.inflation.annual.1,
              model.inflation.annual       = model.inflation.annual,
              model.gdp.y.o.y              = model.gdp.y.o.y,
              model.gdp.annual             = model.gdp.annual,
              model.base.inflation.quarter = model.base.inflation.quarter,
              model.base.gdp.quarter       = model.base.gdp.quarter))
}


compute_decompositions <- function(Model,KF.result){
  
  indic.fact <- c(1:Model$m)
  
  ## Use p.tm1 (data) and compute T.tm1 perfectly modeled
  first.date <- head(observables.with.dates$date,1) - months(3)
  
  if(area=="US"){
    p.tm1 <- DATA.US.all %>% filter(date == first.date) %>% inner_join(DATA.G.US, by="date")
    T.tm1.hcpi.initial <- p.tm1$US.log.cpi - KF.result$xi.tT[1,(Model$m+1):(Model$m+Model$m)]%*%Model$delta.c[1:Model$m,1]
    T.tm1.gdp.initial  <- p.tm1$US.log.gdp - KF.result$xi.tT[1,(Model$m+1):(Model$m+Model$m)]%*%Model$delta.c[1:Model$m,2]
  } 
  
  if(area=="EA"){
    p.tm1 <- DATA %>% filter(date == "1998-10-15") %>% inner_join(DATA.G, by="date")
    T.tm1.hcpi.initial <- log(p.tm1$EA.hcpi.deseasonalized)*100 - KF.result4$xi.tT[1,(Model.final$m+1):(Model.final$m+Model.final$m)]%*%Model.final$delta.c[1:Model.final$m,1]
    T.tm1.gdp.initial <- log(p.tm1$EA.gdp.deseasonalized)*100 - KF.result4$xi.tT[1,(Model.final$m+1):(Model.final$m+Model.final$m)]%*%Model.final$delta.c[1:Model.final$m,2]
  }
  
  bothPositive <- function(row) {
    all(row > 0)
  }
  
  # Apply the function to each row of the matrix
  result1 <- apply(Model$delta.c, 1, bothPositive)
  result2 <- apply(Model$delta.t, 1, bothPositive)
  demand.factor.indic <- which(result1)
  supply.factor.indic <- which(!result1)
  
  # Index for supply and demand
  supply.indic.cycle <- supply.factor.indic
  demand.indic.cycle <- demand.factor.indic
  supply.indic.trend <- supply.factor.indic
  demand.indic.trend <- demand.factor.indic
  
  # Construct the delta.c and delta.t capturing demand and supply effects
  delta.c.supply <- Model$delta.c
  delta.c.supply[demand.indic.cycle,] <- 0
  delta.c.demand <- Model$delta.c
  delta.c.demand[supply.indic.cycle,] <- 0
  delta.t.supply <- Model$delta.t
  delta.t.supply[demand.indic.trend,] <- 0
  delta.t.demand <- Model$delta.t
  delta.t.demand[supply.indic.trend,] <- 0
  
  # Computes cycle parts associated to demand and supply factors
  Cycle.inflation.supply <- KF.result$xi.tT[,indic.fact]%*%delta.c.supply[indic.fact,1]
  Cycle.inflation.demand <- KF.result$xi.tT[,indic.fact]%*%delta.c.demand[indic.fact,1]
  Cycle.gdp.supply <- KF.result$xi.tT[,indic.fact]%*%delta.c.supply[indic.fact,2]
  Cycle.gdp.demand <- KF.result$xi.tT[,indic.fact]%*%delta.c.demand[indic.fact,2]
  
  # Computes delta trend parts associated to demand and supply factors
  Delta.Trend.inflation.supply <- KF.result$xi.tT[,indic.fact]%*%delta.t.supply[indic.fact,1]
  Delta.Trend.inflation.demand <- KF.result$xi.tT[,indic.fact]%*%delta.t.demand[indic.fact,1]
  Delta.Trend.gdp.supply <- KF.result$xi.tT[,indic.fact]%*%delta.t.supply[indic.fact,2]
  Delta.Trend.gdp.demand <- KF.result$xi.tT[,indic.fact]%*%delta.t.demand[indic.fact,2]
  
  ## Use p.tm1 (data) and compute T.tm1 prefectly modeled
  T.tm1.hcpi.initial.supply <-  - KF.result$xi.tT[1,(Model$m+1):(2*Model$m)]%*%delta.c.supply[indic.fact,1]
  T.tm1.hcpi.initial.demand <-  - KF.result$xi.tT[1,(Model$m+1):(2*Model$m)]%*%delta.c.demand[indic.fact,1]
  if(area=="US"){
    T.tm1.hcpi.initial.bis <- T.tm1.hcpi.initial - log(p.tm1$US.cpi)*100
  }
  if(area=="EA"){
    T.tm1.hcpi.initial.bis <- T.tm1.hcpi.initial - log(p.tm1$EA.hcpi.deseasonalized)*100
  }
  T.tm1.gdp.initial.supply <-  - KF.result$xi.tT[1,(Model$m+1):(2*Model$m)]%*%delta.c.supply[indic.fact,2]
  T.tm1.gdp.initial.demand <-  - KF.result$xi.tT[1,(Model$m+1):(2*Model$m)]%*%delta.c.demand[indic.fact,2]
  if(area=="US"){
    T.tm1.gdp.initial.bis <- T.tm1.gdp.initial - log(p.tm1$US.gdp)*100
  }
  if(area=="EA"){
    T.tm1.gdp.initial.bis <- T.tm1.gdp.initial - log(p.tm1$EA.gdp.deseasonalized)*100
  }
  
  # Compute T.t for the first period
  Trend.inflation.supply <- T.tm1.hcpi.initial.supply + Delta.Trend.inflation.supply[1]
  Trend.inflation.demand <- T.tm1.hcpi.initial.demand + Delta.Trend.inflation.demand[1]
  Trend.gdp.supply <- T.tm1.gdp.initial.supply + Delta.Trend.gdp.supply[1]
  Trend.gdp.demand <- T.tm1.gdp.initial.demand + Delta.Trend.gdp.demand[1]
  
  # Loop to compute the trend at all periods
  for(i in 2:length(Delta.Trend.inflation.supply)){
    Trend.inflation.supply[i] <- Trend.inflation.supply[i-1] + Delta.Trend.inflation.supply[i]
    Trend.inflation.demand[i] <- Trend.inflation.demand[i-1] + Delta.Trend.inflation.demand[i]
    Trend.gdp.supply[i] <- Trend.gdp.supply[i-1] + Delta.Trend.gdp.supply[i]
    Trend.gdp.demand[i] <- Trend.gdp.demand[i-1] + Delta.Trend.gdp.demand[i]
  }
  
  # Check
  if(area=="US"){
    Trend.inflation.bis <-Trend.inflation.supply + Trend.inflation.demand + log(p.tm1$US.cpi)*100 + seq(1,length(Trend.inflation),by=1)*Model$pi.bar[1]
    round(Trend.inflation.bis,5) == round(Trend.inflation,5)
    Trend.gdp.bis <- Trend.gdp.supply + Trend.gdp.demand + log(p.tm1$US.gdp)*100 + seq(1,length(Trend.gdp),by=1)*Model$pi.bar[2]
    round(Trend.gdp.bis,5) == round(Trend.gdp,5)
  }
  if(area=="EA"){
    Trend.inflation.bis <-Trend.inflation.supply + Trend.inflation.demand + log(p.tm1$EA.hcpi.deseasonalized)*100 + seq(1,length(Trend.inflation),by=1)*Model.final$pi.bar[1]
    round(Trend.inflation.bis,5) == round(Trend.inflation,5)
    Trend.gdp.bis <- Trend.gdp.supply + Trend.gdp.demand + log(p.tm1$EA.gdp.deseasonalized)*100 + seq(1,length(Trend.gdp),by=1)*Model.final$pi.bar[2]
    round(Trend.gdp.bis,5) == round(Trend.gdp,5)
  }
  
  # Create data base for cycle decomposition (Inflation)
  Cycle.inflation.decomposition <- data.frame(vec.dates, Cycle.inflation=Cycle.inflation/100, Cycle.inflation.supply=Cycle.inflation.supply/100, Cycle.inflation.demand=Cycle.inflation.demand/100) %>%
    gather(group, value, Cycle.inflation.supply:Cycle.inflation.demand)
  
  # Create data base for cycle decomposition (GDP)
  Cycle.gdp.decomposition <- data.frame(vec.dates,Cycle.gdp=Cycle.gdp/100, Cycle.gdp.supply=Cycle.gdp.supply/100, Cycle.gdp.demand=Cycle.gdp.demand/100) %>%
    gather(group, value, Cycle.gdp.supply:Cycle.gdp.demand)
  
  #diff.test<-log(p.tm1$EA.hcpi.deseasonalized)*100 + seq(1,length(Trend.inflation),by=1)*Model$pi.bar[1]
  # Create data base for trend decomposition (Inflation)
  if(area=="US"){
    Trend.inflation.decomposition <- data.frame(vec.dates,Trend.inflation=Trend.inflation/100, Trend.inflation.supply=Trend.inflation.supply/100, Trend.inflation.demand=Trend.inflation.demand/100) %>%
      mutate(Trend.inflation.supply.demand=Trend.inflation.supply + Trend.inflation.demand,
             rho.effect = (Trend.inflation-Trend.inflation.supply.demand -log(p.tm1$US.cpi)*100)/100 )  %>%
      gather(group, value, Trend.inflation.supply:Trend.inflation.demand)
    
    # Create data base for trend decomposition (GDP)
    Trend.gdp.decomposition <- data.frame(vec.dates,Trend.gdp=Trend.gdp/100, Trend.gdp.supply=Trend.gdp.supply/100, Trend.gdp.demand=Trend.gdp.demand/100) %>%
      mutate(Trend.gdp.supply.demand=Trend.gdp.supply + Trend.gdp.demand,
             rho.effect = (Trend.gdp-Trend.gdp.supply.demand -log(p.tm1$US.gdp)*100 )/100)  %>%
      gather(group, value, Trend.gdp.supply:Trend.gdp.demand)
  }
  
  if(area=="EA"){
    # Create data base for trend decomposition (Inflation)
    Trend.inflation.decomposition <- data.frame(vec.dates,Trend.inflation=Trend.inflation/100, Trend.inflation.supply=Trend.inflation.supply/100, Trend.inflation.demand=Trend.inflation.demand/100) %>%
      mutate(Trend.inflation.supply.demand=Trend.inflation.supply + Trend.inflation.demand,
             rho.effect = (Trend.inflation-Trend.inflation.supply.demand -log(p.tm1$EA.hcpi.deseasonalized)*100)/100 )  %>%
      gather(group, value, Trend.inflation.supply:Trend.inflation.demand) 
    
    # Create data base for trend decomposition (GDP)
    Trend.gdp.decomposition <- data.frame(vec.dates,Trend.gdp=Trend.gdp/100, Trend.gdp.supply=Trend.gdp.supply/100, Trend.gdp.demand=Trend.gdp.demand/100) %>%
      mutate(Trend.gdp.supply.demand=Trend.gdp.supply + Trend.gdp.demand,
             rho.effect = (Trend.gdp-Trend.gdp.supply.demand -log(p.tm1$EA.gdp.deseasonalized)*100 )/100)  %>%
      gather(group, value, Trend.gdp.supply:Trend.gdp.demand) 
  }

  
  # Y.O.Y Inflation
  Inflation.q.o.q.supply <- Delta.Trend.inflation.supply + Cycle.inflation.supply - KF.result4$xi.tT[,(Model.final$m+1):(2*Model.final$m)]%*%delta.c.supply[indic.fact,1]
  Inflation.q.o.q.demand <- Delta.Trend.inflation.demand + Cycle.inflation.demand - KF.result4$xi.tT[,(Model.final$m+1):(2*Model.final$m)]%*%delta.c.demand[indic.fact,1]
  cbind(Model.final$pi.bar[1] + Inflation.q.o.q.demand+Inflation.q.o.q.supply, model.base.inflation.quarter)
  
  Inflation.q.o.q.decomposition <- data.frame(vec.dates, Inflation.q.o.q=model.base.inflation.quarter, Inflation.q.o.q.supply=Inflation.q.o.q.supply, Inflation.q.o.q.demand=Inflation.q.o.q.demand) %>%
    gather(group, value, Inflation.q.o.q.supply:Inflation.q.o.q.demand)
  
  Inflation.y.o.y.supply <- rollsum(Inflation.q.o.q.supply,4, align = "right", fill = NA)
  Inflation.y.o.y.demand <- rollsum(Inflation.q.o.q.demand,4, align = "right", fill = NA)
  cbind(4*Model.final$pi.bar[1] + Inflation.y.o.y.supply + Inflation.y.o.y.demand,model.inflation.annual)
  
  #delta.inf.demand <- matrix(c(Model.final$delta.q[,1],rep(0,q)),Model.final$n+Model.final$q,1)
  #KF.result4$xi.tT%*%delta.inf.demand
  #delta.inf.demand[c(supply.indic.cycle,Model.final$m+supply.indic.cycle),] <- 0    
  
  Inflation.y.o.y.decomposition <- data.frame(vec.dates, Inflation.y.o.y=model.inflation.annual, Inflation.y.o.y.supply=Inflation.y.o.y.supply, Inflation.y.o.y.demand=Inflation.y.o.y.demand) %>%
    gather(group, value, Inflation.y.o.y.supply:Inflation.y.o.y.demand)
  
  # Annual GDP growth
  gdp.q.o.q.supply <- Delta.Trend.gdp.supply + Cycle.gdp.supply - KF.result4$xi.tT[,(Model.final$m+1):(2*Model.final$m)]%*%delta.c.supply[indic.fact,2]
  gdp.q.o.q.demand <- Delta.Trend.gdp.demand + Cycle.gdp.demand - KF.result4$xi.tT[,(Model.final$m+1):(2*Model.final$m)]%*%delta.c.demand[indic.fact,2]
  cbind(Model.final$pi.bar[2] + gdp.q.o.q.demand+gdp.q.o.q.supply, model.base.gdp.quarter)
  
  gdp.q.o.q.decomposition <- data.frame(vec.dates, gdp.q.o.q=model.base.gdp.quarter, gdp.q.o.q.supply=gdp.q.o.q.supply, gdp.q.o.q.demand=gdp.q.o.q.demand) %>%
    gather(group, value, gdp.q.o.q.supply:gdp.q.o.q.demand)
  
  gdp.annual.supply <- rollsum(gdp.q.o.q.supply,4, align = "right", fill = NA)
  gdp.annual.demand <- rollsum(gdp.q.o.q.demand,4, align = "right", fill = NA)
  #gdp.annual.supply <- roll_sum(gdp.q.o.q.supply, width = 7, weights = c(1,2,3,4,3,2,1))/4
  #gdp.annual.demand <- roll_sum(gdp.q.o.q.demand, width = 7, weights = c(1,2,3,4,3,2,1))/4
  cbind(4*Model.final$pi.bar[2] + gdp.annual.supply + gdp.annual.demand, model.gdp.annual)
  
  gdp.annual.decomposition <- data.frame(vec.dates, gdp.annual=model.gdp.annual, gdp.annual.supply=gdp.annual.supply, gdp.annual.demand=gdp.annual.demand) %>%
    gather(group, value, gdp.annual.supply:gdp.annual.demand)
  
  return(list(
    Cycle.inflation.supply        = Cycle.inflation.supply,
    Cycle.inflation.demand        = Cycle.inflation.demand,
    Cycle.gdp.supply              = Cycle.gdp.supply,
    Cycle.gdp.demand              = Cycle.gdp.demand,
    Inflation.y.o.y.supply        = Inflation.y.o.y.supply,
    Inflation.y.o.y.demand        = Inflation.y.o.y.demand,
    gdp.annual.demand             = gdp.annual.demand,
    gdp.annual.supply             = gdp.annual.supply,
    Cycle.inflation.decomposition = Cycle.inflation.decomposition,
    Cycle.gdp.decomposition       = Cycle.gdp.decomposition,
    Trend.inflation.decomposition = Trend.inflation.decomposition,
    Trend.gdp.decomposition       = Trend.gdp.decomposition,
    Inflation.q.o.q.decomposition = Inflation.q.o.q.decomposition,
    gdp.q.o.q.decomposition      = gdp.q.o.q.decomposition,
    Inflation.y.o.y.decomposition = Inflation.y.o.y.decomposition,
    gdp.annual.decomposition      = gdp.annual.decomposition
  ))
}


make.model.implied.plot <- function(date,
                                    type.var.implied, horizon,
                                    display.implied=TRUE,
                                    display.mixture=TRUE,
                                    display.legend=TRUE){
  
  if(type.var.implied=="infl"){
    
    if(area=="US"){
    
    if(date <= "1984-10-15"){
      period <- "1981_1985"
    } else if(date <= "1991-10-15"){
      period <- "1985_1992"
    } else if(date <= "2013-10-15"){
      period <- "1992_2014"
    } else if(date > "2013-10-15"){
      period <- "post_2014"
    }
    
    area.var <- 1
    
    eval(parse(text = gsub(" ","", paste('survey.bins.distri <- US.SPF.DISTRI.',
                                         paste(horizon,sep = ""), 'Q.', paste(period,sep = ""), ' %>% left_join(survey.DATA.US.with.param %>% dplyr::select(date, contains("',
                                         paste(horizon,sep = ""), 'Q" )), by="date")',
                                         sep="")))) 
    
    eval(parse(text = gsub(" ","", paste('nbr.bins <- dim(US.SPF.DISTRI.',
                                         paste(horizon,sep = ""), 'Q.', paste(period,sep = ""),')[2]', sep="")))) 
    }
    
    if(area=="EA"){
      
      period <- paste0(year(date), "Q", quarter(date))
      
      if(horizon %in% c(1:4)){
        horizon.calendar <- "cy"
      }
      if(horizon %in% c(5:8)){
        horizon.calendar <- "ny"
      }
      
      area.var <- 1
      
      eval(parse(text = gsub(" ","", paste('survey.bins.distri <- EA.SPF.DISTRI.',
                                           paste(horizon.calendar,sep = ""), '.',paste(period,sep = ""), '.avg %>% left_join(survey.DATA.new.with.param %>% dplyr::select(date, contains(".',
                                           paste(horizon,sep = ""), 'Q" )), by="date")',
                                           sep="")))) 
      
      eval(parse(text = gsub(" ","", paste('nbr.bins <- dim(EA.SPF.DISTRI.',
                                           paste(horizon.calendar,sep = ""), '.',paste(period,sep = ""),'.avg)[2]', sep="")))) 
      
    }
    
  } else if(type.var.implied=="gdp"){
    
    if(area=="US"){
      
    if(date <= "1991-10-15"){
      period <- "1981_1992"
    } else if(date <= "2009-01-15"){
      period <- "1992_2009"
    } else if(date <= "2020-01-15"){
      period <- "2009_2020"
    } else if(date > "2020-01-15"){
      period <- "post_2020"
    }
    
    area.var <- 2
    
    eval(parse(text = gsub(" ","", paste('survey.bins.distri <- US.G.SPF.DISTRI.',
                                         paste(horizon,sep = ""), 'Q.', paste(period,sep = ""), ' %>% left_join(survey.DATA.US.G.with.param %>% dplyr::select(date, contains(".',
                                         paste(horizon,sep = ""), 'Q" )), by="date")',
                                         sep="")))) 
    
    eval(parse(text = gsub(" ","", paste('nbr.bins <- dim(US.G.SPF.DISTRI.',
                                         paste(horizon,sep = ""), 'Q.', paste(period,sep = ""),')[2]', sep="")))) 
    }
    
    if(area=="EA"){
      
      period <- paste0(year(date), "Q", quarter(date))
      
      if(horizon %in% c(1:4)){
        horizon.calendar <- "cy"
      }
      if(horizon %in% c(5:8)){
        horizon.calendar <- "ny"
      }
      
      area.var <- 2
      
      eval(parse(text = gsub(" ","", paste('survey.bins.distri <- EA.G.SPF.DISTRI.',
                                           paste(horizon.calendar,sep = ""), '.',paste(period,sep = ""), '.avg %>% left_join(survey.DATA.G.new.with.param %>% dplyr::select(date, contains(".',
                                           paste(horizon,sep = ""), 'Q" )), by="date")',
                                           sep="")))) 
      
      eval(parse(text = gsub(" ","", paste('nbr.bins <- dim(EA.G.SPF.DISTRI.',
                                           paste(horizon.calendar,sep = ""), '.',paste(period,sep = ""),'.avg)[2]', sep="")))) 
      
    }
    
  } else{
    stop("Choose a correct type of variable, i.e. infl or gdp")
  }
  
  # ============================================================================
  # ============================================================================
  freq <- 4
  horizon.y <- horizon-(freq-1)
  # ============================================================================
  # ============================================================================
  
  all.data.new <- survey.bins.distri[survey.bins.distri$date==as.Date(date),]
  Date <- as.Date(all.data.new[,1])
  
  if(area=="US"){
    start.bin <- 2
  }
  
  if(area=="EA"){
    start.bin <- 4 
  }
  
  if(is.na(sum(all.data.new[,start.bin:nbr.bins]*100))){
    stop("There is no observations available. Choose a valid date")
  }
  
  data1 <- all.data.new[,start.bin:nbr.bins]*100
  if(area=="US"){
    data1 <- data1[length(data1):1] # reverse order
  }
  time <- which(Date==as.matrix(vec.dates))
  
  param1 <- as.numeric(all.data.new[,grep("param", colnames(as.data.frame(all.data.new)))])
  max.bin <- as.numeric(all.data.new[,grep("max.bin", colnames(as.data.frame(all.data.new)))])
  min.bin <- as.numeric(all.data.new[,grep("min.bin", colnames(as.data.frame(all.data.new)))])
  x1 <- seq(min.bin-0.5,max.bin+0.5,by=0.05)
  
  q <- round(as.numeric(format(Date, "%m"))/3) + 1
  t <- format(Date, "%Y")
  
  lower.bound.interval <- round(sort(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",colnames(data1))))))*2)/2
  lower.bound.interval[1] <- min.bin
  upper.bound.interval <- round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(data1)))))*2)/2 
  upper.bound.interval[length(upper.bound.interval)] <- max.bin
  all.intervals <- upper.bound.interval - lower.bound.interval
  
  max.bin.fit <- max(upper.bound.interval)
  min.bin.fit <- min(lower.bound.interval)
  
  cdf.values.new <- c(0,cumsum(c(data1)))
  x.new <- c(lower.bound.interval[1], upper.bound.interval)      
  min.sigma <- max(min(diff(x.new))/2,0.25)
  
  max.bin.fit <- x.new[min(which(cdf.values.new == max(cdf.values.new)))]
  min.bin.fit <- x.new[max(which(cdf.values.new == min(cdf.values.new)))]
  
  if(substr(colnames(all.data.new)[nbr.bins+1], 1, 3)=="PDS"){# PDS
    mean.class <- c(0.75,seq(1.25, 2.75, 0.5),3.5)
    break2 <- sort(c(c(0.5,4), seq(1, 3, 0.5), seq(0.95, 2.95, by =0.5)))
  }else if(substr(colnames(all.data.new)[nbr.bins+1], 1, 3)=="SPF"){# US SPF
    mean.class <- lower.bound.interval + all.intervals/2
    break2 <- sort(c(upper.bound.interval, lower.bound.interval, lower.bound.interval[-order(lower.bound.interval)[1]] - 0.0))
    break2 <- break2[!duplicated(break2)]
  }
  
  if(area.var <= 1){
    plot.fit.survey.distribution.mixture(data1,x1,param1,
                                         min.bin   = min.bin.fit,
                                         max.bin   = max.bin.fit,
                                         min.sigma = min.sigma,
                                         mean.class,
                                         break2,
                                         xtitle = "Inflation rate %",
                                         display.mixture = display.mixture,
                                         cex = .9)
    if(display.implied){
      lines(PDF.x.all[area.var,], PDF.all[time,,horizon.y,area.var], col="grey35",
            lwd=2)
    }
    legend("topleft",legend=paste(t,"Q",as.numeric(q), sep=""), cex=1)
    legend("topright",legend=bquote(E[t]*"("*pi[t*","*t*"+"*.(horizon)]*")"),
           cex=1, bty = "n")
  } else{
    plot.fit.survey.distribution.mixture(data1,
                                         x1,
                                         param1,
                                         min.bin = min.bin.fit,
                                         max.bin = max.bin.fit,
                                         min.sigma = min.sigma,
                                         mean.class,
                                         break2,
                                         xtitle = "GDP growth rate %",
                                         display.mixture = display.mixture,
                                         cex = .9)
    if(display.implied){
      lines(PDF.x.all[area.var,], PDF.all[time,,horizon.y,area.var], col="grey35",
            lwd=2)
    }
    legend("topleft",
           legend=paste(t,"Q",as.numeric(q), sep=""), cex=1)
    legend("topright",
           legend=bquote(E[t]*"("*Delta*y[t*","*t*"+"*.(horizon)]*")"),
           cex=1, bty = "n")
  }
  
  if(display.implied & display.mixture){
    #plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
    if(display.legend){
      legend("bottom",inset = 0, title="Variables:", lty=c('blank','solid','solid'),
             c("Survey (observed)", 'Survey ("Mixture" smoothed)', 'Modeled') ,
             col = c("grey","black", "grey35"), fill=c("grey", "white", "white"), border = c("grey", "white", "white"),
             cex=0.6, horiz = TRUE)
    }
  }
  
  if(!display.implied & display.mixture){
    #plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
    if(display.legend){
      legend("bottom",inset = 0, title="Variables:", lty=c('blank','solid'),
             c("Survey (observed)", 'Survey ("Mixture" smoothed)') ,
             col = c("grey","black"), fill=c("grey", "white"), border = c("grey", "white"),
             cex=0.65, horiz = TRUE)
    }
  }
  
  if(display.implied & !display.mixture){
    #plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
    if(display.legend){
      legend("bottom",inset = 0, title="Variables:", lty=c('blank','solid'),
             c("Survey (observed)", 'Modeled') ,
             col = c("grey", "grey35"), fill=c("grey", "white"), border = c("grey", "white"),
             cex=0.65, horiz = TRUE)
    }
  }
  
  if(!display.implied & !display.mixture){
    #plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
    if(display.legend){
      legend("bottom",inset = 0, title="Variables:", lty=c('blank'),
             c("Survey (observed)") ,
             col = c("grey"), fill=c("grey"), border = c("grey"),
             cex=0.65, horiz = TRUE)
    }
  }
  
}




compute.model.implied <- function(date, type.var.implied, horizon, PDF.all, PDF.x.all){
  
  if(type.var.implied=="infl"){
    
    if(date <= "1984-10-15"){
      period <- "1981_1985"
    } else if(date <= "1991-10-15"){
      period <- "1985_1992"
    } else if(date <= "2013-10-15"){
      period <- "1992_2014"
    } else if(date > "2013-10-15"){
      period <- "post_2014"
    }
    
    area.var <- 1
    
    eval(parse(text = gsub(" ","", paste('survey.bins.distri <- US.SPF.DISTRI.',
                                         paste(horizon,sep = ""), 'Q.', paste(period,sep = ""), ' %>% left_join(survey.DATA.US.with.param %>% dplyr::select(date, contains("',
                                         paste(horizon,sep = ""), 'Q" )), by="date")',
                                         sep="")))) 
    
    eval(parse(text = gsub(" ","", paste('nbr.bins <- dim(US.SPF.DISTRI.',
                                         paste(horizon,sep = ""), 'Q.', paste(period,sep = ""),')[2]', sep="")))) 
    
  } else if(type.var.implied=="gdp"){
    
    if(date <= "1991-10-15"){
      period <- "1981_1992"
    } else if(date <= "2009-01-15"){
      period <- "1992_2009"
    } else if(date <= "2020-01-15"){
      period <- "2009_2020"
    } else if(date > "2020-01-15"){
      period <- "post_2020"
    }
    
    area.var <- 2
    
    eval(parse(text = gsub(" ","", paste('survey.bins.distri <- US.G.SPF.DISTRI.',
                                         paste(horizon,sep = ""), 'Q.', paste(period,sep = ""), ' %>% left_join(survey.DATA.US.G.with.param %>% dplyr::select(date, contains(".',
                                         paste(horizon,sep = ""), 'Q" )), by="date")',
                                         sep="")))) 
    
    eval(parse(text = gsub(" ","", paste('nbr.bins <- dim(US.G.SPF.DISTRI.',
                                         paste(horizon,sep = ""), 'Q.', paste(period,sep = ""),')[2]', sep="")))) 
    
  } else{
    stop("Choose a correct type of variable, i.e. infl or gdp")
  }
  
  # ============================================================================
  # ============================================================================
  freq <- 4
  horizon.y <- horizon-(freq-1)
  # ============================================================================
  # ============================================================================
  
  all.data.new <- survey.bins.distri[survey.bins.distri$date==as.Date(date),]
  Date <- as.Date(all.data.new[,1])
  
  if(is.na(sum(all.data.new[,2:nbr.bins]*100))){
    stop("There is no observations available. Choose a valid date")
  }
  
  data1 <- all.data.new[,2:nbr.bins]*100
  data1 <- data1[length(data1):1] # reverse order
  time <- which(Date==as.matrix(vec.dates))
  
  param1 <- as.numeric(all.data.new[,grep("param", colnames(as.data.frame(all.data.new)))])
  max.bin <- as.numeric(all.data.new[,grep("max.bin", colnames(as.data.frame(all.data.new)))])
  min.bin <- as.numeric(all.data.new[,grep("min.bin", colnames(as.data.frame(all.data.new)))])
  x1 <- PDF.x.all[area.var,]  #seq(min.bin-0.5,max.bin+0.5,by=0.05)
  
  q <- round(as.numeric(format(Date, "%m"))/3) + 1
  t <- format(Date, "%Y")
  
  lower.bound.interval <- round(sort(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",colnames(data1))))))*2)/2
  lower.bound.interval[1] <- min.bin
  upper.bound.interval <- round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(data1)))))*2)/2 
  upper.bound.interval[length(upper.bound.interval)] <- max.bin
  all.intervals <- upper.bound.interval - lower.bound.interval
  
  max.bin.fit <- max(upper.bound.interval)
  min.bin.fit <- min(lower.bound.interval)
  
  cdf.values.new <- c(0,cumsum(c(data1)))
  x.new <- c(lower.bound.interval[1], upper.bound.interval)      
  min.sigma <- max(min(diff(x.new))/2,0.25)
  
  max.bin.fit <- x.new[min(which(cdf.values.new == max(cdf.values.new)))]
  min.bin.fit <- x.new[max(which(cdf.values.new == min(cdf.values.new)))]
  
  pdf.mixture <- pdf.mixture(x1,param1,
                     min.bin   = min.bin.fit,
                     max.bin   = max.bin.fit,
                     min.sigma = min.sigma)
  cdf.mixture <- cdf.mixture(x1,param1,
                             min.bin   = min.bin.fit,
                             max.bin   = max.bin.fit,
                             min.sigma = min.sigma)
    
  x.implied <- PDF.x.all[area.var,]  
  pdf.implied <- PDF.all[time,,horizon.y,area.var]
  cdf.implied <- cumsum(pdf.implied)/sum(pdf.implied)
  
  return(list(
    x.implied = x.implied,
    pdf.implied = pdf.implied,
    cdf.implied = cdf.implied,
    x.mixture = x1,
    pdf.mixture = pdf.mixture,
    cdf.mixture = cdf.mixture
  ))
  
}

# Function to calculate the Total Variation Distance
# This is the absolute area between the two curves
total_variation_distance <- function(p, q) {
  sum(0.5*abs(p - q), na.rm = T)
}

# Function to calculate the  Wasserstein distance
# This is the x distance between two cdf
wasserstein_distance <- function(cdf1, cdf2, bins){
  return(sum(abs(cdf1 - cdf2)*diff(bins)))
}

# Function to calculate the Kullback-Leibler divergence between two probability distributions
KL_divergence <- function(p, q) {
  q[(q==0)] <- NA
  log_p_q <- log(p / q)
  log_p_q[(log_p_q==Inf | log_p_q==-Inf)] <- NA
  return(sum(p * log_p_q, na.rm = TRUE))
}

hellinger_distance <- function(p, q) {
  sqrt(0.5 * sum((sqrt(p) - sqrt(q))^2, na.rm = TRUE))
}

# Kolmogorov Metric Function with Precomputed CDFs
Kolmogorov_Metric <- function(cdf1, cdf2) {
  # Compute the maximum absolute difference between the two CDFs
  return(max(abs(cdf1 - cdf2)))
}

ks_p_value <- function(D, n) {
  return(2 * exp(-2 * (D^2) * n))
}

build_data_to_consider <- function(var) {
  var_quo <- enquo(var)
  
  higher.order.moments.decomposition %>%
    dplyr::select(date, group, value = !!var_quo) %>%
    full_join(
      fitted.obs.var.k3rd %>%
        dplyr::select(date, fit.value = !!var_quo),
      by = "date"
    )
}
