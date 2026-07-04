# ============================================
#        * * * * PAPER INFLATION - PHD * * * *
# ============================================
#           
# --------------------------------------------
#           Adrien Jean-Paul TSCHOPP 
# ============================================
#      GAUSSIAN MIXTURE DISTRIBUTION ESTIMATE EA
# --------------------------------------------
#                 * ALL *
# ============================================

#===========================================================================
# This program looks for parameterizations of incomplete Beta distributions that
# yield the same CDF as SPFs
#
# This program has to be run after load_data_US.R
#===========================================================================

dates_nny_to_keep <- unique((SPF.individual.PE %>% filter(horizon>2 & horizon <=3  & integer.horiz==FALSE))$date)

quarters_nny_to_keep <- paste0(year(dates_nny_to_keep), "Q", quarter(dates_nny_to_keep))
#quarters_nny_to_drop <- yq.infl[!(yq.infl %in% quarters_nny_to_keep)]
indices_to_keep_nny <- (yq.infl %in% quarters_nny_to_keep)

quarters_5y_to_drop <- c("1999Q2","1999Q3","1999Q4",
                         "2000Q2","2000Q3","2000Q4")
indices_to_keep <- !(yq.gdp %in% quarters_5y_to_drop)

# These are the names of the surveys used:
names.all.surveys <- c( ########## Modify
  paste0("EA.G.SPF.DISTRI.cy.",  yq.gdp,".avg"),
  paste0("EA.G.SPF.DISTRI.ny.",  yq.gdp,".avg"),
  paste0("EA.G.SPF.DISTRI.nny.",  yq.gdp[indices_to_keep_nny],".avg"),
  paste0("EA.G.SPF.DISTRI.5y.",  yq.gdp[indices_to_keep],".avg")
)

# The following names are going to be used as prefix in the 
# big DATA base.
short.names.all.surveys <- c( ########## Modify
  rep("SPF.EA.G.cy.beta",length(yq.gdp)),
  rep("SPF.EA.G.ny.beta",length(yq.gdp)),
  rep("SPF.EA.G.nny.beta",length(yq.gdp[indices_to_keep_nny])),
  rep("SPF.EA.G.5y.beta",length(yq.gdp[indices_to_keep]))
)

# Define the horizon to be considered
horizon.spf <- c(str_extract(names.all.surveys, "(?<=EA\\.G\\.SPF\\.DISTRI\\.)([^\\.]+)")) # SPF

# Define the horizon to be considered
horizon <- c(horizon.spf) # SPF

survey.DATA.G <- data.frame(date=unique(SPF.individual.PE.1y$date)) # This is a dataframe that will contain the pe and stdv series (eventually, it will be merged with DATA)
survey.DATA.G.with.param <- data.frame(date=unique(SPF.individual.PE.1y$date)) # This is a dataframe that will contain the pe and stdv series (eventually, it will be merged with DATA)

# Indicator to know if the study should be reversed
indic.US.SPF <- rep(0,length(names.all.surveys))
indic.US.SPF[grepl("SPF.US", names.all.surveys)] <- 1 # The US SPF are in reverse order

start.bin <- 4

# Initialize the count to 0
count.survey <- 0
count.by.survey <- 0

# Loop to compute the parameters for all the beta distributions
# Loop for each survey (horzion)
for(survey in names.all.surveys){
  
  # Increment count.survey by 1
  count.survey <- count.survey + 1
  if(count.survey > 1){
    if(horizon[count.survey] != horizon[count.survey-1]){
      count.by.survey <- 1
    } else{count.by.survey <- count.by.survey+1}
  } else{count.by.survey <- count.by.survey+1}
  
  # Indicate which survey the computer is working on and the number of bins
  print(paste(survey,": working on fitting beta distributions",sep=""))
  print(dim(survey.DATA.G))
  
  # Create survey.bins.distri that contains the considered survey
  eval(parse(text = gsub(" ","",paste("survey.bins.distri <- ",
                                      survey, sep=""))))
  
  # Matrices storing moments and parameters
  matrix.res <- matrix(NA,dim(survey.bins.distri)[1],7) # will contain averages, stdev, variance, skewness, cumulant of order 3, kurtosis and cumulant of order 4 of beta distributions
  matrix.param <- matrix(NA,dim(survey.bins.distri)[1],5) # will contain param of beta distributions
  #matrix.max.min <- matrix(NA,dim(survey.bins.distri)[1],2) #
  matrix.max.min.bin <- matrix(NA,dim(survey.bins.distri)[1],2) #
  
  # Loop for each time period
  i <- 1
  #for(i in 1:dim(survey.bins.distri)[1]){
  
  # keep considered time period and remove data for computation
  y <- survey.bins.distri[i,1]
  indic_id <- survey.bins.distri[i,2]
  indic_PE <- survey.bins.distri[i,3]
  distri <- survey.bins.distri[i,start.bin:dim(survey.bins.distri)[2]]
  
  # Condition if NA
  if(sum(distri,na.rm=TRUE)>.9){
    
    # Define x (bins length) 
    if(grepl("SPF", survey)){# EA SPF
      
      # define lower bound, upper bound and mean class
      lower.bound.interval.initial <- round(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",colnames(distri)))))*2)/2
      upper.bound.interval.initial <- round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(distri)))))*2)/2 
      all.intervals.initial <- upper.bound.interval.initial - lower.bound.interval.initial
      med.classes.US.initial <- lower.bound.interval.initial + all.intervals.initial/2
      
      
      # Extract the right forecast horizon
      numeric_vec <- as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", horizon[count.survey]))
      forecast_horizon <- ifelse(grepl("Q", horizon[count.survey]), ceiling(numeric_vec / 4), numeric_vec)
      # Extract PE data
      
      data.PE <- indic_PE
      
      # Point estimate based on the bins
      PE.estimate <- sum((distri)*(med.classes.US.initial))
      PE.estimate.upper <- sum((distri)*(upper.bound.interval.initial))
      PE.estimate.lower <- sum((distri)*(lower.bound.interval.initial))
      
      VAR.estimate.lower <- sum((distri)*(lower.bound.interval.initial)^2) - PE.estimate.lower^2 
      VAR.estimate.upper <- sum((distri)*(upper.bound.interval.initial)^2) - PE.estimate.upper^2 
      
      # Modify first or last bin if one of the two if > 15% (take the larger of the two)
      if(max(distri[1], distri[length(distri)])>0.1 & !is.na(data.PE)){
        
        print(data.PE)
        # Identify if we need to change the first or the last bins
        max.prob <- max(distri[1], distri[length(distri)])
        indic.pos.max.prob <- match(max.prob, distri)
        
        # Compute the interval that allows to fit perfectly the PE
        new.min.max.interval <- (data.PE - sum(((distri)*(med.classes.US.initial))[-indic.pos.max.prob]))/max.prob
        
        if(indic.pos.max.prob <= 1){ #Condition for smallest bins !!!! WARNING CHANGE THAT FOR US !!!!  
          
          # Compute the value 
          new.min.max <- 2*new.min.max.interval - upper.bound.interval.initial[indic.pos.max.prob]
          
          # Compute the value such that the interval should at least one time the interval already defined and max 4 interval already defined
          min.max.initial <- lower.bound.interval.initial[indic.pos.max.prob]
          new.min.max.corr <- min(min.max.initial, max(upper.bound.interval.initial[indic.pos.max.prob]-min(5*all.intervals.initial[indic.pos.max.prob],5),new.min.max))
          print(c(min.max.initial,new.min.max.corr))
          
          # Compute med.classes.US taking into account the new value
          lower.bound.interval <- lower.bound.interval.initial
          lower.bound.interval[indic.pos.max.prob] <- new.min.max.corr 
          upper.bound.interval <- upper.bound.interval.initial
          
          PE.estimate.lower <- sum((distri)*(lower.bound.interval))
          
          # Burst the last bin if it's bigger than the next-to-last one
          increment.interval <- upper.bound.interval-lower.bound.interval
          if(distri[length(distri)]/increment.interval[length(distri)] > distri[length(distri)-1]/increment.interval[length(distri)-1]){
            print("Change max")
            
            new.increment.interval <- as.numeric(distri[length(distri)]/(distri[length(distri)-1]/increment.interval[length(distri)-1]))
            upper.bound.interval[length(distri)] <- lower.bound.interval[length(distri)] + new.increment.interval
            
          }
          
          x <- sort(upper.bound.interval)
          
        } else{ #Condition for largest bins  
          
          # Compute the value 
          new.min.max <- 2*new.min.max.interval - lower.bound.interval.initial[indic.pos.max.prob]
          
          # Compute the value such that the interval should at least one time the interval already defined and max 5 interval already defined or 5 or new.min.max
          min.max.initial <- upper.bound.interval.initial[indic.pos.max.prob]
          new.min.max.corr <- max(min.max.initial, min(lower.bound.interval.initial[indic.pos.max.prob]+min(5*all.intervals.initial[indic.pos.max.prob],5),new.min.max))
          print(c(min.max.initial,new.min.max.corr))
          
          # Compute med.classes.US taking into account the new value
          upper.bound.interval <- upper.bound.interval.initial
          upper.bound.interval[indic.pos.max.prob] <- new.min.max.corr
          lower.bound.interval <- lower.bound.interval.initial
          
          PE.estimate.upper <- sum((distri)*(upper.bound.interval))
          
          # Burst the first bin if it's bigger than the second one
          increment.interval <- upper.bound.interval-lower.bound.interval
          if(distri[1]/increment.interval[1] > distri[2]/increment.interval[1]){
            print("Change min")
            
            new.increment.interval <- as.numeric(distri[1]/(distri[2]/increment.interval[2]))
            lower.bound.interval[1] <- upper.bound.interval[1] - new.increment.interval
            
          }
          
          x <- sort(upper.bound.interval)
          
        }
        
      } else{ #Condition to consider initial value
        upper.bound.interval <- upper.bound.interval.initial
        lower.bound.interval <- lower.bound.interval.initial
        x <- sort(round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(survey.bins.distri)[-c(1:(start.bin-1))]))))*2)/2) 
      }
      
      # If there is too much difference between the two measure
      if(!is.na(data.PE) & (PE.estimate.upper*0.975<data.PE | PE.estimate.lower*1.025>data.PE)){
        data.PE<-NA
      }
      
    }
    
    if(indic.US.SPF[count.survey]==1){# reverse order
      distri <- distri[length(distri):1]
    }
    
    # Define cdf.values (sum of probabilities)
    #x <- seq(-1.5,8,by=.5)
    cdf.values <- cumsum(c(distri))
    
    # Account for data transformation
    # if((indic.US.SPF[count.survey]==1)&
    #    (as.integer(format(survey.bins.distri[i,1],"%Y"))<2014)
    # ){# In that case, the US data had been divided on two bins
    #   x <- x[seq(1,length(x),by=2)]
    #   cdf.values <- cdf.values[seq(1,length(cdf.values),by=2)]
    # }
    
    # Check if potential problems with dimension
    if(length(cdf.values)!=length(x)){
      print("pbm ici")
      stop()
    }
    
    # Compute the results 
    #res.fit <- fit.cdf(x,cdf.values,param.0=c(5,5,-1,6))
    #### For individuals who only have one bin use triangular distribution
    
    # Add min of lower intervall in x 
    # Pdf value is 0 for this x
    # Cdf value is 0 for this x
    x <- c(min(lower.bound.interval),x)
    cdf.values <- c(0,cdf.values)
    pdf.values <- c(0,c(distri))
    min.sigma <- max(min(diff(x))/2,0.25) #min(min(diff(x))/2)
    
    max.bin <- max(upper.bound.interval)
    min.bin <- min(lower.bound.interval)
    
    max.bin.fit <- x[min(which(cdf.values == max(cdf.values)))]
    min.bin.fit <- x[max(which(cdf.values == min(cdf.values)))]
    
    res.fit <- fit.cdf.mixture.PE(x,cdf.values,param.0=c(0,.1,.2,0,3),min.bin=min.bin.fit,max.bin=max.bin.fit,min.sigma=min.sigma, PE=NA)
    moments <- moments.mixture.PE(res.fit$param,min.bin=min.bin,max.bin=max.bin,min.sigma=min.sigma,PE=NA)
    
    if(!is.na(data.PE) & abs(moments$Mean/data.PE-1)>0.1){
      #res.fit <- fit.cdf.mixture.PE(x,cdf.values,param.0=c(0,.1,.2,0,3),min.bin=min.bin,max.bin=max.bin,min.sigma=min.sigma,PE=data.PE)
      #moments <- moments.mixture.PE(res.fit$param,min.bin=min.bin,max.bin=max.bin,min.sigma=min.sigma,PE=data.PE)
    }
    
    cdf.values.used <- cdf.values
    x.used <- x
    #print(moments)
    
    # Remove first value of x
    x <- x[-1]
    cdf.values <- cdf.values[-1]
    pdf.values <- pdf.values[-1] 
    
    # Store the results in the matrices
    matrix.res[i,1] <- moments$Mean
    matrix.res[i,2] <- sqrt(moments$Variance)
    matrix.res[i,3] <- moments$Variance
    matrix.res[i,4] <- moments$Skewness
    matrix.res[i,5] <- make.cumulant.until.order.4(moments$Mean,moments$Variance,moments$Skewness,moments$Kurtosis)$k.3
    matrix.res[i,6] <- moments$Kurtosis
    matrix.res[i,7] <- make.cumulant.until.order.4(moments$Mean,moments$Variance,moments$Skewness,moments$Kurtosis)$k.4
    
    matrix.param[i,] <- res.fit$param
    #matrix.max.min[i,] <- c(res.fit$max.c,res.fit$min.d) #
    matrix.max.min.bin[i,] <- c(max(upper.bound.interval),min(lower.bound.interval))
    
  }
  #}
  
  # Store Moments of the incomplete beta distribution.
  name.pe.series <- paste(short.names.all.surveys[count.survey],".pe", sep="")
  name.stdev.series <- paste(short.names.all.surveys[count.survey],".stdev", sep="")
  name.var.series <- paste(short.names.all.surveys[count.survey],".var", sep="")
  name.skew.series <- paste(short.names.all.surveys[count.survey],".skew", sep="")
  name.k3rd.series <- paste(short.names.all.surveys[count.survey],".k3rd", sep="")
  name.kurtosis.series <- paste(short.names.all.surveys[count.survey],".kurtosis", sep="")
  name.k4th.series <- paste(short.names.all.surveys[count.survey],".k4th", sep="")
  
  # Store Parameterization of the incomplete beta distribution.
  name.beta.param1.series <- paste(short.names.all.surveys[count.survey],".param1", sep="")
  name.beta.param2.series <- paste(short.names.all.surveys[count.survey],".param2", sep="")
  name.beta.param3.series <- paste(short.names.all.surveys[count.survey],".param3", sep="")
  name.beta.param4.series <- paste(short.names.all.surveys[count.survey],".param4", sep="")
  name.beta.param5.series <- paste(short.names.all.surveys[count.survey],".param5", sep="")
  
  
  # Store Parameterization of the incomplete beta distribution.
  #name.max.c.series <- paste(short.names.all.surveys[count.survey],".max.c", sep="")
  #name.min.d.series <- paste(short.names.all.surveys[count.survey],".min.d", sep="")
  
  name.max.bin.series <- paste(short.names.all.surveys[count.survey],".max.bin", sep="")
  name.min.bin.series <- paste(short.names.all.surveys[count.survey],".min.bin", sep="")
  
  # Add results to survey.bin.distri
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.pe.series,
                                      "<- matrix.res[,1]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.stdev.series,
                                      "<- matrix.res[,2]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.var.series,
                                      "<- matrix.res[,3]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.skew.series,
                                      "<- matrix.res[,4]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.k3rd.series,
                                      "<- matrix.res[,5]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.kurtosis.series,
                                      "<- matrix.res[,6]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.k4th.series,
                                      "<- matrix.res[,7]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.beta.param1.series,
                                      "<- matrix.param[,1]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.beta.param2.series,
                                      "<- matrix.param[,2]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.beta.param3.series,
                                      "<- matrix.param[,3]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.beta.param4.series,
                                      "<- matrix.param[,4]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.beta.param5.series,
                                      "<- matrix.param[,5]",sep=""))))
  # eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.max.c.series,
  #                                     "<- matrix.max.min[,1]",sep=""))))
  # eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.min.d.series,
  #                                     "<- matrix.max.min[,2]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.max.bin.series,
                                      "<- matrix.max.min.bin[,1]",sep=""))))
  eval(parse(text = gsub(" ","",paste("survey.bins.distri$",name.min.bin.series,
                                      "<- matrix.max.min.bin[,2]",sep=""))))
  
  
  # Remove dates to survey.bins.distri
  eval(parse(text = gsub(" ","",paste('survey.stats <- subset(survey.bins.distri,select = c("date","',
                                      name.pe.series,'","',
                                      name.stdev.series,'","',
                                      name.var.series,'","',
                                      name.skew.series,'","',
                                      name.k3rd.series,'","',
                                      name.kurtosis.series,'","',
                                      name.k4th.series,'"))',sep=""))))
  
  eval(parse(text = gsub(" ","",paste('survey.stats.with.param <- subset(survey.bins.distri,select = c("date","',
                                      name.pe.series,'","',
                                      name.stdev.series,'","',
                                      name.var.series,'","',
                                      name.skew.series,'","',
                                      name.k3rd.series,'","',
                                      name.kurtosis.series,'","',
                                      name.k4th.series,'","',
                                      name.beta.param1.series,'","',
                                      name.beta.param2.series,'","',
                                      name.beta.param3.series,'","',
                                      name.beta.param4.series,'","',
                                      name.beta.param5.series,'","',
                                      #name.max.c.series,'","',
                                      #name.min.d.series,'","',
                                      name.max.bin.series,'","',
                                      name.min.bin.series,'"))',sep=""))))
  
  
  # Merge all data
  survey.DATA.G            <- survey.DATA.G  %>%
    full_join(survey.stats, by=c("date"))
  
  ## Extract variables finishing by .x if the column already existed
  var_names <- survey.DATA.G %>%
    dplyr::select(ends_with(".x")) %>%
    names() %>%
    sub("\\.x$", "", .)
  
  ## combine variable finishing by .x and .y and keep only one removing .x
  for (col in var_names) {
    cols_x <- paste0(col, ".x")
    cols_y <- paste0(col, ".y")
    survey.DATA.G <- survey.DATA.G %>% mutate(!!paste0(col) := coalesce(!!sym(cols_x), !!sym(cols_y)))
  }
  survey.DATA.G <- survey.DATA.G %>% dplyr::select(-ends_with(".x"), -ends_with(".y"))
  
  # Merge all data
  survey.DATA.G.with.param            <- survey.DATA.G.with.param  %>%
    full_join(survey.stats.with.param, by=c("date"))
  
  ## Extract variables finishing by .x if the column already existed
  var_names <- survey.DATA.G.with.param %>%
    dplyr::select(ends_with(".x")) %>%
    names() %>%
    sub("\\.x$", "", .)
  
  ## combine variable finishing by .x and .y and keep only one removing .x
  for (col in var_names) {
    cols_x <- paste0(col, ".x")
    cols_y <- paste0(col, ".y")
    survey.DATA.G.with.param <- survey.DATA.G.with.param %>% mutate(!!paste0(col) := coalesce(!!sym(cols_x), !!sym(cols_y)))
  }
  survey.DATA.G.with.param <- survey.DATA.G.with.param %>% dplyr::select(-ends_with(".x"), -ends_with(".y"))
  
  # Plots evolutions of the moments for all time periods
  ## Parameters in order to diplay well the chart
  # r <- 2
  # m <- matrix(seq(1,r+1),nrow = (r+1),ncol = 1,byrow = TRUE)
  # layout(mat = m, heights = c(rep(0.8/r,r),0.2))
  # par(mar=c(2, 3, 1, 1))
  
  # ============================================  
  # Plot the distribution for the last 6 observations.
  
  #survey.bins.distri.without.na <- na.omit(survey.bins.distri) 
  survey.bins.distri.without.na <- survey.bins.distri[rowSums(is.na(survey.bins.distri[ , 3:dim(survey.bins.distri)[2]])) <= 10, ]
  max.m <- dim(survey.bins.distri.without.na)[1]
  #seq.m <- seq(max.m-max(5,max.m-1),max.m, by=1)
  seq.m <- seq(max.m-min(5,max.m-1),max.m, by=1)
  
  if(count.by.survey==1 | ((count.by.survey-1) %% 6)==0){
    ## Parameters in order to diplay well the chart
    m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
    #m <- matrix(c(1,2,3,4,5,6),nrow = 2,ncol = 3,byrow = TRUE)
    
    
    layout(mat = m, heights = c(0.84/2,0.84/2,0.16))
    
    par(mar=c(4, 4.5, 1, 1))
  }
  #for (m in seq.m) {
  m <- 1
  l.survey <- length(x)
  
  all.data.new <- survey.bins.distri.without.na[m,]
  Date <- all.data.new[,1]
  ID <- all.data.new[,2]
  data1 <- all.data.new[,start.bin:(l.survey+start.bin-1)]*100
  
  if(indic.US.SPF[count.survey]==1){# reverse order
    data1 <- data1[length(data1):1]
  }
  
  param1 <- as.numeric(all.data.new[,grep("param", colnames(as.data.frame(all.data.new)))])
  #max.c1 <- as.numeric(all.data.new[,grep("max.c", colnames(as.data.frame(all.data.new)))])
  #min.d1 <- as.numeric(all.data.new[,grep("min.d", colnames(as.data.frame(all.data.new)))])
  max.bin <- as.numeric(all.data.new[,grep("max.bin", colnames(as.data.frame(all.data.new)))])
  min.bin <- as.numeric(all.data.new[,grep("min.bin", colnames(as.data.frame(all.data.new)))])
  #x1 <- seq(-2,7,by=0.05)
  #x1 <- seq(min(x)-2,max(x),by=0.05)
  x1 <- seq(min.bin-0.5,max.bin+0.5,by=0.05)
  
  q <- format(Date, "%m")
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
  
  if(substr(short.names.all.surveys[count.survey], 1, 3)=="SPF"){# US SPF
    #mean.class <- c(-1,seq(0.25, 7.75, 0.5),9)
    #break2 <- sort(c(c(-2,10), seq(0, 8, 0.5), seq(0.45, 7.95, by =0.5)))
    
    mean.class <- lower.bound.interval + all.intervals/2
    #break2 <- sort(c(upper.bound.interval, lower.bound.interval, lower.bound.interval[-order(lower.bound.interval)[1:2]] - 0.05))
    break2 <- sort(c(upper.bound.interval, lower.bound.interval, lower.bound.interval[-order(lower.bound.interval)[1]] - 0.00))
    break2 <- break2[!duplicated(break2)]
  }
  
  
  plot.fit.survey.distribution.mixture(data1,x1,param1,min.bin=min.bin.fit,max.bin=max.bin.fit,min.sigma=min.sigma, mean.class, break2, xtitle="GDP growth rate %")
  legend("topleft",legend=paste(t,":",q, sep=""), cex=0.6)
  legend("topright",legend=paste(horizon[count.survey], sep=""), cex=0.6)
  
  
  #}
  
  if(count.by.survey!= 1 & (count.by.survey %% 6)==0){
    par(mar=c(1, 4, 1, 1))
    plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
    legend("bottom",inset = 0, title="Variables:", lty=c('blank','solid'),
           c("Survey (observed)", 'Survey ("Beta" smoothed)') ,
           col = c("grey","black"), fill=c("grey", "white"), border = c("grey", "white"),
           cex=1, horiz = TRUE)}
  
  
}



all_objects <- ls()
# Find all objects that start with "EA.SPF.DISTRI"
data.to.delete <- all_objects[grep("^EA\\.G\\.SPF\\.DISTRI", all_objects)]
# Exclude objects that end with ".avg"
data.to.delete <- data.to.delete[!grepl("\\.avg$", data.to.delete)]
#data.to.delete.1 <- all_objects[grep("^EA.G.SPF.DISTRI.1y", all_objects)]
#data.to.delete.2 <- all_objects[grep("^EA.G.SPF.DISTRI.2y", all_objects)]
#data.to.delete.5 <- all_objects[grep("^EA.G.SPF.DISTRI.5y", all_objects)]

rm(list = data.to.delete)
#rm(list = data.to.delete.1)
#rm(list = data.to.delete.2)
#rm(list = data.to.delete.5)


