# ============================================
#        * * * * MRT * * * *
# ============================================
#           
# --------------------------------------------
#           Adrien Jean-Paul TSCHOPP 
# ============================================
# GAUSSIAN MIXTURE DISTRIBUTION ESTIMATE US GROWTH
# --------------------------------------------
#                 * ALL *
# ============================================

# 
# This program looks for parameterizations of incomplete Beta distributions that
# yield the same CDF as SPFs
#
# This program has to be run after load.data.R
#

if (!exists("save_mixture_diagnostic_plots")) {
  save_mixture_diagnostic_plots <- FALSE
}
if (isTRUE(save_mixture_diagnostic_plots)) {
  dir.create("graphs/mixture_diagnostics", recursive = TRUE, showWarnings = FALSE)
  pdf("graphs/mixture_diagnostics/gdp_distribution_smoothing_diagnostics.pdf")
} else {
  pdf(NULL)
}

# Define the different dates at which the bins changed
dates.all.surveys <- c(
  "1981_1992",
  "1992_2009",
  "2009_2020",
  "post_2020"
  )

# Define the horizon to be considered
horizon <- c(paste0(sort(c(rep(seq(1,8,1),length(dates.all.surveys)),rep(seq(9,16,1),length(dates.all.surveys[3:4])))),"Q", sep="")) # SPF


# These are the names of the surveys used:
names.all.surveys <- c(
  paste0("US.G.SPF.DISTRI.1Q", ".", dates.all.surveys),
  paste0("US.G.SPF.DISTRI.2Q", ".", dates.all.surveys),
  paste0("US.G.SPF.DISTRI.3Q", ".", dates.all.surveys),
  paste0("US.G.SPF.DISTRI.4Q", ".", dates.all.surveys),
  paste0("US.G.SPF.DISTRI.5Q", ".", dates.all.surveys),
  paste0("US.G.SPF.DISTRI.6Q", ".", dates.all.surveys),
  paste0("US.G.SPF.DISTRI.7Q", ".", dates.all.surveys),
  paste0("US.G.SPF.DISTRI.8Q", ".", dates.all.surveys),
  paste0("US.G.SPF.DISTRI.9Q", ".", dates.all.surveys[3:4]),
  paste0("US.G.SPF.DISTRI.10Q", ".", dates.all.surveys[3:4]),
  paste0("US.G.SPF.DISTRI.11Q", ".", dates.all.surveys[3:4]),
  paste0("US.G.SPF.DISTRI.12Q", ".", dates.all.surveys[3:4]),
  paste0("US.G.SPF.DISTRI.13Q", ".", dates.all.surveys[3:4]),
  paste0("US.G.SPF.DISTRI.14Q", ".", dates.all.surveys[3:4]),
  paste0("US.G.SPF.DISTRI.15Q", ".", dates.all.surveys[3:4]),
  paste0("US.G.SPF.DISTRI.16Q", ".", dates.all.surveys[3:4])
)

# The following names are going to be used as prefix in the 
# big DATA base.
short.names.all.surveys <- c(
  rep("SPF.US.G.1Q.beta", length(dates.all.surveys)),
  rep("SPF.US.G.2Q.beta", length(dates.all.surveys)),
  rep("SPF.US.G.3Q.beta", length(dates.all.surveys)),
  rep("SPF.US.G.4Q.beta", length(dates.all.surveys)),
  rep("SPF.US.G.5Q.beta", length(dates.all.surveys)),
  rep("SPF.US.G.6Q.beta", length(dates.all.surveys)),
  rep("SPF.US.G.7Q.beta", length(dates.all.surveys)),
  rep("SPF.US.G.8Q.beta", length(dates.all.surveys)),
  rep("SPF.US.G.9Q.beta", length(dates.all.surveys[3:4])),
  rep("SPF.US.G.10Q.beta", length(dates.all.surveys[3:4])),
  rep("SPF.US.G.11Q.beta", length(dates.all.surveys[3:4])),
  rep("SPF.US.G.12Q.beta", length(dates.all.surveys[3:4])),
  rep("SPF.US.G.13Q.beta", length(dates.all.surveys[3:4])),
  rep("SPF.US.G.14Q.beta", length(dates.all.surveys[3:4])),
  rep("SPF.US.G.15Q.beta", length(dates.all.surveys[3:4])),
  rep("SPF.US.G.16Q.beta", length(dates.all.surveys[3:4]))
)

survey.DATA.US.G <- data.frame(date=DATA.G.US$date) # This is a dataframe that will contain the pe and stdv series (eventually, it will be merged with DATA)
survey.DATA.US.G.with.param <- data.frame(date=DATA.G.US$date) # This is a dataframe that will contain the pe and stdv series (eventually, it will be merged with DATA)

# Indicator to know if the study should be reversed
indic.US.G.SPF <- rep(0,length(names.all.surveys))
indic.US.G.SPF[grepl("SPF", names.all.surveys)] <- 1 # The US SPF are in reverse order
#indic.US.SPF[5:length(indic.US.SPF)] <- 1 # The US SPF are in reverse order


# Initialize the count to 0
count.survey <- 0

# Loop to compute the parameters for all the beta distributions
# Loop for each survey (horzion)
for(survey in names.all.surveys){
  
  # Increment count.survey by 1
  count.survey <- count.survey + 1
  
  # Indicate which survey the computer is working on and the number of bins
  print(paste(survey,": working on fitting beta distributions",sep=""))
  print(dim(survey.DATA.US.G))
  
  # Create survey.bins.distri that contains the considered survey
  eval(parse(text = gsub(" ","",paste("survey.bins.distri <- ",
                                      survey, sep=""))))
  
  # Matrices storing moments and parameters
  matrix.res <- matrix(NA,dim(survey.bins.distri)[1],7) # will contain averages, stdev, variance, skewness, cumulant of order 3, kurtosis and cumulant of order 4 of beta distributions
  matrix.param <- matrix(NA,dim(survey.bins.distri)[1],5) # will contain param of beta distributions
  #matrix.max.min <- matrix(NA,dim(survey.bins.distri)[1],2) #
  matrix.max.min.bin <- matrix(NA,dim(survey.bins.distri)[1],2) #
  
  # Loop for each time period
  for(i in 1:(dim(survey.bins.distri)[1])){
    
    # keep considered time period and remove data for computation
    y <- survey.bins.distri[i,1]
    distri <- survey.bins.distri[i,2:dim(survey.bins.distri)[2]]
    
    # Condition if NA
    if(sum(distri,na.rm=TRUE)>.9){
      
      # Define x (bins length) 
      #if(dim(survey.bins.distri)[2]==7){# PDS
      if(grepl("PDS", survey)){# PDS  
        lower.bound.interval <- 0.5  
        upper.bound.interval <- 7
        
        # define lower bound, upper bound and mean class
        lower.bound.interval.initial <- round(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",colnames(distri)))))*2)/2
        upper.bound.interval.initial <- round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(distri)))))*2)/2 
        all.intervals.initial <- upper.bound.interval.initial - lower.bound.interval.initial
        med.classes.US.initial <- lower.bound.interval.initial + all.intervals.initial/2
        
        x <- c(1,1.5,2,2.5,3,7)
        data.PE <- data.frame(date=y,value=NA)
        PE.estimate <- sum((distri)*(med.classes.US.initial))
        PE.estimate.upper <- sum((distri)*(upper.bound.interval.initial))
        PE.estimate.lower <- sum((distri)*(lower.bound.interval.initial))
        VAR.estimate.lower <- sum((distri)*(lower.bound.interval.initial)^2) - PE.estimate.lower^2 
        VAR.estimate.lower <- 0
        
        
        #}else if(dim(survey.bins.distri)[2]<20){# US SPF
      }else if(grepl("SPF", survey)){# US SPF
        
        # define lower bound, upper bound and mean class
        lower.bound.interval.initial <- round(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",colnames(distri)))))*2)/2
        upper.bound.interval.initial <- round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(distri)))))*2)/2 
        all.intervals.initial <- upper.bound.interval.initial - lower.bound.interval.initial
        med.classes.US.initial <- lower.bound.interval.initial + all.intervals.initial/2
        
        # Extract the right forecast horizon
        numeric_vec <- as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", horizon[count.survey]))
        forecast_horizon <- ifelse(grepl("Q", horizon[count.survey]), ceiling(numeric_vec / 4), numeric_vec)
        # Extract PE data
        data.PE <- SPF_US_G_PE_aggregate_reshape %>% filter(date %in% c(y) & (Forecast_Horizon==forecast_horizon))
        
        # Point estimate based on the bins
        PE.estimate <- sum((distri)*(med.classes.US.initial))
        PE.estimate.upper <- sum((distri)*(upper.bound.interval.initial))
        PE.estimate.lower <- sum((distri)*(lower.bound.interval.initial))
        
        VAR.estimate.lower <- sum((distri)*(lower.bound.interval.initial)^2) - PE.estimate.lower^2 
        
        # Modify first or last bin if one of the two if > 10% (take the larger of the two)
        if(max(distri[1], distri[length(distri)])>0.1 & !is.na(data.PE$value)){
          
          print(data.PE)
          # Identify if we need to change the first or the last bins
          max.prob <- max(distri[1], distri[length(distri)])
          indic.pos.max.prob <- match(max.prob, distri)
          
          # Compute the interval that allows to fit perfectly the PE
          new.min.max.interval <- (data.PE$value - sum(((distri)*(med.classes.US.initial))[-indic.pos.max.prob]))/max.prob
          
          if(indic.pos.max.prob > 1){ #Condition for smallest bins  
            
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
            VAR.estimate.lower <- sum((distri)*(lower.bound.interval)^2) - PE.estimate.lower^2 
            
            # Burst the largest bin if it's bigger than the second largest
            increment.interval <- upper.bound.interval-lower.bound.interval
            if(distri[1]/increment.interval[1] > distri[2]/increment.interval[1]){
              print("Change max")
              
              new.increment.interval <- as.numeric(distri[1]/(distri[2]/increment.interval[2]))
              lower.bound.interval[1] <- upper.bound.interval[1] - new.increment.interval
              
            }
            
            x <- sort(upper.bound.interval)
            
          } else{ #Condition for largest bins  
            
            # Compute the value 
            new.min.max <- 2*new.min.max.interval - lower.bound.interval.initial[indic.pos.max.prob]
            
            # Compute the value such that the interval should at least one time the interval already defined and max 4 interval already defined
            min.max.initial <- upper.bound.interval.initial[indic.pos.max.prob]
            new.min.max.corr <- max(min.max.initial, min(lower.bound.interval.initial[indic.pos.max.prob]+min(5*all.intervals.initial[indic.pos.max.prob],5),new.min.max))
            print(c(min.max.initial,new.min.max.corr))
            
            # Compute med.classes.US taking into account the new value
            upper.bound.interval <- upper.bound.interval.initial
            upper.bound.interval[indic.pos.max.prob] <- new.min.max.corr
            lower.bound.interval <- lower.bound.interval.initial
            
            PE.estimate.lower <- sum((distri)*(lower.bound.interval))
            VAR.estimate.lower <- sum((distri)*(lower.bound.interval)^2) - PE.estimate.lower^2 
            
            # Burst the smallest bin if it's bigger than second smallest one
            increment.interval <- upper.bound.interval-lower.bound.interval
            
            if(distri[length(distri)]/increment.interval[length(distri)] > distri[length(distri)-1]/increment.interval[length(distri)-1]){
              print("Change min")
              
              new.increment.interval <- as.numeric(distri[length(distri)]/(distri[length(distri)-1]/increment.interval[length(distri)-1]))
              upper.bound.interval[length(distri)] <- lower.bound.interval[length(distri)] + new.increment.interval
              
            }
            
            x <- sort(upper.bound.interval)
            
          }
          
        } else{ #Condition to consider initial value
          upper.bound.interval <- upper.bound.interval.initial
          lower.bound.interval <- lower.bound.interval.initial
          x <- sort(round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(survey.bins.distri)[-1]))))*2)/2) 
        }
        
        # x <- c(seq(0,8,by=.5),10)
        # x <- sort(round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(survey.bins.distri)[-1]))))*2)/2) 
        
      }
      
      # If there is too much difference between the two measure
      # if(!is.na(data.PE$value) & (PE.estimate.upper*0.99<data.PE$value | PE.estimate.lower*1.01>data.PE$value)){
      #   data.PE$value<-NA
      # }
      
      if(indic.US.G.SPF[count.survey]==1){# reverse order
        distri <- distri[length(distri):1]
      }
      
      # Define cdf.values (sum of probabilities)
      #x <- seq(-1.5,8,by=.5)
      cdf.values <- cumsum(c(distri))
      
      # Check if potential problems with dimension
      if(length(cdf.values)!=length(x)){
        print("pbm ici")
        stop()
      }
      
      # Compute the results 
      #res.fit <- fit.cdf(x,cdf.values,param.0=c(5,5,-1,6))
      #res.fit <- fit.cdf.PE(x,cdf.values,param.0=c(5,5,-1,6),PE=data.PE$value, max.c=min(lower.bound.interval))
      
      # Add min of lower intervall in x 
      # Pdf value is 0 for this x
      # Cdf value is 0 for this x
      x <- c(min(lower.bound.interval),x)
      cdf.values <- c(0,cdf.values)
      pdf.values <- c(0,c(distri))
      min.sigma <- max(min(diff(x))/2,0.25)

      print(data.PE$date)
      max.bin.fit <- max(upper.bound.interval)
      min.bin.fit <- min(lower.bound.interval)
      print("first")
      print(c(min.bin.fit,max.bin.fit))

      max.bin.fit <- x[min(which(cdf.values == max(cdf.values)))]
      min.bin.fit <- x[max(which(cdf.values == min(cdf.values)))]
      print("second")
      print(c(min.bin.fit,max.bin.fit))
      
      res.fit <- fit.cdf.mixture.PE(x,cdf.values,param.0=c(0,.1,.2,0,3),min.bin=min.bin.fit,max.bin=max.bin.fit,min.sigma=min.sigma, PE=NA)
      moments <- moments.mixture.PE(res.fit$param,min.bin=min.bin.fit,max.bin=max.bin.fit,min.sigma=min.sigma,PE=NA)
      
      if(!is.na(data.PE$value) & max(distri[1], distri[length(distri)])>0.1 & abs(moments$Mean/data.PE$value-1)>0.2){
        #res.fit <- fit.cdf.mixture.PE(x,cdf.values,param.0=c(0,.1,.2,0,3),min.bin=min.bin.fit,max.bin=max.bin.fit,min.sigma=min.sigma,PE=data.PE$value)
        #moments <- moments.mixture.PE(res.fit$param,min.bin=min.bin.fit,max.bin=max.bin.fit,min.sigma=min.sigma,PE=data.PE$value)
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
  }
  
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
  survey.DATA.US.G            <- merge(survey.DATA.US.G,survey.stats,
                                       by="date",all=TRUE)
  
  ## Extract variables finishing by .x if the column already existed
  var_names <- survey.DATA.US.G %>%
    dplyr::select(ends_with(".x")) %>%
    names() %>%
    sub("\\.x$", "", .)
  
  ## combine variable finishing by .x and .y and keep only one removing .x
  for (col in var_names) {
    cols_x <- paste0(col, ".x")
    cols_y <- paste0(col, ".y")
    survey.DATA.US.G <- survey.DATA.US.G %>% mutate(!!paste0(col) := coalesce(!!sym(cols_x), !!sym(cols_y)))
  }
  survey.DATA.US.G <- survey.DATA.US.G %>% dplyr::select(-ends_with(".x"), -ends_with(".y"))
  
  # Merge all data
  survey.DATA.US.G.with.param <- merge(survey.DATA.US.G.with.param,survey.stats.with.param,
                                       by="date",all=TRUE)
  ## Extract variables finishing by .x if the column already existed
  var_names <- survey.DATA.US.G.with.param %>%
    dplyr::select(ends_with(".x")) %>%
    names() %>%
    sub("\\.x$", "", .)
  
  ## combine variable finishing by .x and .y and keep only one removing .x
  for (col in var_names) {
    cols_x <- paste0(col, ".x")
    cols_y <- paste0(col, ".y")
    survey.DATA.US.G.with.param <- survey.DATA.US.G.with.param %>% mutate(!!paste0(col) := coalesce(!!sym(cols_x), !!sym(cols_y)))
  }
  survey.DATA.US.G.with.param <- survey.DATA.US.G.with.param %>% dplyr::select(-ends_with(".x"), -ends_with(".y"))
  
  
  # Plots evolutions of the moments for all time periods
  ## Parameters in order to diplay well the chart
  r <- 2
  m <- matrix(seq(1,r+1),nrow = (r+1),ncol = 1,byrow = TRUE)
  layout(mat = m, heights = c(rep(0.8/r,r),0.2))
  par(mar=c(2, 3, 1, 1))
  
  # extract observed data of the considered survey
  if(substr(short.names.all.surveys[count.survey], 1, 3)=="PDS"){# PDS
    
    if(grepl("avg",short.names.all.surveys[count.survey])){
      eval(parse(text = gsub(" ","", paste('survey.date  <- DATA.G.US$date', sep="")))) 
      eval(parse(text = gsub(" ","", paste('survey.expectation  <- DATA.G.US$PDS.US.G.',
                                           horizon[count.survey], '.pe',sep="")))) 
      eval(parse(text = gsub(" ","", paste('survey.standard.deviation <- DATA.US$PDS.US.G.',
                                           paste(horizon[count.survey], '.stdv' , sep = ""),
                                           sep="")))) 
    } else {
      eval(parse(text = gsub(" ","", paste('survey.date  <- US.G.PDS.data$date', sep="")))) 
      eval(parse(text = gsub(" ","", paste('survey.expectation  <- US.G.PDS.data$PDS.pe.',
                                           horizon[count.survey], sep="")))) 
      eval(parse(text = gsub(" ","", paste('survey.standard.deviation <- US.G.PDS.data$PDS.disagreement.',
                                           paste(horizon[count.survey],sep = ""),
                                           sep=""))))
    }
  }else if(substr(short.names.all.surveys[count.survey], 1, 3)=="SPF"){# US SPF
    eval(parse(text = gsub(" ","", paste('survey.date  <- US.G.SPF.data$date', sep="")))) 
    eval(parse(text = gsub(" ","", paste('survey.expectation  <- US.G.SPF.data$SPF.US.G.pe.',
                                         horizon[count.survey], sep="")))) 
    eval(parse(text = gsub(" ","", paste('survey.standard.deviation <- US.G.SPF.data$SPF.US.G.disagreement.',
                                         horizon[count.survey], sep="")))) 
  }
  
  
  # ============================================
  #                      PLOTS 
  # ============================================
  
  # ============================================
  # Plot time series
  
  if(grepl("post", survey)){
    
    indic.serie <- colnames(survey.DATA.US.G) == paste(short.names.all.surveys[count.survey],".pe",sep="")
    # Plot Evolution of each survey
    ## MEAN AND VARIANCE
    plot(survey.DATA.US.G$date,survey.DATA.US.G[,indic.serie], xlab="Time", ylab="Inflation rate (%)",
         col="blue", lwd="2", main=names.all.surveys[count.survey],
         ylim = c(min(survey.DATA.US.G[,indic.serie], na.rm = T) -0.15, max(survey.DATA.US.G[,indic.serie], na.rm = T)+0.15))
    lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G[,indic.serie])],survey.DATA.US.G[,indic.serie][!is.na(survey.DATA.US.G[,indic.serie])],
          col="blue", lwd="2")
    points(survey.date, survey.expectation, pch=19)
    ## AND VARIANCE
    indic.serie <- colnames(survey.DATA.US.G) == paste(short.names.all.surveys[count.survey],".var",sep="")
    plot(survey.DATA.US.G$date,survey.DATA.US.G[,indic.serie], xlab="Time", ylab="Inflation rate (%)",
         col="red", lwd="2",
         ylim = c(min(survey.DATA.US.G[,indic.serie], na.rm = T) -0.25, max(survey.DATA.US.G[,indic.serie], na.rm = T)+0.25))
    lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G[,indic.serie])],survey.DATA.US.G[,indic.serie][!is.na(survey.DATA.US.G[,indic.serie])],
          col="red", lwd="2")
    points(survey.date, survey.standard.deviation^2, pch=19)
    
    par(mar=c(0.1, 4, 1, 1))
    plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
    legend("bottom",inset = 0, title="Variable(s):",
           c("Expected Probability", "Variance", "Observed Moment"),
           fill=c("blue", "red", "black"), cex=1.1, horiz = TRUE)
    
    par(mar=c(2, 3, 1, 1))
    
    ## SKEWNESS AND CUMULANT OF ORDER 3
    indic.serie <- colnames(survey.DATA.US.G) == paste(short.names.all.surveys[count.survey],".skew",sep="")
    plot(survey.DATA.US.G$date,survey.DATA.US.G[,indic.serie], xlab="Time", ylab="Inflation rate (%)",
         col="blue", lwd="2", main=names.all.surveys[count.survey],
         ylim = c(min(survey.DATA.US.G[,indic.serie], na.rm = T) -0.15, max(survey.DATA.US.G[,indic.serie], na.rm = T)+0.15))
    lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G[,indic.serie])],survey.DATA.US.G[,indic.serie][!is.na(survey.DATA.US.G[,indic.serie])],
          col="blue", lwd="2")
    indic.serie <- colnames(survey.DATA.US.G) == paste(short.names.all.surveys[count.survey],".k3rd",sep="")
    plot(survey.DATA.US.G$date,survey.DATA.US.G[,indic.serie], xlab="Time", ylab="Inflation rate (%)",
         col="red", lwd="2",
         ylim = c(min(survey.DATA.US.G[,indic.serie], na.rm = T) -0.05, max(survey.DATA.US.G[,indic.serie], na.rm = T)+0.05))
    lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G[,indic.serie])],survey.DATA.US.G[,indic.serie][!is.na(survey.DATA.US.G[,indic.serie])],
          col="red", lwd="2")
    
    par(mar=c(0.1, 4, 1, 1))
    plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
    legend("bottom",inset = 0, title="Variable(s):",
           c("Skewness", "Cumulant of order 3"),
           fill=c("blue", "red"), cex=1.1, horiz = TRUE)
    
    ## KURTOSIS AND CUMULANT OF ORDER 4
    indic.serie <- colnames(survey.DATA.US.G) == paste(short.names.all.surveys[count.survey],".kurtosis",sep="")
    plot(survey.DATA.US.G$date,survey.DATA.US.G[,indic.serie], xlab="Time", ylab="Inflation rate (%)",
         col="blue", lwd="2", main=names.all.surveys[count.survey],
         ylim = c(min(survey.DATA.US.G[,indic.serie], na.rm = T) -0.15, max(survey.DATA.US.G[,indic.serie], na.rm = T)+0.15))
    lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G[,indic.serie])],survey.DATA.US.G[,indic.serie][!is.na(survey.DATA.US.G[,indic.serie])],
          col="blue", lwd="2")
    indic.serie <- colnames(survey.DATA.US.G) == paste(short.names.all.surveys[count.survey],".k4th",sep="")
    plot(survey.DATA.US.G$date,survey.DATA.US.G[,indic.serie], xlab="Time", ylab="Inflation rate (%)",
         col="red", lwd="2",
         ylim = c(min(survey.DATA.US.G[,indic.serie], na.rm = T) -0.05, max(survey.DATA.US.G[,indic.serie], na.rm = T)+0.07))
    lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G[,indic.serie])],survey.DATA.US.G[,indic.serie][!is.na(survey.DATA.US.G[,indic.serie])],
          col="red", lwd="2")
    
    par(mar=c(0.1, 4, 1, 1))
    plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
    legend("bottom",inset = 0, title="Variable(s):",
           c("Kurtosis", "Cumulant of order 4"),
           fill=c("blue", "red"), cex=1.1, horiz = TRUE)
    
  } 
  # ============================================  
  # Plot the distribution for the last 6 observations.
  
  multiple <- 4
  survey.bins.distri.without.na <- na.omit(survey.bins.distri) 
  max.m <- dim(survey.bins.distri.without.na)[1]
  n.row <- ceiling(max.m/multiple)
  #seq.m <- seq(max.m-min(5,max.m-1),max.m, by=1)
  
  rest <- n.row*multiple - max.m #max.m %% multiple 
  seq.m <- seq(1,max.m, by=1)
  seq.m.rest <- seq(1,max.m-rest, by=1)
  
  
  
  ## Parameters in order to diplay well the chart
  m <- matrix(c(seq.m,rep(0,rest), rep(max.m+1,multiple)),nrow = (n.row+1), ncol = multiple,byrow = TRUE)
  layout(mat = m, heights = c(rep(0.9/n.row,n.row),0.1))
  
  if(max.m < 6){
    
    if(max.m < 4){
      m <- matrix(c(seq(1,max.m,1),rep(max.m+1,min(length(c(seq(1,max.m,1))),3))),
                  nrow = (ceiling(max.m/3)+1),ncol = min(length(c(seq(1,max.m,1))),3),byrow = TRUE)
      
      layout(mat = m, heights = c(0.42,0.58))
      
    }else if(max.m == 4){m <- matrix(c(1,2,3,4,5,5), ncol=2,nrow = 3,byrow = TRUE)
    layout(mat = m, heights = c(0.84/2,0.84/2,0.16))
    }else{ m <- matrix(c(1,2,3,4,5,5,6,6,6), ncol=3,nrow = 3,byrow = TRUE )
    layout(mat = m, heights = c(0.84/2,0.84/2,0.16))}
    
    
  }
  #layout.show(7)
  par(mar=c(4, 4.5, 1, 1))
  
  for (m in seq.m) {
    
    l.survey <- length(x)
    
    all.data.new <- survey.bins.distri.without.na[m,]
    Date <- all.data.new[,1]
    data1 <- all.data.new[,2:(l.survey+1)]*100
    
    if(indic.US.G.SPF[count.survey]==1){# reverse order
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
    
    if(substr(short.names.all.surveys[count.survey], 1, 3)=="PDS"){# PDS
      mean.class <- c(0.75,seq(1.25, 2.75, 0.5),3.5)
      #break2 <- sort(c(c(0.5,4), seq(1, 3, 0.5), seq(1.45, 2.95, by =0.5)))
      break2 <- sort(c(c(0.5,4), seq(1, 3, 0.5), seq(0.95, 2.95, by =0.5)))
      
    }else if(substr(short.names.all.surveys[count.survey], 1, 3)=="SPF"){# US SPF
      #mean.class <- c(-1,seq(0.25, 7.75, 0.5),9)
      #break2 <- sort(c(c(-2,10), seq(0, 8, 0.5), seq(0.45, 7.95, by =0.5)))
      
      mean.class <- lower.bound.interval + all.intervals/2
      #break2 <- sort(c(upper.bound.interval, lower.bound.interval, lower.bound.interval[-order(lower.bound.interval)[1:2]] - 0.05))
      break2 <- sort(c(upper.bound.interval, lower.bound.interval, lower.bound.interval[-order(lower.bound.interval)[1]] - 0.0))
      break2 <- break2[!duplicated(break2)]
    }
    
    
    plot.fit.survey.distribution.mixture(data1,x1,param1,min.bin=min.bin.fit,max.bin=max.bin.fit,min.sigma=min.sigma, mean.class, break2, xtitle="GDP growth rate %")
    abline(v=all.data.new[,(l.survey+2)], lty=2)
    legend("topleft",legend=paste(t,":",q, sep=""), cex=0.6)  
    legend("topright",legend=paste(horizon[count.survey], sep=""), cex=0.6)

    
  }
  
  par(mar=c(1, 4, 1, 1))
  plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
  legend("bottom",inset = 0, title="Variables:", lty=c('blank','solid'),
         c("Survey (observed)", 'Survey ("Mixture" smoothed)') ,
         col = c("grey","black"), fill=c("grey", "white"), border = c("grey", "white"),
         cex=1, horiz = TRUE)
  
  
}


m <- matrix(c(1,2,3,4,5,5), ncol=2,nrow = 3,byrow = TRUE )
layout(mat = m, heights = c(0.84/2,0.84/2,0.16))
par(mar=c(2, 3, 1, 1))
plot(survey.DATA.US.G$date, survey.DATA.US.G$SPF.US.G.4Q.beta.skew, col="blue", lwd="2", xlab="", ylab="", main="Skewness", ylim=c(-2,1.5), xlim=as.Date(c("1980-01-15","2024-12-15")))
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.4Q.beta.skew)],survey.DATA.US.G$SPF.US.G.4Q.beta.skew[!is.na(survey.DATA.US.G$SPF.US.G.4Q.beta.skew)],
      col="blue", lwd="2", main="Skewness")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.skew)],survey.DATA.US.G$SPF.US.G.5Q.beta.skew[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.skew)], col="red", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.skew)],survey.DATA.US.G$SPF.US.G.5Q.beta.skew[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.skew)],
      col="red", lwd="2", main="Skewness")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.skew)],survey.DATA.US.G$SPF.US.G.6Q.beta.skew[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.skew)], col="black", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.skew)],survey.DATA.US.G$SPF.US.G.6Q.beta.skew[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.skew)],
      col="black", lwd="2", main="Skewness")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.skew)],survey.DATA.US.G$SPF.US.G.7Q.beta.skew[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.skew)], col="purple", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.skew)],survey.DATA.US.G$SPF.US.G.7Q.beta.skew[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.skew)],
      col="purple", lwd="2", main="Skewness")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.skew)],survey.DATA.US.G$SPF.US.G.8Q.beta.skew[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.skew)], col="darkgreen", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.skew)],survey.DATA.US.G$SPF.US.G.8Q.beta.skew[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.skew)],
      col="darkgreen", lwd="2", main="Skewness")
abline(h=0, lty=2)



plot(survey.DATA.US.G$date, survey.DATA.US.G$SPF.US.G.4Q.beta.k3rd, col="blue", lwd="2", xlab="", ylab="", main="3rd Cumulant", ylim=c(-12,8), xlim=as.Date(c("1980-01-15","2024-12-15")))
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.4Q.beta.k3rd)],survey.DATA.US.G$SPF.US.G.4Q.beta.k3rd[!is.na(survey.DATA.US.G$SPF.US.G.4Q.beta.k3rd)],
      col="blue", lwd="2")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.k3rd)],survey.DATA.US.G$SPF.US.G.5Q.beta.k3rd[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.k3rd)], col="red", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.k3rd)],survey.DATA.US.G$SPF.US.G.5Q.beta.k3rd[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.k3rd)],
      col="red", lwd="2")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.k3rd)],survey.DATA.US.G$SPF.US.G.6Q.beta.k3rd[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.k3rd)], col="black", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.k3rd)],survey.DATA.US.G$SPF.US.G.6Q.beta.k3rd[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.k3rd)],
      col="black", lwd="2")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.k3rd)],survey.DATA.US.G$SPF.US.G.7Q.beta.k3rd[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.k3rd)], col="purple", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.k3rd)],survey.DATA.US.G$SPF.US.G.7Q.beta.k3rd[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.k3rd)],
      col="purple", lwd="2")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.k3rd)],survey.DATA.US.G$SPF.US.G.8Q.beta.k3rd[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.k3rd)], col="darkgreen", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.k3rd)],survey.DATA.US.G$SPF.US.G.8Q.beta.k3rd[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.k3rd)],
      col="darkgreen", lwd="2")
abline(h=0, lty=2)



plot(survey.DATA.US.G$date, survey.DATA.US.G$SPF.US.G.4Q.beta.kurtosis, col="blue", lwd="2", xlab="", ylab="", main="Kurtosis", ylim=c(2,14), xlim=as.Date(c("1980-01-15","2024-12-15")))
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.4Q.beta.kurtosis)],survey.DATA.US.G$SPF.US.G.4Q.beta.kurtosis[!is.na(survey.DATA.US.G$SPF.US.G.4Q.beta.kurtosis)],
      col="blue", lwd="2", main="Kurtosis")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.kurtosis)],survey.DATA.US.G$SPF.US.G.5Q.beta.kurtosis[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.kurtosis)], col="red", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.kurtosis)],survey.DATA.US.G$SPF.US.G.5Q.beta.kurtosis[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.kurtosis)],
      col="red", lwd="2")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.kurtosis)],survey.DATA.US.G$SPF.US.G.6Q.beta.kurtosis[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.kurtosis)], col="black", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.kurtosis)],survey.DATA.US.G$SPF.US.G.6Q.beta.kurtosis[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.kurtosis)],
      col="black", lwd="2")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.kurtosis)],survey.DATA.US.G$SPF.US.G.7Q.beta.kurtosis[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.kurtosis)], col="purple", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.kurtosis)],survey.DATA.US.G$SPF.US.G.7Q.beta.kurtosis[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.kurtosis)],
      col="purple", lwd="2")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.kurtosis)],survey.DATA.US.G$SPF.US.G.8Q.beta.kurtosis[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.kurtosis)], col="darkgreen", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.kurtosis)],survey.DATA.US.G$SPF.US.G.8Q.beta.kurtosis[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.kurtosis)],
      col="darkgreen", lwd="2")
abline(h=3, lty=2)




plot(survey.DATA.US.G$date, survey.DATA.US.G$SPF.US.G.4Q.beta.k4th, col="blue", lwd="2", xlab="", ylab="", main="4th Cumulant", ylim=c(-2,50), xlim=as.Date(c("1980-01-15","2024-12-15")))
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.4Q.beta.k4th)],survey.DATA.US.G$SPF.US.G.4Q.beta.k4th[!is.na(survey.DATA.US.G$SPF.US.G.4Q.beta.k4th)],
      col="blue", lwd="2")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.k4th)],survey.DATA.US.G$SPF.US.G.5Q.beta.k4th[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.k4th)], col="red", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.k4th)],survey.DATA.US.G$SPF.US.G.5Q.beta.k4th[!is.na(survey.DATA.US.G$SPF.US.G.5Q.beta.k4th)],
      col="red", lwd="2")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.k4th)],survey.DATA.US.G$SPF.US.G.6Q.beta.k4th[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.k4th)], col="black", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.k4th)],survey.DATA.US.G$SPF.US.G.6Q.beta.k4th[!is.na(survey.DATA.US.G$SPF.US.G.6Q.beta.k4th)],
      col="black", lwd="2")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.k4th)],survey.DATA.US.G$SPF.US.G.7Q.beta.k4th[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.k4th)], col="purple", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.k4th)],survey.DATA.US.G$SPF.US.G.7Q.beta.k4th[!is.na(survey.DATA.US.G$SPF.US.G.7Q.beta.k4th)],
      col="purple", lwd="2")
points(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.k4th)],survey.DATA.US.G$SPF.US.G.8Q.beta.k4th[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.k4th)], col="darkgreen", lwd="2")
lines(survey.DATA.US.G$date[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.k4th)],survey.DATA.US.G$SPF.US.G.8Q.beta.k4th[!is.na(survey.DATA.US.G$SPF.US.G.8Q.beta.k4th)],
      col="darkgreen", lwd="2")
abline(h=0, lty=2)

par(mar=c(1, 4, 1, 1))
plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom",inset = 0, title="Variables:", lty=c('solid','solid','solid','solid'),
       c("4Q", '5Q', '6Q', '7Q', '8Q') ,
       col = c("blue", "red", "black", "purple", "darkgreen"), fill=c("blue", "red", "black", "purple", "darkgreen"),
       cex=1, horiz = TRUE)

dev.off()
