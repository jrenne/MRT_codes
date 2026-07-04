# =======================================================
# Automatic treatment of ECB SPF inflation data
# =======================================================
#1999Q1-2000Q3: T0_0:F3_5 (0.5)
#2000Q4: T0_0:F4_0 (0.5)
#2001Q1:2008Q2: T0_0:F3_5 (0.5)
#2008Q3-2009Q1: T0_0:F4_0 (0.5)
#2009Q2-2009Q4: TN2_0:F4_0 (0.5)
#2010Q1-2020Q1: TN1_0:F4_0 (0.5)
#2020Q2:2020Q4: TN4_0:F4_0 (0.5)
#2021Q1-2022Q2: TN1_0:F4_0 (0.5)
#2022Q3-...: TN1_0:F5_0 (0.5)
all.info.delete.data <- NULL

all.files <- list.files("./data/EA/raw/SPF_ECB")
all.csv.files <- all.files[grep("\\.csv$", all.files, ignore.case = TRUE)]
first.year    <- min(as.numeric(substr(all.csv.files, 1, 4)))
first.quarter <- as.numeric(substr(min(substr(all.csv.files, 1, 6)), 6, 6))
date.with.pb.5y <- c("1999-02-15", "2000-02-15", "2001-02-15")

last.year     <- max(as.numeric(substr(all.csv.files, 1, 4)))
last.quarter  <- as.numeric(substr(max(substr(all.csv.files, 1, 6)), 6, 6))

nb.years <- last.year - first.year + 1
# Compute vector of quarters:

years <- matrix(as.character(first.year:last.year),ncol=4,nrow=nb.years)
quarters <- t(matrix(as.character(1:4),ncol = nb.years,nrow=4))

yq <- array(0,c(nb.years,4,2))
yq[,,1] <- years
yq[,,2] <- quarters
yq <- apply(yq,c(2,1),function(x){paste(x[1],x[2],sep="Q")})
yq.infl <- yq #paste(yq,"_infl",sep="")
if(first.quarter>1){
  yq.infl <- yq.infl[first.quarter:length(yq.infl)]
}
if(last.quarter<4){
  yq.infl <- yq.infl[1:(length(yq.infl)-(4-last.quarter))]
}

count <- 0
for(spf.file in yq.infl){
  count <- count+1
  
  #================================
  # Load data:
  data <- read.table(paste("./data/EA/raw/SPF_ECB/",spf.file,".csv",sep=""), 
                     sep=",", quote="\"",header=TRUE, skip =1, na.strings = c("")) %>%
    slice(1:(which(is.na(TARGET_PERIOD))[1]-1)) %>%
    mutate_at(vars(-c("TARGET_PERIOD")), as.numeric) %>%
    dplyr::select(-matches('^X'))
  
  
  #================================
  # Detect first and last class:
  initial.classes <- colnames(data)[4:dim(data)[2]]
  increment <- diff(round(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",initial.classes))))*2)/2)
  min.classes.lower <- min(round(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",initial.classes))))*2)/2)-max(1,increment)
  min.classes.upper <- min(round(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",initial.classes))))*2)/2)-0.1
  
  # New for 2024Q4 and ...
  if(count>103){
    min.classes.lower <- min(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",initial.classes)))))-max(1,increment)+0.1
    min.classes.upper <- min(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",initial.classes))))*2/2)
  }
  
  
  colnames(data)[4] <- paste0(
    "F",gsub("\\.","_",gsub("-","N",sprintf("%.1f", min.classes.lower))), 
    "T",gsub("\\.","_",gsub("-","N",sprintf("%.1f", min.classes.upper)))
  )
  max.classes.upper <- max(round(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",initial.classes))))*2)/2)+max(1,increment)-0.1
  max.classes.lower <- max(round(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",initial.classes))))*2)/2)
  
  
  # New for 2024Q4 and ...
  if(count>103){
    max.classes.upper <- max(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",initial.classes))))*2/2)+max(1,increment)-0.1
    max.classes.lower <- max(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",initial.classes))))*2/2)
  }
  
  colnames(data)[dim(data)[2]] <- paste0(
    "F",gsub("\\.","_",gsub("-","N",sprintf("%.1f", max.classes.lower))), 
    "T",gsub("\\.","_",gsub("-","N",sprintf("%.1f", max.classes.upper)))
  )
  classes <- colnames(data)[4:dim(data)[2]]
  nb.classes <- length(colnames(data)[4:dim(data)[2]])
  
  #Add column that indicates (=1) if distribution is provided:
  #data$check.not.na <- 1*(apply(data[,4:dim(data)[2]],1,function(x){sum(is.na(x))})<(nb.classes-1))
  data$check.not.na <- 1*(apply(data[,4:dim(data)[2]],1,function(x){sum(is.na(x))})<(nb.classes))
  # If a distri is provided, replace NAs by 0s:
  index.not.na <- which(data$check.not.na==1)
  aux <- data[index.not.na,4:dim(data)[2]]
  aux[is.na(aux)] <- 0
  data[index.not.na,4:dim(data)[2]] <- aux
  
  
  #REMOVE DATA WITH NO DISTRIBUTION.... DELETE THAT FOR THE MOMENT
  #data <- data[index.not.na,1:dim(data)[2]]
  #Add column that indicates (=1) if distribution sum to 100%
  data$check.100 <- apply(data[,4:(dim(data)[2]-1)],1,sum)
  
  # Replace TARGET_PERIOD by common code
  ## ex: Jan by 0.08
  current.quarter <- as.numeric(substr(spf.file,1,4)) + as.numeric(substr(spf.file,6,6))/4 - 1/4
  data$current.quarter <- current.quarter
  data$TARGET_PERIOD_OLD <- data$TARGET_PERIOD
  data$TARGET_PERIOD <- gsub("Jan",".08",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Apr",".33",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Jul",".58",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Oct",".83",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Feb",".16",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("May",".41",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Aug",".66",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Nov",".91",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Mar",".24",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Jun",".49",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Sep",".74",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Dec",".99",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Q1",".24",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Q2",".49",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Q3",".74",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- gsub("Q4",".99",data$TARGET_PERIOD)
  data$TARGET_PERIOD <- as.numeric(data$TARGET_PERIOD)
  
  data$TARGET_PERIOD_OLD <- gsub("Jan","-01-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Apr","-04-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Jul","-07-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Oct","-10-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Feb","-02-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("May","-05-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Aug","-08-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Nov","-11-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Mar","-03-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Jun","-06-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Sep","-09-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Dec","-12-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Q1","-03-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Q2","-06-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Q3","-09-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- gsub("Q4","-12-15",data$TARGET_PERIOD_OLD)
  data$TARGET_PERIOD_OLD <- as.Date(paste(data$TARGET_PERIOD_OLD, "12", "31", sep = "-"))
  # Use ifelse to modify only the years
  data$TARGET_PERIOD_OLD <- ifelse(grepl("^\\d{4}$", as.character(data$TARGET_PERIOD_OLD)), 
                                   paste(data$TARGET_PERIOD_OLD, "12", "31", sep = "-"), 
                                   as.character(data$TARGET_PERIOD_OLD))
  data$TARGET_PERIOD_OLD <- as.Date(data$TARGET_PERIOD_OLD)
  
  
  
  
  # Detect those lines where the horizon is an integer number of years:
  ## integer.horizon = 1 for one year ahead and two year ahead (rolling horizon)
  data$integer.horiz <- (data$TARGET_PERIOD - as.integer(data$TARGET_PERIOD))>0
  data$TARGET_PERIOD[data$TARGET_PERIOD==as.integer(data$TARGET_PERIOD)] <- 
    data$TARGET_PERIOD[data$TARGET_PERIOD==as.integer(data$TARGET_PERIOD)]+1
  data$horizon <- data$TARGET_PERIOD - current.quarter
  # added to get back original target period
  data$TARGET_PERIOD[data$TARGET_PERIOD==as.integer(data$TARGET_PERIOD)] <- 
    data$TARGET_PERIOD[data$TARGET_PERIOD==as.integer(data$TARGET_PERIOD)]-1
  data$horizon[data$horizon<=0] <- data$horizon[data$horizon<=0] + 1
  data$horizon[data$integer.horiz] <- round(data$horizon[data$integer.horiz])
  
  
  # build vector of class' designations:
  inf.classes <- round(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",classes))))*2)/2
  sup.classes <- round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",classes))))*2)/2 
  all.intervals <- sup.classes  - inf.classes
  med.classes <- inf.classes + all.intervals/2
  
  # Eliminate inconsistent forecasts
  data$mean <- (data.matrix(data[,4:(4+nb.classes-1)])%*%matrix(med.classes,ncol=1))/100
  data$mean.inf <- (data.matrix(data[,4:(4+nb.classes-1)])%*%matrix(inf.classes,ncol=1))/100
  data$mean.sup <- (data.matrix(data[,4:(4+nb.classes-1)])%*%matrix(sup.classes,ncol=1))/100
  # data$point.mean.diff <- (abs(data$POINT/data$mean-1)>.25)
  data$point.mean.diff <- (abs(data$POINT/data$mean-1)>.25| data$POINT<data$mean.inf | data$POINT>data$mean.sup)
  
  # Remove non-rational forecasters and forecasters with extreme bins on both sides
  # data <- data %>%
  #    filter(point.mean.diff==FALSE | is.na(point.mean.diff)) 
  
  data$extreme.bins.both.side <- (data[,4]> 15 & data[,(4+nb.classes-1)]>15 )
  #print(sum(data$extreme.bins.both.side, na.rm=T))
  
  # Function to check for zeros between positive bins
  check_zeros_in_between <- function(row) {
    
    if(sum(row > 0, na.rm=T)>1){
      
      first.obs <- which(row > 0)[1]
      last.obs <- tail(which(row > 0),1)
      
      # Remove leading and trailing zeros
      row_no_edges <- row[c(first.obs:last.obs)]
      
      # Check for zeros in between positive bins
      if (any(row_no_edges == 0)) {
        return(TRUE)  # Return TRUE if a zero is found between positive values
      } else {
        return(FALSE)  # Return FALSE otherwise
      }
    } else{
      return(FALSE)
    }
    
  }
  
  # Apply the function row-wise and create a new column 'has_zero_in_between'
  data <- data %>%
    rowwise() %>%
    mutate(has_zero_in_between = check_zeros_in_between(c_across(4:(4+nb.classes-1)))) %>%
    ungroup()  # Remove row-wise grouping
  
  # # Remove data with check.100 below 90
  # data <- data %>%
  #   filter(check.100>90 | is.na(check.100))
  # 
  # # Remove data with extreme bins 
  # data <- data %>%
  #   filter(extreme.bins.both.side==FALSE | is.na(extreme.bins.both.side)) %>%
  #   filter(has_zero_in_between==FALSE | is.na(has_zero_in_between))  # %>%
  #dplyr::select(-has_zero_in_between)
  
  data <- as.data.frame(data)
  #===========================================================================
  #===========================================================================
  # Add column with month with Date format (mid month date):
  #Account for change of SPF publication data
  #Prior to Jan 2015: Feb, May, Aug, Nov and Post Jan 2015: Jan, Apr, Jul, Oct
  data$date <- NaN
  data$year <- as.integer(data$current.quarter)
  aux <- as.factor(sprintf("%4.5f", data$current.quarter))
  indices.before.2015 <- which(data$year<2015)
  indices.after.2015 <- which(data$year>=2015)
  aux1 <- aux
  aux2 <- aux
  aux1 <- as.factor(gsub(".00000","/02/15",aux1))
  aux1 <- as.factor(gsub(".25000","/05/15",aux1))
  aux1 <- as.factor(gsub(".50000","/08/15",aux1))
  aux1 <- as.factor(gsub(".75000","/11/15",aux1))
  aux2 <- as.factor(gsub(".00000","/01/15",aux2))
  aux2 <- as.factor(gsub(".25000","/04/15",aux2))
  aux2 <- as.factor(gsub(".50000","/07/15",aux2))
  aux2 <- as.factor(gsub(".75000","/10/15",aux2))
  aux1 <- as.Date(aux1,"%Y/%m/%d")
  aux2 <- as.Date(aux2,"%Y/%m/%d")
  
  if(sum(data$year<2015)==dim(data)[1]){
    aux <- aux1
  } else{
    aux <- aux2
  }
  # aux <- c(aux1[1:max(indices.before.2015)],
  #          aux2[min(indices.after.2015):max(indices.after.2015)])
  
  data$date <- aux
  
  data$date <- aux
  
  ## PRINT INFO ABOUT FORECASTER TO REMOVE
  # print(paste0("nbr of forecaster removed (extreme bin both side) : ",
  #              sum(data$extreme.bins.both.side, na.rm = T)))
  # print(paste0("nbr of forecaster removed (zero inbetween) : ",
  #              sum(data$has_zero_in_between, na.rm = T)))
  # print(paste0("nbr of forecaster removed (no PE) : ",
  #              sum(is.na(data$POINT))))
  # print(paste0("nbr of forecaster removed (no histogram data) : ",
  #              sum(data$check.not.na==0)))
  # print(paste0("nbr of forecaster removed (no PE and no histogram data) : ",
  #              sum(is.na(data$POINT) | data$check.not.na==0)))
  
  # Compute the number of remove observation per survey
  data$uniqueHorizon <- paste0(data$horizon,data$integer.horiz)
  
  for(h in unique(data$uniqueHorizon)){ #c(1,2,5)
    data.before <- data %>% 
      filter(uniqueHorizon==h)
    
    # if(h < 3){
    #   data.before <- data %>% 
    #     filter(horizon==h,integer.horiz==TRUE) 
    # } else{
    #   data.before <- data %>% 
    #     filter(horizon>4) 
    # }
    
    data.after <- data %>%
      filter(uniqueHorizon==h) %>% 
      #filter(extreme.bins.both.side==FALSE | is.na(extreme.bins.both.side)) %>%
      #filter(has_zero_in_between==FALSE | is.na(has_zero_in_between)) %>%
      filter(check.not.na ==1) %>%
      filter(check.100 >= 90)
    # if(h < 3){
    #   data.after <- data.after %>% 
    #     filter(horizon==h,integer.horiz==TRUE) 
    # } else{
    #   data.after <- data.after %>% 
    #     filter(horizon>4) 
    # }
    
    # create new row
    new.row <- data.frame(
      date = data.before$date[1],
      horizon = data.before$horizon[1],
      inter.horizon = data.before$integer.horiz[1],
      no100 = sum(data.before$check.100 < 90, na.rm=T),
      extreme_bins_both_side = sum(data.before$extreme.bins.both.side, na.rm = TRUE),
      has_zero_in_between = sum(data.before$has_zero_in_between, na.rm = TRUE),
      no_POINT = sum(is.na(data.before$POINT)),
      no_histogramm = sum(data.before$check.not.na == 0),
      no_POINT_no_histogramm = sum(is.na(data.before$POINT) | data.before$check.not.na == 0),
      total_removed = dim(data.before)[1] - dim(data.after)[1],
      total_after_remove = dim(data.after)[1],
      total_before_remove = dim(data.before)[1],
      irrelevant_forecaster_all = sum(data.before$point.mean.diff, na.rm=T),
      irrelevant_forecaster_staying = sum(data.after$point.mean.diff, na.rm=T)
    )
    
    # Append the new row
    all.info.delete.data <- rbind(all.info.delete.data, new.row)
    
  }
  
  data <- data %>%
    dplyr::select(-uniqueHorizon)
  
  # Remove data with check.100 below 90
  data <- data %>%
    filter(check.100>90 | is.na(check.100))
  
  # Remove data with extreme bins 
  data <- data %>%
    filter(extreme.bins.both.side==FALSE | is.na(extreme.bins.both.side)) %>%
    filter(has_zero_in_between==FALSE | is.na(has_zero_in_between))  # %>%
  #dplyr::select(-has_zero_in_between)
  
  # create database    
  eval(parse(text = gsub(" ","",paste("EA.SPF.DISTRI.", spf.file, "<- data %>% 
                                      dplyr::select(date, TARGET_PERIOD, horizon, FCT_SOURCE, POINT, everything()) %>%
                                      rename(Forecaster_ID= FCT_SOURCE)", sep=""))))
  
  ### ROLLING FORECASTS
  eval(parse(text = gsub(" ","",paste("EA.SPF.DISTRI.1y.", spf.file, "<- data %>% 
                                      filter(horizon==1, integer.horiz==TRUE) %>% 
                                      dplyr::select(date, FCT_SOURCE, POINT, TARGET_PERIOD_OLD, matches('[0-9]')) %>%
                                      dplyr::select(-check.100) %>%
                                      mutate_at(vars(-date, -FCT_SOURCE, -POINT, - TARGET_PERIOD_OLD), ~ . / 100) %>%
                                      rename(Forecaster_ID= FCT_SOURCE,
                                        date_target= TARGET_PERIOD_OLD)", sep=""))))
  
  eval(parse(text = gsub(" ","",paste("EA.SPF.DISTRI.2y.", spf.file, "<- data %>% 
                                      filter(horizon==2, integer.horiz==TRUE) %>% 
                                      dplyr::select(date, FCT_SOURCE, POINT, TARGET_PERIOD_OLD, matches('[0-9]')) %>%
                                      dplyr::select(-check.100) %>%
                                      mutate_at(vars(-date, -FCT_SOURCE, -POINT, - TARGET_PERIOD_OLD), ~ . / 100) %>%
                                      rename(Forecaster_ID= FCT_SOURCE,
                                        date_target= TARGET_PERIOD_OLD)", sep=""))))
  
  ### CALENDAR FORECASTS
  eval(parse(text = gsub(" ","",paste("EA.SPF.DISTRI.cy.", spf.file, "<- data %>% 
                                      filter(horizon <= 1, integer.horiz==FALSE) %>% 
                                      dplyr::select(date, FCT_SOURCE, POINT, TARGET_PERIOD_OLD, matches('[0-9]')) %>%
                                      dplyr::select(-check.100) %>%
                                      mutate_at(vars(-date, -FCT_SOURCE, -POINT, - TARGET_PERIOD_OLD), ~ . / 100) %>%
                                      rename(Forecaster_ID= FCT_SOURCE,
                                        date_target= TARGET_PERIOD_OLD)", sep=""))))
  
  eval(parse(text = gsub(" ","",paste("EA.SPF.DISTRI.ny.", spf.file, "<- data %>% 
                                      filter(horizon > 1 & horizon <= 2, integer.horiz==FALSE) %>% 
                                      dplyr::select(date, FCT_SOURCE, POINT, TARGET_PERIOD_OLD, matches('[0-9]')) %>%
                                      dplyr::select(-check.100) %>%
                                      mutate_at(vars(-date, -FCT_SOURCE, -POINT, - TARGET_PERIOD_OLD), ~ . / 100) %>%
                                      rename(Forecaster_ID= FCT_SOURCE,
                                        date_target= TARGET_PERIOD_OLD)", sep=""))))
  
  eval(parse(text = gsub(" ","",paste("EA.SPF.DISTRI.nny.", spf.file, "<- data %>% 
                                      filter(horizon>2 & horizon <= 3, integer.horiz==FALSE) %>% 
                                      dplyr::select(date, FCT_SOURCE, POINT, TARGET_PERIOD_OLD, matches('[0-9]')) %>%
                                      dplyr::select(-check.100) %>%
                                      mutate_at(vars(-date, -FCT_SOURCE, -POINT, - TARGET_PERIOD_OLD), ~ . / 100) %>%
                                      rename(Forecaster_ID= FCT_SOURCE,
                                        date_target= TARGET_PERIOD_OLD)", sep=""))))
  
  
  if(data$date[1] %in% date.with.pb.5y){
    
    eval(parse(text = gsub(" ","",paste("EA.SPF.DISTRI.5y.", spf.file, "<- data %>% 
                                      filter(horizon>4, integer.horiz==FALSE) %>% 
                                      dplyr::select(date, FCT_SOURCE, POINT, TARGET_PERIOD_OLD, matches('[0-9]')) %>%
                                      dplyr::select(-check.100) %>%
                                      mutate_at(vars(-date, -FCT_SOURCE, -POINT, - TARGET_PERIOD_OLD), ~ . / 100) %>%
                                      rename(Forecaster_ID= FCT_SOURCE,
                                        date_target= TARGET_PERIOD_OLD)", sep=""))))
    
    
  } else{
    
    eval(parse(text = gsub(" ","",paste("EA.SPF.DISTRI.5y.", spf.file, "<- data %>% 
                                      filter(horizon>4) %>% 
                                      dplyr::select(date, FCT_SOURCE, POINT, TARGET_PERIOD_OLD, matches('[0-9]')) %>%
                                      dplyr::select(-check.100) %>%
                                      mutate_at(vars(-date, -FCT_SOURCE, -POINT, - TARGET_PERIOD_OLD), ~ . / 100) %>%
                                      rename(Forecaster_ID= FCT_SOURCE,
                                        date_target= TARGET_PERIOD_OLD)", sep=""))))
  }
  
  
  if(count==1){
    SPF.individual.PE.1y <- distinct(data, date, FCT_SOURCE)
  } else{
    SPF.individual.PE.1y <- rbind(SPF.individual.PE.1y,distinct(data, date, FCT_SOURCE))
  }
  
  # SPF.individual.PE.1y <- data %>% 
  #   dplyr::select(date, POINT, horizon, FCT_SOURCE, integer.horiz) %>% 
  #   tidyr::spread(key=horizon, value =POINT)
  
  if(count==1){
    
    SPF.individual.PE <- data %>% 
      #filter((horizon==1 & integer.horiz==TRUE) | (horizon==2 & integer.horiz==TRUE) | (horizon>4)) %>%
      dplyr::select(date, FCT_SOURCE, POINT, horizon, TARGET_PERIOD_OLD, integer.horiz) %>%
      rename(date_target=TARGET_PERIOD_OLD,
             Forecaster_ID=FCT_SOURCE) #%>%
    #mutate(horizon = ifelse(horizon > 2, 5, horizon))
    #mutate(horizon = round(horizon))
    
  } else{
    
    SPF.individual.PE <- rbind(SPF.individual.PE, data %>% 
                                 #filter((horizon==1 & integer.horiz==TRUE) | (horizon==2 & integer.horiz==TRUE) | (horizon>4)) %>%
                                 dplyr::select(date, FCT_SOURCE, POINT, horizon,TARGET_PERIOD_OLD, integer.horiz) %>%
                                 rename(date_target=TARGET_PERIOD_OLD,
                                        Forecaster_ID=FCT_SOURCE) #%>%
                               #mutate(horizon = ifelse(horizon > 2, 5, horizon))
    )
    
  }
  
  
  
  
}

