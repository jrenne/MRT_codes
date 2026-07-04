#===========================================================================
# This program builds a big data.frame object, called DATA,
#   with all data at the monthly frequency
#===========================================================================
# In the end, we have that in:
## FOR SPF EA DATA
## - SPF.data contains the aggregate PE provided by ECB with:
###    - SPF.EA.1y.pe, SPF.EA.2y.pe, SPF.EA.5y.pe,
###    - contains data with calculation (pe, stdev, 3rd.cum) made 
###      based on the medium classes. This is not a very interesting data 
###      as we do the same in a more sophisticated way with Beta distributions.
## - EA.SPF.DISTRI.1y (5y) contains the aggregate bins of the SPF survey
##   of EA. Useful for Beta distribution
## - DATA contains inflation data.

#===========================================================================
# load euro area SPFs
# NB: This requires a previous use of "build.data.frame.R"
#===========================================================================

library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(extraDistr)
library(writexl)
library(eurostat)
library(mFilter)

# Download HICP data directly for the Euro Area (EA zone)
hcpi <- get_eurostat("prc_hicp_midx", filters = list(geo = "EA"), time_format = "date") %>%
  filter(geo == "EA",
         coicop == "CP00",
         unit == "I15") %>%
  mutate(date = time +14) %>%
  dplyr::select(date, values) %>%
  rename(EA.hcpi=values)

ts_hcpi <- ts(hcpi$EA.hcpi, frequency=12, start=c(1996,1))
stl_hcpi <- stl(ts_hcpi, s.window="periodic")
plot(stl_hcpi)
seasonal_stl_hcpi <- stl_hcpi$time.series[,1]
trend_stl_hcpi <- stl_hcpi$time.series[,2]
random_stl_hcpi <- stl_hcpi$time.series[,3]
hcpi$EA.hcpi.deseasonalized <- trend_stl_hcpi + random_stl_hcpi

ts.test <- ts(hcpi$EA.hcpi.deseasonalized[seq(1,length(hcpi$EA.hcpi.deseasonalized),by=3)], frequency=4, start=c(1996,1))
#decomposition trend cycle with HP filter
hcpi.decomposition <- hpfilter(log(trend_stl_hcpi + random_stl_hcpi), drift=FALSE)
plot(hcpi.decomposition$cycle)

hcpi$EA.log.hcpi.deseasonalized.cycle <- as.vector(hcpi.decomposition$cycle*100)
hcpi$EA.log.hcpi.deseasonalized.trend <- as.vector(hcpi.decomposition$trend*100)

plot(hcpi$EA.hcpi, type="l", col="red")
plot(hcpi$EA.hcpi.deseasonalized)  

# Add data to the big data.frame:
DATA <- as.data.frame(hcpi)
T.DATA <- dim(DATA)[1]

# Last month:
LAST.MONTH <- as.numeric(tail(format(DATA$date, "%m"), n=1))
LAST.YEAR <- as.numeric(tail(format(DATA$date, "%Y"), n=1))
#### CHANGE
#LAST.MONTH <- 9 
#LAST.YEAR <- 2022

# Most recent available ECB SPF:
last.year.4.EA.SPF     <- LAST.YEAR
last.quarter.4.EA.SPF  <- trunc((LAST.MONTH-1)/3)+1

# Select sample size:
first.month <- as.Date("15-12-1998","%d-%m-%Y")
#last.month  <- as.Date("15-8-2020","%d-%m-%Y")


DATA$EA.infl <- NA
DATA$EA.infl[13:T.DATA] <- log(DATA$EA.hcpi[13:T.DATA]/DATA$EA.hcpi[1:(T.DATA-12)])*100
plot(DATA$date,DATA$EA.infl,type='l')


DATA$EA.infl.deseasonalized <- NA
DATA$EA.infl.deseasonalized[13:T.DATA] <- log(DATA$EA.hcpi.deseasonalized[13:T.DATA]/DATA$EA.hcpi.deseasonalized[1:(T.DATA-12)])*100
plot(DATA$date,DATA$EA.infl.deseasonalized,type='l')

DATA$EA.infl.monthly <- NA
DATA$EA.infl.monthly[2:T.DATA] <- log(DATA$EA.hcpi[2:T.DATA]/DATA$EA.hcpi[1:(T.DATA-1)])*100
plot(DATA$date,DATA$EA.infl.monthly,type='l')

DATA$EA.infl.monthly.deseasonalized <- NA
DATA$EA.infl.monthly.deseasonalized[2:T.DATA] <- log(DATA$EA.hcpi.deseasonalized[2:T.DATA]/DATA$EA.hcpi.deseasonalized[1:(T.DATA-1)])*100
plot(DATA$date,DATA$EA.infl.monthly.deseasonalized,type='l')

DATA$EA.infl.quarterly <- NA
DATA$EA.infl.quarterly[4:T.DATA] <- log(DATA$EA.hcpi[4:T.DATA]/DATA$EA.hcpi[1:(T.DATA-3)])*100
plot(DATA$date,DATA$EA.infl.quarterly,type='l')

DATA$EA.infl.quarterly.deseasonalized <- NA
DATA$EA.infl.quarterly.deseasonalized[4:T.DATA] <- log(DATA$EA.hcpi.deseasonalized[4:T.DATA]/DATA$EA.hcpi.deseasonalized[1:(T.DATA-3)])*100
plot(DATA$date,DATA$EA.infl.quarterly.deseasonalized,type='l')


DATA <- DATA %>%
  mutate(EA.hcpi.annual = rollsum(EA.hcpi, k=12, fill=NA, align='right'),
         EA.hcpi.annual.deseasonalized = rollsum(EA.hcpi.deseasonalized, k=12, fill=NA, align='right')) %>%
  mutate(EA.infl.annual = log(EA.hcpi.annual/lag(EA.hcpi.annual,12))*100,
         EA.infl.annual.deseasonalized  = log(EA.hcpi.annual.deseasonalized/lag(EA.hcpi.annual,12))*100) # error with lag(EA.hcpi.annual.deseasonalized), but doesn't change the results

DATA.quarter <- DATA %>% 
  mutate(date=floor_date(date, "quarter")) %>% # Create Quarter column
  group_by(date) %>%                                       # Group by Quarter
  summarise(
    EA.infl= mean(EA.infl, na.rm = TRUE),
  ) %>% 
  mutate(
    date = case_when(
      date < as.Date("2016-01-01") ~ date %m+% months(1) + 14,
      date >= as.Date("2016-01-01") ~ date + 14)
  )


#===========================================================================
# load euro area SPFs
# NB: This requires a previous use of "build.data.frame.R"
#===========================================================================

source('./RScripts/EA/build.data.frame_new_new.R', echo=TRUE)

dates_nny_to_keep <- unique((SPF.individual.PE %>% filter(horizon>2 & horizon <=3  & integer.horiz==FALSE))$date)

quarters_nny_to_keep <- paste0(year(dates_nny_to_keep), "Q", quarter(dates_nny_to_keep))
indices_to_keep_nny <- (yq.infl %in% quarters_nny_to_keep)

quarters_5y_to_drop <- c("1999Q2","1999Q3","1999Q4",
                         "2000Q2","2000Q3","2000Q4")
indices_to_keep <- !(yq.infl %in% quarters_5y_to_drop)

# These are the names of the surveys used:
names.all.surveys <- c( ########## Modify
  paste0("EA.SPF.DISTRI.1y.",  yq.infl),
  paste0("EA.SPF.DISTRI.2y.",  yq.infl),
  paste0("EA.SPF.DISTRI.cy.",  yq.infl),
  paste0("EA.SPF.DISTRI.ny.",  yq.infl),
  paste0("EA.SPF.DISTRI.nny.",  yq.infl[indices_to_keep_nny]),
  paste0("EA.SPF.DISTRI.5y.",  yq.infl[indices_to_keep])
)

# Define the horizon to be considered
horizon <- c(str_extract(names.all.surveys, "(?<=EA\\.SPF\\.DISTRI\\.)([^\\.]+)")) # SPF


# Loop to compute the parameters for all the beta distributions
# Loop for each survey (horzion)
for(survey in names.all.surveys){
  
  
  eval(parse(text = gsub(" ","",paste(survey,".avg <-",survey, 
                                      "%>% dplyr::select(-date_target) %>%
                                      summarise_all(~ mean(., na.rm = TRUE))",sep=""))))


}
