#===========================================================================
# This program builds a big data.frame object, called DATA,
#   with all data at the monthly frequency
#===========================================================================

# In the end, we have that in:
## FOR SPF EA DATA
## - SPF.data.G contains the aggregate PE provided by ECB with:
###    - SPF.EA.G.1y.pe, SPF.EA.G.2y.pe, SPF.EA.G.5y.pe,
###    - contains data with calculation (pe, stdev, 3rd.cum) made 
###      based on the medium classes. This is not a very interesting data 
###      as we do the same in a more sophisticated way with Beta distributions.
## - EA.SPF.DISTRI.1y (5y) contains the aggregate bins of the SPF survey
##   of EA. Useful for Beta distribution
## - DATA.G contains gdp data.

#===========================================================================
# Unzip files
#===========================================================================

zipF2<- "./data/EA/raw/SPF_individual_forecasts.zip"
outDir2<-"./data/EA/raw"
unzip(zipF2,exdir=outDir2)

#===========================================================================
# Load euro area GROWTH data
#===========================================================================

gdp_raw <- read.table("./data/EA/raw/namq_10_gdp_page_linear.csv.gz",
                      sep=",", quote="\"", header = T)

gdp <- gdp_raw  %>%
  filter(geo == "EA") %>%
  mutate(date = as.Date(as.yearqtr(paste(substr(TIME_PERIOD, 1, 4), substr(TIME_PERIOD, 6, 7)))) + 14) %>%
  dplyr::select(date, OBS_VALUE) %>%
  rename(EA.gdp=OBS_VALUE)

ts_gdp <- ts(gdp$EA.gdp, frequency=4, start=c(1995,1))
stl_gdp <- stl(ts_gdp, s.window="periodic")
plot(stl_gdp)
seasonal_stl_gdp <- stl_gdp$time.series[,1]
trend_stl_gdp <- stl_gdp$time.series[,2]
random_stl_gdp <- stl_gdp$time.series[,3]
gdp$EA.gdp.deseasonalized <- as.vector(trend_stl_gdp + random_stl_gdp)

#decomposition trend cycle with HP filter
gdp.decomposition <- hpfilter(log(trend_stl_gdp + random_stl_gdp), drift=FALSE)
plot(gdp.decomposition$cycle)

gdp$EA.log.gdp.deseasonalized.cycle <- as.vector(gdp.decomposition$cycle*100)
gdp$EA.log.gdp.deseasonalized.trend <- as.vector(gdp.decomposition$trend*100)


plot(gdp$EA.gdp, type="l", col="red")
plot(gdp$EA.gdp.deseasonalized)  

# Add data to the big data.frame:
DATA.G <- as.data.frame(gdp)
T.DATA.G <- dim(DATA.G)[1]

# Last month:
LAST.MONTH.G <- as.numeric(tail(format(DATA.G$date, "%m"), n=1))
LAST.YEAR.G <- as.numeric(tail(format(DATA.G$date, "%Y"), n=1))
### CHANGE
#LAST.MONTH.G <- 4
#LAST.YEAR.G <- 2022

# Most recent available ECB SPF:
last.year.4.EA.G.SPF     <- LAST.YEAR.G
last.quarter.4.EA.G.SPF  <- trunc(LAST.MONTH.G/3)+1

# Select sample size:
first.month <- as.Date("15-12-1998","%d-%m-%Y")
#last.month  <- as.Date("15-8-2020","%d-%m-%Y")

DATA.G$EA.growth <- NA
DATA.G$EA.growth[5:T.DATA.G] <- log(DATA.G$EA.gdp[5:T.DATA.G]/DATA.G$EA.gdp[1:(T.DATA.G-4)])*100

DATA.G$EA.growth.deseasonalized <- NA
DATA.G$EA.growth.deseasonalized[5:T.DATA.G] <- log(DATA.G$EA.gdp.deseasonalized[5:T.DATA.G]/DATA.G$EA.gdp.deseasonalized[1:(T.DATA.G-4)])*100

DATA.G$EA.growth.quarterly <- NA
DATA.G$EA.growth.quarterly[2:T.DATA.G] <- log(DATA.G$EA.gdp[2:T.DATA.G]/DATA.G$EA.gdp[1:(T.DATA.G-1)])*100
plot(DATA.G$date[!is.na(DATA.G$EA.growth.quarterly)],DATA.G$EA.growth.quarterly[!is.na(DATA.G$EA.growth.quarterly)],type='l')

DATA.G$EA.growth.quarterly.deseasonalized <- NA
DATA.G$EA.growth.quarterly.deseasonalized[2:T.DATA.G] <- log(DATA.G$EA.gdp.deseasonalized[2:T.DATA.G]/DATA.G$EA.gdp.deseasonalized[1:(T.DATA.G-1)])*100
plot(DATA.G$date[!is.na(DATA.G$EA.growth.quarterly.deseasonalized)],DATA.G$EA.growth.quarterly.deseasonalized[!is.na(DATA.G$EA.growth.quarterly.deseasonalized)],type='l')

DATA.G <- DATA.G %>%
  mutate(EA.gdp.annual = rollsum(EA.gdp, k=4, fill=NA, align='right'),
         EA.gdp.annual.deseasonalized = rollsum(EA.gdp.deseasonalized, k=4, fill=NA, align='right')) %>%
  mutate(EA.growth.annual = log(EA.gdp.annual/lag(EA.gdp.annual,4))*100,
         EA.growth.annual.deseasonalized  = log(EA.gdp.annual.deseasonalized/lag(EA.gdp.annual.deseasonalized ,4))*100)

plot(DATA.G$date[!is.na(DATA.G$EA.growth.annual)],DATA.G$EA.growth.annual[!is.na(DATA.G$EA.growth.annual)],type='l')

#===========================================================================
# load euro area SPFs
# NB: This requires a previous use of "build.data.frame.R"
#===========================================================================

source('./RScripts/EA/build.data.frame.GDP.growth_new_new.R', echo=TRUE)


# Look for all dates present in all.data.gdp:
dates.in.SPF <- levels(as.factor(SPF.G.individual.PE.1y$date))

# 1999Q1, 2000Q1 and 2001Q1 contains duplicate

dates_nny_to_keep <- unique((SPF.individual.PE %>% filter(horizon>2 & horizon <=3  & integer.horiz==FALSE))$date)

quarters_nny_to_keep <- paste0(year(dates_nny_to_keep), "Q", quarter(dates_nny_to_keep))
#quarters_nny_to_drop <- yq.infl[!(yq.infl %in% quarters_nny_to_keep)]
indices_to_keep_nny <- (yq.infl %in% quarters_nny_to_keep)

quarters_5y_to_drop <- c("1999Q2","1999Q3","1999Q4",
                         "2000Q2","2000Q3","2000Q4")
indices_to_keep <- !(yq.gdp %in% quarters_5y_to_drop)

# These are the names of the surveys used:
names.all.surveys <- c( ########## Modify
  paste0("EA.G.SPF.DISTRI.1y.",  yq.gdp),
  paste0("EA.G.SPF.DISTRI.2y.",  yq.gdp),
  paste0("EA.G.SPF.DISTRI.cy.",  yq.gdp),
  paste0("EA.G.SPF.DISTRI.ny.",  yq.gdp),
  paste0("EA.G.SPF.DISTRI.nny.",  yq.gdp[indices_to_keep_nny]),
  paste0("EA.G.SPF.DISTRI.5y.",  yq.gdp[indices_to_keep])
)

# Define the horizon to be considered
horizon <- c(str_extract(names.all.surveys, "(?<=EA\\.G\\.SPF\\.DISTRI\\.)([^\\.]+)")) # SPF


# Loop to compute the parameters for all the beta distributions
# Loop for each survey (horzion)
for(survey in names.all.surveys){
  
  
  eval(parse(text = gsub(" ","",paste(survey,".avg <-",survey, 
                                      "%>% dplyr::select(-date_target) %>%
                                      summarise_all(~ mean(., na.rm = TRUE))",sep=""))))
  
  
}

start.date <- DATA.G$date[1] 
end.date <- DATA.G$date[dim(DATA.G)[1]] 
vec.dates.monthly <- data.frame("date"=seq(start.date, end.date, by="months"))

DATA.G <- vec.dates.monthly %>%
  full_join(DATA.G, by ="date")
