#===========================================================================
# This program builds a big data.frame object, called DATA,
#   with all data at the monthly frequency
#===========================================================================
library(readxl)

# In the end, we have that in:
## FOR FED PHILADELPHIA
## - DATA.US contains PE provided by FED Philadelphia with:
###    - SPF.US.1y.pe, SPF.US.10y.pe, SPF.US.5y.pe,
###    - US.cpi
###    - US.SPF.data (NOT USEFUL => DELETE IT FOR THE MODEL).
## - data.SPF.US contains the bins for individual forecasters.
## - US.SPF.DISTRI.1Q (8Q) contains the aggregate bins of the SPF survey
##   of FED Philadelphia. Useful for Beta distribution
## - US.SPF.data contains data with calculation (pe, stdev, 3rd.cum) made 
### based on the medium classes. This is not a very interesting database 
### as we do the same in a more sophisticated way with Beta distribution. 

## FOR PDS SURVEY:
##  - PDS.5y.avg: Bins for Beta distribution horizon 5y 
### A very close database is US.PDS.DISTRI.5. This database contains
### the same survey but the probabilities don't necessary add up to one.
### Better to use PDS.5y5y.avg as sum up to 1.
#### DATA.US$PDS.US.5y.pe contains the PE based on the medium classes
##  - PDS.5y5y.avg: Bins for Beta distribution horizon 5y but 5y in 5y:
### A very close database is US.PDS.DISTRI.5_10. This database contains
### the same survey but the probabilities don't necessary add up to one.
### Better to use PDS.5y5y.avg as sum up to 1.
#### DATA.US$PDS.US.5y5y.pe contains the PE based on the medium classes
##  - US.PDS.DISTRI.5Y: Broken serie... Don't not use it
##  - US.PDS.DISTRI.5_10Y: Broken serie... Don't not use it
##  - US:PDS.data: contains calculated PDS that for the two above series
###   Don't use it.

##### IN THE END, ONE SHOULD USE:
####### DATA.US: - SPF.US.1y.pe, SPF.US.10y.pe, SPF.US.5y.pe, US.cpi
####### US.SPF.DISTRI.1Q (8Q), for Beta, then in survey.DATA (after beta calculation)
####### PDS.5y.avg, for Beta, then in survey.DATA (after beta calculation)
####### PDS.5y5y.avg, for Beta, then in survey.DATA (after beta calculation)

#===========================================================================
# Load US inflation
#===========================================================================

CPI_raw <- read.csv("./data/US/raw/CPIAUCSL.csv")

CPI <- CPI_raw %>%
  rename(date = DATE) %>%
  mutate(date = as.Date(date, format="%Y-%m-%d") + 14)

US.CPI.aux <- data.frame(date=CPI$date,US.cpi=CPI$CPIAUCSL)
US.CPI.q.aux <- data.frame(date=CPI$date,US.cpi=CPI$CPIAUCSL)

US.CPI.q.aux <- US.CPI.q.aux %>% 
  mutate(quarter = as.yearqtr(date)) %>%  # Convert dates to quarters
  group_by(quarter) %>%                   # Group data by quarter
  dplyr::summarize(US.cpi.q = mean(US.cpi)) %>%      # Calculate the mean for each quarter
  mutate(date = as.Date(quarter) + months(1) + 14) %>%   
  dplyr::select(date,US.cpi.q)

inflation.us.q.q <- US.CPI.q.aux %>%
  mutate(US.infl = log(US.cpi.q/lag(US.cpi.q,1))*100) %>%
  dplyr::select(date,US.infl)

##### HP FILTER TO EXTRACT TREND CYCLE COMPONENT

ts_cpi_us <- ts(US.CPI.aux$US.cpi, frequency=12, start=c(1947,1))

#decomposition trend cycle with HP filter
cpi.decomposition.us <- hpfilter(log(ts_cpi_us), drift=FALSE)
plot(cpi.decomposition.us$cycle)

US.CPI.aux$US.log.cpi <- as.vector(log(ts_cpi_us)*100)
US.CPI.aux$US.log.cpi.cycle <- as.vector(cpi.decomposition.us$cycle*100)
US.CPI.aux$US.log.cpi.trend <- as.vector(cpi.decomposition.us$trend*100)

# Add data to the big data.frame:
DATA.US.all <- US.CPI.aux %>%
  mutate(US.infl = log(US.cpi/lag(US.cpi,12))*100,
         US.infl.quarterly = log(US.cpi/lag(US.cpi,3))*100,
         US.infl.monthly = log(US.cpi/lag(US.cpi,1))*100)

DATA.US <- US.CPI.aux %>%
  mutate(US.infl = log(US.cpi/lag(US.cpi,12))*100,
         US.infl.quarterly = log(US.cpi/lag(US.cpi,3))*100,
         US.infl.monthly = log(US.cpi/lag(US.cpi,1))*100) %>%
  filter(date >= "1981-07-15")
  #filter(date >= "1990-01-15")

inflation.us.q.m <- DATA.US %>%
  dplyr::select(date, US.infl.quarterly) %>%
  filter(format(date, "%m") %in% c("01","04","07","10"))

plot(inflation.us.q.m, type="l")
lines(inflation.us.q.q, col="red")


#Plot Inflation US data
plot(DATA.US$date, DATA.US$US.infl.monthly,type='l')
lines(DATA.US$date, DATA.US$US.infl.quarterly, col="red")


#===========================================================================
# load US SPFs
#===========================================================================

SPF_US_aggr.Mean_CPI5YR_Level <-read_excel("./data/US/raw/Mean_CPI5YR_Level.xlsx", 
                                           progress = readxl_progress(), .name_repair = "unique",
                                           na = "#N/A") 

SPF_US_aggr.Inflation <-read_excel("./data/US/raw/Inflation.xlsx", 
                                   progress = readxl_progress(), .name_repair = "unique",
                                   na = "#N/A")

SPF_US_aggr.PE <- merge(SPF_US_aggr.Inflation,SPF_US_aggr.Mean_CPI5YR_Level, by=c("YEAR", "QUARTER")) %>%
  mutate(date = as.Date(paste(YEAR, QUARTER*3-1,"15",sep="-"),format="%Y-%m-%d")) %>%
  dplyr::select(YEAR, QUARTER, date, INFCPI1YR,	INFCPI10YR,	CPI5YR) 

SPF_US_aggr.PE <- rename.vars(SPF_US_aggr.PE,from=c("INFCPI1YR","INFCPI10YR","CPI5YR"),
                              to=c("SPF.US.1y.pe","SPF.US.10y.pe","SPF.US.5y.pe"))
SPF_US_aggr.PE$date <- as.Date(SPF_US_aggr.PE$date,"%m/%d/%Y")
SPF_US_aggr.PE$date <- as.Date(
  paste("15-",
        as.integer(format(SPF_US_aggr.PE$date,"%m"))-1,
        "-",
        as.integer(format(SPF_US_aggr.PE$date,"%y")),sep=""),
  "%d-%m-%y")

# DATA.US: containes inflation data and PE data.
DATA.US <- merge(DATA.US,SPF_US_aggr.PE,by="date",all=TRUE)


###### PE 1-2 year ahead PGDP
# Contains the level only
SPF_US_PE_1_2y_individual <- read_excel("./data/US/raw/Individual_PGDP.xlsx", 
                               progress = readxl_progress(), .name_repair = "unique",
                               na = "#N/A", col_type="numeric") %>%
  mutate(date = as.Date(
    paste("15-",
          QUARTER*3-2,
          "-",
          YEAR,sep=""),
    "%d-%m-%Y")) 

###### VINTAGE PGPD
PGDP_vintages <- read_excel("./data/US/raw/PQvQd.xlsx", 
                            progress = readxl_progress(), .name_repair = "unique",
                            na = "#N/A")

PGDP_vintages_last_4 <- data.frame(date=unique(SPF_US_PE_1_2y_individual$date),
                                   PGDPm5=NA,PGDPm4=NA,PGDPm3=NA,PGDPm2=NA,PGDPm1=NA,PGDP0=NA,PGDP1_bis=NA)
indic.1st.vintage <- grep("P68Q4",colnames(PGDP_vintages))

j=1
for(i in indic.1st.vintage:dim(PGDP_vintages)[2]){
  
  PGDP_vintages_last_4[j,2:8] <- tail(PGDP_vintages[,i][!is.na(PGDP_vintages[,i])],7)
  j = j+1
}
PGDP_vintages_last_4 <- PGDP_vintages_last_4 %>% drop_na()

SPF_US_PE_1_2y_aggregate <- SPF_US_PE_1_2y_individual %>%
  group_by(date) %>%
  dplyr::summarize(YEAR = mean(YEAR, na.rm=TRUE),
                   QUARTER = mean(QUARTER, na.rm=TRUE),
                   PGDP1 = mean(PGDP1, na.rm=TRUE),
                   PGDP2 = mean(PGDP2, na.rm=TRUE),
                   PGDP3 = mean(PGDP3, na.rm=TRUE),
                   PGDP4 = mean(PGDP4, na.rm=TRUE),
                   PGDP5 = mean(PGDP5, na.rm=TRUE),
                   PGDP6 = mean(PGDP6, na.rm=TRUE),
                   PGDPA = mean(PGDPA, na.rm=TRUE),
                   PGDPB = mean(PGDPB, na.rm=TRUE)) %>%
  full_join(PGDP_vintages_last_4, by="date") %>%
  mutate(PGDPA_check = case_when(QUARTER == 1 ~ (PGDP2 + PGDP3 + PGDP4 + PGDP5)/4,
                                 QUARTER == 2 ~ (PGDP1 + PGDP2 + PGDP3 + PGDP4)/4,
                                 QUARTER == 3 ~ (PGDP0 + PGDP1 + PGDP2 + PGDP3)/4,
                                 QUARTER == 4 ~ (PGDPm1 + PGDP0 + PGDP1 + PGDP2)/4),
         PGDPAm1 = case_when(QUARTER == 1 ~ (PGDPm2 + PGDPm1 + PGDP0 + PGDP1_bis)/4,
                             QUARTER == 2 ~ (PGDPm3 + PGDPm2 + PGDPm1 + PGDP0)/4,
                             QUARTER == 3 ~ (PGDPm4 + PGDPm3 + PGDPm2 + PGDPm1)/4,
                             QUARTER == 4 ~ (PGDPm5 + PGDPm4 + PGDPm3 + PGDPm2)/4),
         PE_1y = log(PGDPA/PGDPAm1)*100,
         PE_2y = log(PGDPB/PGDPA)*100,
         diff_A = abs(PGDPA - PGDPA_check)) %>%
  mutate(#PE_1y = ifelse(is.na(diff_A) | diff_A > 1, PE_2y, PE_1y),
         PE_1y = replace(PE_1y, which(diff_A>1), NA),
         PE_2y = replace(PE_2y, which(diff_A>1), NA)) # Problem with 1985-01-15, 1990-01-15, 1986-01-15 => remove them

SPF_US_PE_1_2y_aggregate_reshape <- SPF_US_PE_1_2y_aggregate %>%
  dplyr::select(date,PE_1y, PE_2y) %>%
  filter(date > "1981-04-15") %>%
  tidyr::gather(key=Forecast_Horizon, value =value, PE_1y:PE_2y) %>%
  mutate(Forecast_Horizon = case_when(Forecast_Horizon == "PE_1y" ~ 1,
                                      Forecast_Horizon == "PE_2y" ~ 2)) 

# SPF_US_PE_1_2y_aggregate_bis <- SPF_US_PE_1_2y_individual %>%
#   group_by(date) %>%
#   full_join(PGDP_vintages_last_4, by="date") %>%
#   mutate(PGDPA_check = case_when(QUARTER == 1 ~ (PGDP2 + PGDP3 + PGDP4 + PGDP5)/4,
#                                  QUARTER == 2 ~ (PGDP1 + PGDP2 + PGDP3 + PGDP4)/4,
#                                  QUARTER == 3 ~ (PGDP0 + PGDP1 + PGDP2 + PGDP3)/4,
#                                  QUARTER == 4 ~ (PGDPm1 + PGDP0 + PGDP1 + PGDP2)/4),
#          PGDPAm1 = case_when(QUARTER == 1 ~ (PGDPm2 + PGDPm1 + PGDP0 + PGDP1_bis)/4,
#                              QUARTER == 2 ~ (PGDPm3 + PGDPm2 + PGDPm1 + PGDP0)/4,
#                              QUARTER == 3 ~ (PGDPm4 + PGDPm3 + PGDPm2 + PGDPm1)/4,
#                              QUARTER == 4 ~ (PGDPm5 + PGDPm4 + PGDPm3 + PGDPm2)/4),
#          PE_1y = log(PGDPA/PGDPAm1)*100,
#          PE_2y = log(PGDPB/PGDPA)*100,
#          diff_A = abs(PGDPA - PGDPA_check)) %>%
#   dplyr::summarize(YEAR = mean(YEAR, na.rm=TRUE),
#                    QUARTER = mean(QUARTER, na.rm=TRUE),
#                    PGDP1 = mean(PGDP1, na.rm=TRUE),
#                    PGDP2 = mean(PGDP2, na.rm=TRUE),
#                    PGDP3 = mean(PGDP3, na.rm=TRUE),
#                    PGDP4 = mean(PGDP4, na.rm=TRUE),
#                    PGDP5 = mean(PGDP5, na.rm=TRUE),
#                    PGDP6 = mean(PGDP6, na.rm=TRUE),
#                    PGDPA = mean(PGDPA, na.rm=TRUE),
#                    PGDPB = mean(PGDPB, na.rm=TRUE),
#                    PE_1y = mean(PE_1y, na.rm=TRUE),
#                    PE_2y = mean(PE_2y, na.rm=TRUE),
#                    diff_A = mean(diff_A, na.rm=TRUE)) %>%
#     mutate(PE_1y = replace(PE_1y, which(diff_A>1), NA),
#            PE_2y = replace(PE_2y, which(diff_A>1), NA)) 

### DEFINE THE NAME OF THE SAMPLE CONSIDERED
US.SPF.names <- c("SPF_US_individual_1981_1985", "SPF_US_individual_1985_1992",	
                    "SPF_US_individual_1992_2014",	"SPF_US_individual_post_2014")

# Dates to remove
target_dates <- as.Date(c("1985-01-15","1990-01-15", "1986-01-15"),"%Y-%m-%d")

# Contains Bins but no point estimate
SPF_US_individual <-read_excel("./data/US/raw/Individual_PRPGDP.xlsx", 
                               progress = readxl_progress(), .name_repair = "unique",
                               na = "#N/A", col_type="numeric") %>%
  mutate(date = as.Date(
    paste("15-",
          QUARTER*3-2,
          "-",
          YEAR,sep=""),
    "%d-%m-%Y"))


SPF_US_individual_1981_1985 <- SPF_US_individual %>%
  filter(date > "1981-04-15" & date <= "1985-01-15") %>%
  rename(Year = YEAR,
         Quarter = QUARTER,
         Forecaster_ID = ID) %>%
  dplyr::select(-INDUSTRY) %>%
  tidyr::gather(key=Survey_Bin, value =value, PRPGDP1:PRPGDP20) %>%
  mutate(Forecast_Horizon = case_when(as.numeric(substr(Survey_Bin,7,8)) < 7 ~ 1,
                                      TRUE ~ 2)) %>% # put 1 if 1-6, 2if 7-12
  mutate(Survey_Bin = case_when(Forecast_Horizon == 2 ~ paste0("PRPGDP", as.numeric(substr(Survey_Bin,7,8))-6),
                                TRUE ~ Survey_Bin)) %>%
  spread(key=Survey_Bin, value) %>%
  mutate(F12_0T13_9 = PRPGDP1, #F12_0
         F10_0T11_9 = PRPGDP2, 
         F8_0T9_9 = PRPGDP3, 
         F6_0T7_9 = PRPGDP4, 
         F4_0T5_9 = PRPGDP5, 
         F2_0T3_9 = PRPGDP6) %>% #T4_0
  mutate(Variable_Forecasted = "PRPGDP") %>%
  dplyr::select(-starts_with("PRPGDP"))  %>%
  dplyr::select(Year,	Quarter,	Forecaster_ID,	Variable_Forecasted,	Forecast_Horizon,
                everything()) %>%
  mutate(across(7:(ncol(.)), ~if_else(date %in% target_dates, NA_real_, .)))


SPF_US_individual_1985_1992 <- SPF_US_individual %>%
  filter(date > "1985-01-15" & date <= "1991-10-15") %>%
  rename(Year = YEAR,
         Quarter = QUARTER,
         Forecaster_ID = ID) %>%
  dplyr::select(-INDUSTRY) %>%
  tidyr::gather(key=Survey_Bin, value =value, PRPGDP1:PRPGDP20) %>%
  mutate(Forecast_Horizon = case_when(as.numeric(substr(Survey_Bin,7,8)) < 7 ~ 1,
                                      TRUE ~ 2)) %>% # put 1 if 1-6, 2if 7-12
  mutate(Survey_Bin = case_when(Forecast_Horizon == 2 ~ paste0("PRPGDP", as.numeric(substr(Survey_Bin,7,8))-6),
                                TRUE ~ Survey_Bin)) %>%
  spread(key=Survey_Bin, value) %>%
  mutate(F10_0T11_9 = PRPGDP1, #F10_0
         F8_0T9_9 = PRPGDP2, 
         F6_0T7_9 = PRPGDP3, 
         F4_0T5_9 = PRPGDP4, 
         F2_0T3_9 = PRPGDP5, 
         F0_0T1_9 = PRPGDP6) %>% #T2_0
  mutate(Variable_Forecasted = "PRPGDP") %>%
  dplyr::select(-starts_with("PRPGDP"))  %>%
  dplyr::select(Year,	Quarter,	Forecaster_ID,	Variable_Forecasted,	Forecast_Horizon,
                everything()) %>%
  mutate(across(7:(ncol(.)), ~if_else(date %in% target_dates, NA_real_, .)))


SPF_US_individual_1992_2014 <- SPF_US_individual %>%
  filter(YEAR > 1991 & YEAR < 2014) %>%
  rename(Year = YEAR,
         Quarter = QUARTER,
         Forecaster_ID = ID) %>%
  dplyr::select(-INDUSTRY) %>%
  tidyr::gather(key=Survey_Bin, value =value, PRPGDP1:PRPGDP20) %>%
  mutate(Forecast_Horizon = case_when(as.numeric(substr(Survey_Bin,7,8)) < 11 ~ 1,
                                      TRUE ~ 2)) %>%
  mutate(Survey_Bin = case_when(Forecast_Horizon == 2 ~ paste0("PRPGDP", as.numeric(substr(Survey_Bin,7,8))-10),
                                TRUE ~ Survey_Bin)) %>%
  spread(key=Survey_Bin, value) %>%
  mutate(F8_0T10_0 = PRPGDP1, # F8_0
         F7_0T7_9 = PRPGDP2, 
         F6_0T6_9 = PRPGDP3, 
         F5_0T5_9 = PRPGDP4, 
         F4_0T4_9 = PRPGDP5, 
         F3_0T3_9 = PRPGDP6, 
         F2_0T2_9 = PRPGDP7, 
         F1_0T1_9 = PRPGDP8, 
         F0_0T0_9 = PRPGDP9, 
         FN2_0TN0_1 = PRPGDP10) %>% #T0_0
  mutate(Variable_Forecasted = "PRPGDP") %>%
  dplyr::select(-starts_with("PRPGDP"))  %>%
  dplyr::select(Year,	Quarter,	Forecaster_ID,	Variable_Forecasted,	Forecast_Horizon,
                everything()) %>%
  mutate(across(7:(ncol(.)), ~if_else(date %in% target_dates, NA_real_, .)))


SPF_US_individual_post_2014 <- SPF_US_individual %>%
  filter(YEAR > 2013) %>%
  rename(Year = YEAR,
         Quarter = QUARTER,
         Forecaster_ID = ID) %>%
  dplyr::select(-INDUSTRY) %>%
  tidyr::gather(key=Survey_Bin, value =value, PRPGDP1:PRPGDP20) %>%
  mutate(Forecast_Horizon = case_when(as.numeric(substr(Survey_Bin,7,8)) < 11 ~ 1,
                                      TRUE ~ 2)) %>%
  mutate(Survey_Bin = case_when(Forecast_Horizon == 2 ~ paste0("PRPGDP", as.numeric(substr(Survey_Bin,7,8))-10),
                                TRUE ~ Survey_Bin)) %>%
  spread(key=Survey_Bin, value) %>%
  mutate(F4_0_T4_9 = PRPGDP1,	F3_5T3_9 = PRPGDP2, #F4_0
         F3_0T3_4 = PRPGDP3,	F2_5T2_9 = PRPGDP4,
         F2_0T2_4 = PRPGDP5,	F1_5T1_9 = PRPGDP6,
         F1_0T1_4 = PRPGDP7,	F0_5T0_9 = PRPGDP8,
         F0_0T0_4  = PRPGDP9,	FN1_0TN0_1 = PRPGDP10) %>% #T0_0
  mutate(Variable_Forecasted = "PRPGDP") %>%
  dplyr::select(-starts_with("PRPGDP")) %>%
  dplyr::select(Year,	Quarter,	Forecaster_ID,	Variable_Forecasted,	Forecast_Horizon,
                everything()) %>%
  mutate(across(7:(ncol(.)), ~if_else(date %in% target_dates, NA_real_, .)))

US.SPF.dates <- SPF_US_individual %>% dplyr::select(date) %>% filter(date > "1981-04-15", !duplicated(date))

US.SPF.data <- data.frame(date=US.SPF.dates,
                          SPF.US.pe.1Q=NaN,SPF.US.pe.2Q=NaN,SPF.US.pe.3Q=NaN,SPF.US.pe.4Q=NaN,
                          SPF.US.pe.5Q=NaN,SPF.US.pe.6Q=NaN,SPF.US.pe.7Q=NaN,SPF.US.pe.8Q=NaN,
                          SPF.US.stdev.1Q=NaN,SPF.US.stdev.2Q=NaN,SPF.US.stdev.3Q=NaN,SPF.US.stdev.4Q=NaN,
                          SPF.US.stdev.5Q=NaN,SPF.US.stdev.6Q=NaN,SPF.US.stdev.7Q=NaN,SPF.US.stdev.8Q=NaN,
                          SPF.US.disagreement.1Q=NaN,SPF.US.disagreement.2Q=NaN,SPF.US.disagreement.3Q=NaN,SPF.US.disagreement.4Q=NaN,
                          SPF.US.disagreement.5Q=NaN,SPF.US.disagreement.6Q=NaN,SPF.US.disagreement.7Q=NaN,SPF.US.disagreement.8Q=NaN,
                          SPF.US.thirdcumulant.1Q=NaN,SPF.US.thirdcumulant.2Q=NaN,SPF.US.thirdcumulant.3Q=NaN,SPF.US.thirdcumulant.4Q=NaN,
                          SPF.US.thirdcumulant.5Q=NaN,SPF.US.thirdcumulant.6Q=NaN,SPF.US.thirdcumulant.7Q=NaN,SPF.US.thirdcumulant.8Q=NaN,
                          SPF.US.fourthcumulant.1Q=NaN,SPF.US.fourthcumulant.2Q=NaN,SPF.US.fourthcumulant.3Q=NaN,SPF.US.fourthcumulant.4Q=NaN,
                          SPF.US.fourthcumulant.5Q=NaN,SPF.US.fourthcumulant.6Q=NaN,SPF.US.fourthcumulant.7Q=NaN,SPF.US.fourthcumulant.8Q=NaN
)

count <- 0
for(spf.sample in US.SPF.names){
  
  # Create data.SPF.US that contains the considered sample
  eval(parse(text = gsub(" ","",paste("data.SPF.US <- ",
                                      spf.sample, sep="")))) 
  
  data.SPF.US$check.na <- !is.na(data.SPF.US[,7])
  data.SPF.US <- data.SPF.US %>%
    dplyr::select(Year,	Quarter,	Forecaster_ID,	Variable_Forecasted,	Forecast_Horizon, date, check.na,
                  everything())
  
  lower.bound.interval.initial <- round(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",colnames(data.SPF.US)[8:dim(data.SPF.US)[2]]))))*2)/2
  upper.bound.interval.initial <- round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(data.SPF.US)[8:dim(data.SPF.US)[2]]))))*2)/2 
  all.intervals.initial <- upper.bound.interval.initial - lower.bound.interval.initial
  
  med.classes.US.initial <- lower.bound.interval.initial + all.intervals.initial/2

  US.SPF.dates.indic <- data.SPF.US %>% dplyr::select(date) %>% filter(!duplicated(date))
  
  nbr.bins.pre <-length(colnames(data.SPF.US))-7
  
  # The distributions will be saved in dataframes US.SPF.DISTRI.XQ
  US.SPF.DISTRI.1Q <- data.frame(date=US.SPF.dates.indic$date)
  US.SPF.DISTRI.1Q <- cbind(US.SPF.DISTRI.1Q,matrix(NaN,dim(US.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.SPF.DISTRI.1Q) <- colnames(data.SPF.US)[c(6,8:dim(data.SPF.US)[2])]
  US.SPF.DISTRI.2Q <- data.frame(date=US.SPF.dates.indic$date)
  US.SPF.DISTRI.2Q <- cbind(US.SPF.DISTRI.2Q,matrix(NaN,dim(US.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.SPF.DISTRI.2Q) <- colnames(data.SPF.US)[c(6,8:dim(data.SPF.US)[2])]
  US.SPF.DISTRI.3Q <- data.frame(date=US.SPF.dates.indic$date)
  US.SPF.DISTRI.3Q <- cbind(US.SPF.DISTRI.3Q,matrix(NaN,dim(US.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.SPF.DISTRI.3Q) <- colnames(data.SPF.US)[c(6,8:dim(data.SPF.US)[2])]
  US.SPF.DISTRI.4Q <- data.frame(date=US.SPF.dates.indic$date)
  US.SPF.DISTRI.4Q <- cbind(US.SPF.DISTRI.4Q,matrix(NaN,dim(US.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.SPF.DISTRI.4Q) <- colnames(data.SPF.US)[c(6,8:dim(data.SPF.US)[2])]
  US.SPF.DISTRI.5Q <- data.frame(date=US.SPF.dates.indic$date)
  US.SPF.DISTRI.5Q <- cbind(US.SPF.DISTRI.5Q,matrix(NaN,dim(US.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.SPF.DISTRI.5Q) <- colnames(data.SPF.US)[c(6,8:dim(data.SPF.US)[2])]
  US.SPF.DISTRI.6Q <- data.frame(date=US.SPF.dates.indic$date)
  US.SPF.DISTRI.6Q <- cbind(US.SPF.DISTRI.6Q,matrix(NaN,dim(US.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.SPF.DISTRI.6Q) <- colnames(data.SPF.US)[c(6,8:dim(data.SPF.US)[2])]
  US.SPF.DISTRI.7Q <- data.frame(date=US.SPF.dates.indic$date)
  US.SPF.DISTRI.7Q <- cbind(US.SPF.DISTRI.7Q,matrix(NaN,dim(US.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.SPF.DISTRI.7Q) <- colnames(data.SPF.US)[c(6,8:dim(data.SPF.US)[2])]
  US.SPF.DISTRI.8Q <- data.frame(date=US.SPF.dates.indic$date)
  US.SPF.DISTRI.8Q <- cbind(US.SPF.DISTRI.8Q,matrix(NaN,dim(US.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.SPF.DISTRI.8Q) <- colnames(data.SPF.US)[c(6,8:dim(data.SPF.US)[2])]
  
  
  for(y in US.SPF.dates.indic$date){
    indic.month <- which(US.SPF.dates.indic$date==c(y))
    count <- count + 1
    print(as.Date(y))
    for(horizon in 1:2){
      data.aux <- subset(data.SPF.US,(date==c(y))&(Forecast_Horizon==horizon)&check.na)
      avg.distri <- apply(data.aux[,8:dim(data.aux)[2]],2,mean,na.rm=T)/100
      
      # Modify first or last bin if one of the two if > 5% (take the larger of the two)
      if(max(avg.distri[1], avg.distri[length(avg.distri)])>0.1 & !is.na(max(avg.distri[1], avg.distri[length(avg.distri)]))){
        
        # Extract PE data
        data.PE <- SPF_US_PE_1_2y_aggregate_reshape %>% filter(date==c(y)& (Forecast_Horizon==horizon))
        # Identify if we need to change the first or the last bins
        max.prob <- max(avg.distri[1], avg.distri[length(avg.distri)])
        indic.pos.max.prob <- match(max.prob, avg.distri)
        
        # Compute the interval that allows to fit perfectly the PE
        new.min.max.interval <- (data.PE$value - sum(((avg.distri)*(med.classes.US.initial))[-indic.pos.max.prob]))/max.prob

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
          all.intervals <- upper.bound.interval.initial - lower.bound.interval
          med.classes.US <- lower.bound.interval + all.intervals/2
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
          all.intervals <- upper.bound.interval - lower.bound.interval.initial
          med.classes.US <- lower.bound.interval.initial + all.intervals/2
        }
        
      } else{ #Condition to consider initial value
        med.classes.US <- med.classes.US.initial
      }
      
      mean.per.forecaster <- apply((data.aux[,8:dim(data.aux)[2]]/100) * (matrix(1,dim(data.aux)[1],1) %*% med.classes.US),1,sum)
      disagreement.aux <- sd(mean.per.forecaster)
      mean.aux  <- sum((avg.distri)*(med.classes.US))
      stdev.aux <- sqrt(sum((avg.distri)*(med.classes.US)^2) - mean.aux^2)
      thirdcumulant.aux  <- sum((avg.distri)*(med.classes.US - mean.aux)^3)
      fourthcumulant.aux <- sum((avg.distri)*(med.classes.US - mean.aux)^4) - 3*(stdev.aux^2)^2
      # Simple average of forecasts??... potential pbm: not the same forecasters for horizons 1 and 2
      if(horizon==1){
        if(data.aux$Quarter[1]==1 & !is.na(data.aux$Quarter[1])){
          US.SPF.data$SPF.US.pe.4Q[count]    <- mean.aux
          US.SPF.data$SPF.US.stdev.4Q[count] <- stdev.aux
          US.SPF.data$SPF.US.thirdcumulant.4Q[count] <- thirdcumulant.aux
          US.SPF.data$SPF.US.fourthcumulant.4Q[count] <- fourthcumulant.aux
          US.SPF.data$SPF.US.disagreement.4Q[count] <- disagreement.aux
          US.SPF.DISTRI.4Q[indic.month,2:dim(US.SPF.DISTRI.1Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==2 & !is.na(data.aux$Quarter[1])){
          US.SPF.data$SPF.US.pe.3Q[count]    <- mean.aux
          US.SPF.data$SPF.US.stdev.3Q[count] <- stdev.aux
          US.SPF.data$SPF.US.thirdcumulant.3Q[count] <- thirdcumulant.aux
          US.SPF.data$SPF.US.fourthcumulant.3Q[count] <- fourthcumulant.aux
          US.SPF.data$SPF.US.disagreement.3Q[count] <- disagreement.aux
          US.SPF.DISTRI.3Q[indic.month,2:dim(US.SPF.DISTRI.1Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==3 & !is.na(data.aux$Quarter[1])){
          US.SPF.data$SPF.US.pe.2Q[count]    <- mean.aux
          US.SPF.data$SPF.US.stdev.2Q[count] <- stdev.aux
          US.SPF.data$SPF.US.thirdcumulant.2Q[count] <- thirdcumulant.aux
          US.SPF.data$SPF.US.fourthcumulant.2Q[count] <- fourthcumulant.aux
          US.SPF.data$SPF.US.disagreement.2Q[count] <- disagreement.aux
          US.SPF.DISTRI.2Q[indic.month,2:dim(US.SPF.DISTRI.1Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==4 & !is.na(data.aux$Quarter[1])){
          US.SPF.data$SPF.US.pe.1Q[count]    <- mean.aux
          US.SPF.data$SPF.US.stdev.1Q[count] <- stdev.aux
          US.SPF.data$SPF.US.thirdcumulant.1Q[count] <- thirdcumulant.aux
          US.SPF.data$SPF.US.fourthcumulant.1Q[count] <- fourthcumulant.aux
          US.SPF.data$SPF.US.disagreement.1Q[count] <- disagreement.aux
          US.SPF.DISTRI.1Q[indic.month,2:dim(US.SPF.DISTRI.1Q)[2]] <- avg.distri
        }
      }else{#horizon == 2
        if(data.aux$Quarter[1]==1 & !is.na(data.aux$Quarter[1])){
          US.SPF.data$SPF.US.pe.8Q[count]    <- mean.aux
          US.SPF.data$SPF.US.stdev.8Q[count] <- stdev.aux
          US.SPF.data$SPF.US.disagreement.8Q[count] <- disagreement.aux
          US.SPF.data$SPF.US.thirdcumulant.8Q[count] <- thirdcumulant.aux
          US.SPF.data$SPF.US.fourthcumulant.8Q[count] <- fourthcumulant.aux
          US.SPF.DISTRI.8Q[indic.month,2:dim(US.SPF.DISTRI.1Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==2 & !is.na(data.aux$Quarter[1])){
          US.SPF.data$SPF.US.pe.7Q[count]    <- mean.aux
          US.SPF.data$SPF.US.stdev.7Q[count] <- stdev.aux
          US.SPF.data$SPF.US.thirdcumulant.7Q[count] <- thirdcumulant.aux
          US.SPF.data$SPF.US.fourthcumulant.7Q[count] <- fourthcumulant.aux
          US.SPF.data$SPF.US.disagreement.7Q[count] <- disagreement.aux
          US.SPF.DISTRI.7Q[indic.month,2:dim(US.SPF.DISTRI.1Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==3 & !is.na(data.aux$Quarter[1])){
          US.SPF.data$SPF.US.pe.6Q[count]    <- mean.aux
          US.SPF.data$SPF.US.stdev.6Q[count] <- stdev.aux
          US.SPF.data$SPF.US.thirdcumulant.6Q[count] <- thirdcumulant.aux
          US.SPF.data$SPF.US.fourthcumulant.6Q[count] <- fourthcumulant.aux
          US.SPF.data$SPF.US.disagreement.6Q[count] <- disagreement.aux
          US.SPF.DISTRI.6Q[indic.month,2:dim(US.SPF.DISTRI.1Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==4 & !is.na(data.aux$Quarter[1])){
          US.SPF.data$SPF.US.pe.5Q[count]    <- mean.aux
          US.SPF.data$SPF.US.stdev.5Q[count] <- stdev.aux
          US.SPF.data$SPF.US.thirdcumulant.5Q[count] <- thirdcumulant.aux
          US.SPF.data$SPF.US.fourthcumulant.5Q[count] <- fourthcumulant.aux
          US.SPF.data$SPF.US.disagreement.5Q[count] <- disagreement.aux
          US.SPF.DISTRI.5Q[indic.month,2:dim(US.SPF.DISTRI.1Q)[2]] <- avg.distri
        }
      }
    }
    print(sum(avg.distri))
  }
  
  eval(parse(text = gsub(" ","",paste("US.SPF.DISTRI.1Q.", str_sub(spf.sample,-9,-1), "<- US.SPF.DISTRI.1Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.SPF.DISTRI.2Q.", str_sub(spf.sample,-9,-1), "<- US.SPF.DISTRI.2Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.SPF.DISTRI.3Q.", str_sub(spf.sample,-9,-1), "<- US.SPF.DISTRI.3Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.SPF.DISTRI.4Q.", str_sub(spf.sample,-9,-1), "<- US.SPF.DISTRI.4Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.SPF.DISTRI.5Q.", str_sub(spf.sample,-9,-1), "<- US.SPF.DISTRI.5Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.SPF.DISTRI.6Q.", str_sub(spf.sample,-9,-1), "<- US.SPF.DISTRI.6Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.SPF.DISTRI.7Q.", str_sub(spf.sample,-9,-1), "<- US.SPF.DISTRI.7Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.SPF.DISTRI.8Q.", str_sub(spf.sample,-9,-1), "<- US.SPF.DISTRI.8Q", sep="")))) 
  
}


# Add data to the big data.frame:
DATA.US <- merge(DATA.US,US.SPF.data,by="date",all=TRUE)




plot(US.SPF.data$date,US.SPF.data$SPF.US.stdev.1Q,ylim=c(0,2))
points(US.SPF.data$date,US.SPF.data$SPF.US.stdev.3Q,col='blue')
points(US.SPF.data$date,US.SPF.data$SPF.US.stdev.8Q,col='red')

#===========================================================================
# load US PDSs
#===========================================================================

# =====================
# 5y in 5y: Result in PDS.5y5y.avg

PDS_CPI <- read.csv("./data/US/PDS/PDS_CPI_avg_5_10.csv")

PDS_CPI$date <- as.Date(PDS_CPI$date,"%m/%d/%y")
PDS_CPI$date <- as.Date(paste(
  "15-",
  ifelse(as.integer(format(PDS_CPI$date,"%d"))>=18,format(PDS_CPI$date,"%m"),
         ifelse(format(PDS_CPI$date,"%m")=="01",12,
                as.character(as.integer(format(PDS_CPI$date,"%m"))-1))),
  "-",
  ifelse(as.integer(format(PDS_CPI$date,"%d"))>=18,format(PDS_CPI$date,"%y"),
         ifelse(format(PDS_CPI$date,"%m")=="01",as.character(as.integer(format(PDS_CPI$date,"%y"))-1),
                format(PDS_CPI$date,"%y"))),
  sep=""),
  "%d-%m-%y")


PDS.5y5y.avg <- PDS_CPI[,1:7]
PDS.5y5y.avg[,2:7] <- PDS.5y5y.avg[,2:7]/100
sum.PDS.5y5y.avg <- apply(PDS.5y5y.avg[,2:7],1,sum) #Check if sum to one
PDS.5y5y.avg[,2:7] <- PDS.5y5y.avg[,2:7] / (matrix(apply(PDS.5y5y.avg[,2:7],1,sum),ncol=1) %*% matrix(1,1,dim(PDS.5y5y.avg)[2])) # Divide by sum of probability to obtain 1

US.PDS.DISTRI.5_10 <- data.frame(date=DATA.US$date)
US.PDS.DISTRI.5_10 <- merge(US.PDS.DISTRI.5_10,PDS_CPI,by="date",all=TRUE)

vector.of.PDS.values <- c(.5,1.25,1.75,2.25,2.75,3.5)
DATA.US$PDS.US.5y5y.pe <- as.matrix(US.PDS.DISTRI.5_10[,2:7]/100) %*% vector.of.PDS.values
DATA.US$PDS.US.5y5y.stdv <- sqrt(
  as.matrix(US.PDS.DISTRI.5_10[,2:7]/100) %*% (vector.of.PDS.values^2)
  - DATA.US$PDS.US.5y5y.pe^2)
points(DATA.US$date,DATA.US$PDS.US.5y5y.pe,pch=4)


# =====================
# 5y: Result in PDS.5y.avg

PDS_CPI <- read.csv("./data/US/PDS/PDS_CPI_avg_5.csv")

PDS_CPI$date <- as.Date(PDS_CPI$date,"%m/%d/%y")
PDS_CPI$date <- as.Date(paste(
  "15-",
  ifelse(as.integer(format(PDS_CPI$date,"%d"))>=18,format(PDS_CPI$date,"%m"),
         ifelse(format(PDS_CPI$date,"%m")=="01",12,
                as.character(as.integer(format(PDS_CPI$date,"%m"))-1))),
  "-",
  ifelse(as.integer(format(PDS_CPI$date,"%d"))>=18,format(PDS_CPI$date,"%y"),
         ifelse(format(PDS_CPI$date,"%m")=="01",as.character(as.integer(format(PDS_CPI$date,"%y"))-1),
                format(PDS_CPI$date,"%y"))),
  sep=""),
  "%d-%m-%y")

PDS.5y.avg <- PDS_CPI[,1:7]
PDS.5y.avg[,2:7] <- PDS.5y.avg[,2:7]/100
sum.PDS.5y.avg <- apply(PDS.5y.avg[,2:7],1,sum)
PDS.5y.avg[,2:7] <- PDS.5y.avg[,2:7] / (matrix(apply(PDS.5y.avg[,2:7],1,sum),ncol=1) %*% matrix(1,1,dim(PDS.5y.avg)[2]))

US.PDS.DISTRI.5 <- data.frame(date=DATA.US$date)
US.PDS.DISTRI.5 <- merge(US.PDS.DISTRI.5,PDS_CPI,by="date",all=TRUE)

vector.of.PDS.values <- c(.5,1.25,1.75,2.25,2.75,3.5)
DATA.US$PDS.US.5y.pe <- as.matrix(US.PDS.DISTRI.5[,2:7]/100) %*% vector.of.PDS.values
DATA.US$PDS.US.5y.stdv <- sqrt(
  as.matrix(US.PDS.DISTRI.5[,2:7]/100) %*% (vector.of.PDS.values^2)
  - DATA.US$PDS.US.5y.pe^2)
#plot(DATA$date,DATA$PDS.US.5y.stdv,pch=4)



# Extract PDS from extended dataset: Incomplete Data
PDSfrom2007 <- read.csv("./data/US/PDS/PDSfrom2007.csv")
PDSfrom2007$DATE <- as.Date(PDSfrom2007$DATE,"%d/%m/%y")

US.PDS.DISTRI.5Y <- data.frame(date=DATA.US$date)
US.PDS.DISTRI.5_10Y <- data.frame(date=DATA.US$date)
US.PDS.DISTRI.5Y    <- cbind(US.PDS.DISTRI.5Y   ,matrix(NaN,dim(US.PDS.DISTRI.5Y)[1],6))
US.PDS.DISTRI.5_10Y <- cbind(US.PDS.DISTRI.5_10Y,matrix(NaN,dim(US.PDS.DISTRI.5Y)[1],6))

PDSfrom2007$date <- as.Date(paste(
  "15-",
  ifelse(as.integer(format(PDSfrom2007$DATE,"%d"))<20,format(PDSfrom2007$DATE,"%m"),
         ifelse(format(PDSfrom2007$DATE,"%m")=="12",1,
                as.character(as.integer(format(PDSfrom2007$DATE,"%m"))+1))),
  "-",
  ifelse(as.integer(format(PDSfrom2007$DATE,"%d"))<20,format(PDSfrom2007$DATE,"%y"),
         ifelse(format(PDSfrom2007$DATE,"%m")=="12",as.character(as.integer(format(PDSfrom2007$DATE,"%y"))+1),
                format(PDSfrom2007$DATE,"%y"))),
  sep=""),
  "%d-%m-%y")

PDSfrom2007$Horizon <- ifelse(PDSfrom2007$Horizon>8,10,5)

PDSfrom2007$sum.proba <- apply(PDSfrom2007[,4:9],1,function(x){sum(x,na.rm=TRUE)})
PDSfrom2007$check.na <- (PDSfrom2007$sum.proba==1)


US.PDS.dates <- levels(as.factor(PDSfrom2007$date))


US.PDS.data <- data.frame(date=as.Date(US.PDS.dates),
                          PDS.pe.5=NaN,PDS.stdev.5=NaN,PDS.disagreement.5=NaN,
                          PDS.pe.5.10=NaN,PDS.stdev.5.10=NaN,PDS.disagreement.5.10=NaN,
                          PDS.L100=NaN,PDS.100.150=NaN,PDS.150.200=NaN,PDS.200.250=NaN,
                          PDS.250.300=NaN,PDS.H300=NaN)

count <- 0
for(y in US.PDS.dates){
  indic.month <- which(DATA.US$date==y)
  count <- count + 1
  for(horizon in c(5,10)){
    data.aux <- subset(PDSfrom2007,(date==y)&(Horizon==horizon)&check.na)
    avg.distri <- apply(data.aux[,4:9],2,function(x){mean(x,na.rm=TRUE)})
    if(dim(data.aux)[1]>0){
      avg.distri[is.na(avg.distri)] <- 0
      mean.per.forecaster <- apply((data.aux[,4:9]) * (matrix(1,dim(data.aux)[1],1) %*% vector.of.PDS.values),1,function(x){sum(x,na.rm=TRUE)})
      disagreement.aux <- sd(mean.per.forecaster)
      mean.aux <- sum((avg.distri)*(vector.of.PDS.values))
      stdev.aux <- sqrt(sum((avg.distri)*(vector.of.PDS.values)^2) - mean.aux^2)
      thirdcumulant.aux <- sum((avg.distri)*(vector.of.PDS.values - mean.aux)^3)
      # Simple average of forecasts??... potential pbm: not the same forecasters for horizons 1 and 2
      if(horizon==5){
        US.PDS.data$PDS.pe.5[count]    <- mean.aux
        US.PDS.data$PDS.stdev.5[count] <- stdev.aux
        US.PDS.data$PDS.disagreement.5[count] <- disagreement.aux
        US.PDS.DISTRI.5Y[indic.month,2:dim(US.PDS.DISTRI.5Y)[2]] <- avg.distri
      }else{#horizon == 10
        US.PDS.data$PDS.pe.5.10[count]    <- mean.aux
        US.PDS.data$PDS.stdev.5.10[count] <- stdev.aux
        US.PDS.data$PDS.disagreement.5.10[count] <- disagreement.aux
        US.PDS.DISTRI.5_10Y[indic.month,2:dim(US.PDS.DISTRI.5Y)[2]] <- avg.distri
        US.PDS.data$PDS.L100[count]    = avg.distri[1]
        US.PDS.data$PDS.100.150[count] = avg.distri[2]
        US.PDS.data$PDS.150.200[count] = avg.distri[3]
        US.PDS.data$PDS.200.250[count] = avg.distri[4]
        US.PDS.data$PDS.250.300[count] = avg.distri[5]
        US.PDS.data$PDS.H300[count]    = avg.distri[6]
      }
    }
  }
}

#Add colnames
colnames(US.PDS.DISTRI.5Y) <- c("date", colnames(PDSfrom2007)[4:9])
colnames(US.PDS.DISTRI.5_10Y) <- c("date", colnames(PDSfrom2007)[4:9])


# Add data to the big data.frame:
DATA.US <- merge(DATA.US,US.PDS.data,by="date",all=TRUE)


plot(SPF_US_PE_1_2y_aggregate$date,SPF_US_PE_1_2y_aggregate$PE_1y, type="l", lwd=2)
points(US.SPF.data$date[!is.na(US.SPF.data$SPF.US.pe.1Q)],US.SPF.data$SPF.US.pe.1Q[!is.na(US.SPF.data$SPF.US.pe.1Q)], col="blue")
points(US.SPF.data$date[!is.na(US.SPF.data$SPF.US.pe.2Q)],US.SPF.data$SPF.US.pe.2Q[!is.na(US.SPF.data$SPF.US.pe.2Q)], col="red")
points(US.SPF.data$date[!is.na(US.SPF.data$SPF.US.pe.3Q)],US.SPF.data$SPF.US.pe.3Q[!is.na(US.SPF.data$SPF.US.pe.3Q)], col="green")
points(US.SPF.data$date[!is.na(US.SPF.data$SPF.US.pe.4Q)],US.SPF.data$SPF.US.pe.4Q[!is.na(US.SPF.data$SPF.US.pe.4Q)], col="purple")

Data.US.PE.1y <- US.SPF.data %>%
  dplyr::select(date, SPF.US.pe.1Q, SPF.US.pe.2Q, SPF.US.pe.3Q, SPF.US.pe.4Q) %>%
  mutate(PE.1y = coalesce(SPF.US.pe.1Q, SPF.US.pe.2Q, SPF.US.pe.3Q, SPF.US.pe.4Q))

diff <- Data.US.PE.1y$PE.1y - SPF_US_PE_1_2y_aggregate$PE_1y[52:dim(SPF_US_PE_1_2y_aggregate)[1]]
mean.diff <- mean(diff[-c(length(diff)-c(0,1,2,3))], na.rm=T)

plot(Data.US.PE.1y$date, diff, type="l")
abline(h=0, col="blue")

