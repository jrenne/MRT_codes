#===============================================================================
# This program builds a big data.frame object, called DATA,
#   with all data at the monthly frequency
#===============================================================================

library(readxl)


# exctract first digit of the col.names.US.G.SPF vector.


### Extract end of the vector for US
#sort(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",col.names.US.SPF[-1])))))

### Extract end of the vector for EA
#round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(EA.SPF.DISTRI.1y)[-1]))))*2)/2 

# In the end, we have that in:
## FOR FED PHILADELPHIA
## - DATA.US.G contains PE provided by FED Philadelphia with:
###    - SPF.US.G.1y.pe, SPF.US.10y.pe, SPF.US.5y.pe,
###    - US.cpi
###    - US.SPF.G.data (NOT USEFUL => DELETE IT FOR THE MODEL).
## - data.SPF.US.G contains the bins for individual forecasters.
## - US.SPF.G.DISTRI.1Q (8Q) contains the aggregate bins of the SPF survey
##   of FED Philadelphia. Useful for Beta distribution
## - US.SPF.G.data contains data with calculation (pe, stdev, 3rd.cum) made 
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

#### GDP
GDP_raw <- read.csv("./data/US/raw/GDPC1.csv")

GDP <- GDP_raw %>%
  rename(date = DATE) %>%
  mutate(date = as.Date(date, format="%Y-%m-%d") + 14)

US.GDP.aux <- data.frame(date=GDP$date,US.gdp=GDP$GDPC1)

##### HP FILTER TO EXTRACT TREND CYCLE COMPONENT

ts_gdp_us <- ts(US.GDP.aux$US.gdp, frequency=4, start=c(1947,1))

#decomposition trend cycle with HP filter
gdp.decomposition.us <- hpfilter(log(ts_gdp_us), drift=FALSE)
plot(gdp.decomposition.us$cycle)

US.GDP.aux$US.log.gdp <- as.vector(log(ts_gdp_us)*100)
US.GDP.aux$US.log.gdp.cycle <- as.vector(gdp.decomposition.us$cycle*100)
US.GDP.aux$US.log.gdp.trend <- as.vector(gdp.decomposition.us$trend*100)

# Add data to the big data.frame:
DATA.G.US <- US.GDP.aux %>%
  mutate(US.gdp.growth = log(US.gdp/lag(US.gdp,4))*100,
         US.gdp.growth.quarterly = log(US.gdp/lag(US.gdp,1))*100) #%>%
  #filter(date >= "1990-01-15")

#### real GNP as survey use it at the begining
GNP_raw <- read.csv("./data/US/raw/GNPC96.csv")

GNP <- GNP_raw %>%
  rename(date = DATE) %>%
  mutate(date = as.Date(date, format="%Y-%m-%d") + 14)

US.GNP.aux <- data.frame(date=GNP$date,US.gnp=GNP$GNPC96) %>%
  mutate(US.gnp.growth = log(US.gnp/lag(US.gnp,4))*100,
         US.gnp.growth.quarterly = log(US.gnp/lag(US.gnp,1))*100)

#### GNP  as survey use it at the begining
GNP.n_raw <- read.csv("./data/US/raw/GNP.csv")

GNP.n <- GNP.n_raw %>%
  rename(date = DATE) %>%
  mutate(date = as.Date(date, format="%Y-%m-%d") + 14)

US.GNP.n.aux <- data.frame(date=GNP.n$date,US.gnp.n=GNP.n$GNP) %>%
  mutate(US.gnp.n.growth = log(US.gnp.n/lag(US.gnp.n,4))*100,
         US.gnp.n.growth.quarterly = log(US.gnp.n/lag(US.gnp.n,1))*100)

  

#Plot Inflation US data
plot(DATA.G.US$date, DATA.G.US$US.gdp.growth.quarterly,type='l', lwd=1.5)
lines(US.GNP.n.aux$date, US.GNP.n.aux$US.gnp.n.growth.quarterly, col="blue")
lines(US.GNP.aux$date, US.GNP.aux$US.gnp.growth.quarterly, col="red")

plot(DATA.G.US$date, DATA.G.US$US.gdp.growth.quarterly-US.GNP.aux$US.gnp.growth.quarterly, type='l')
plot(DATA.G.US$date, DATA.G.US$US.gdp.growth.quarterly-US.GNP.n.aux$US.gnp.n.growth.quarterly, type='l',col="red")

if(length(DATA.G.US$US.gdp.growth.quarterly)==length(US.GNP.aux$US.gnp.growth.quarterly)){

cor(DATA.G.US$US.gdp.growth.quarterly[-1],US.GNP.aux$US.gnp.growth.quarterly[-1])
cor(DATA.G.US$US.gdp.growth.quarterly[-1],US.GNP.n.aux$US.gnp.n.growth.quarterly[-1])

}
#===========================================================================
# load US SPFs
#===========================================================================



###### PE 1-2 year ahead RGDP
# Contains the level only
SPF_US_G_PE_individual <- read_excel("./data/US/raw/Individual_RGDP.xlsx", 
                                        progress = readxl_progress(), .name_repair = "unique",
                                        na = "#N/A", col_type="numeric") %>%
  mutate(date = as.Date(
    paste("15-",
          QUARTER*3-2,
          "-",
          YEAR,sep=""),
    "%d-%m-%Y")) 

###### VINTAGE PGPD
RGDP_vintages <- read_excel("./data/US/raw/ROUTPUTQvQd.xlsx", 
                            progress = readxl_progress(), .name_repair = "unique",
                            na = "#N/A")

RGDP_vintages_last_4 <- data.frame(date=unique(SPF_US_G_PE_individual$date),
                                   RGDPm5=NA,RGDPm4=NA,RGDPm3=NA,RGDPm2=NA,RGDPm1=NA,RGDP0=NA,RGDP1_bis=NA)
indic.1st.vintage <- grep("ROUTPUT68Q4",colnames(RGDP_vintages))

j=1
for(i in indic.1st.vintage:dim(RGDP_vintages)[2]){
  
  RGDP_vintages_last_4[j,2:8] <- tail(RGDP_vintages[,i][!is.na(RGDP_vintages[,i])],7)
  j = j+1
}
RGDP_vintages_last_4 <- RGDP_vintages_last_4 %>% drop_na()

SPF_US_G_PE_aggregate <- SPF_US_G_PE_individual %>%
  group_by(date) %>%
  dplyr::summarize(YEAR = mean(YEAR, na.rm=TRUE),
                   QUARTER = mean(QUARTER, na.rm=TRUE),
                   RGDP1 = mean(RGDP1, na.rm=TRUE),
                   RGDP2 = mean(RGDP2, na.rm=TRUE),
                   RGDP3 = mean(RGDP3, na.rm=TRUE),
                   RGDP4 = mean(RGDP4, na.rm=TRUE),
                   RGDP5 = mean(RGDP5, na.rm=TRUE),
                   RGDP6 = mean(RGDP6, na.rm=TRUE),
                   RGDPA = mean(RGDPA, na.rm=TRUE),
                   RGDPB = mean(RGDPB, na.rm=TRUE),
                   RGDPC = mean(RGDPC, na.rm=TRUE),
                   RGDPD = mean(RGDPD, na.rm=TRUE)) %>%
  full_join(RGDP_vintages_last_4, by="date") %>%
  mutate(RGDPA_check = case_when(QUARTER == 1 ~ (RGDP2 + RGDP3 + RGDP4 + RGDP5)/4,
                                 QUARTER == 2 ~ (RGDP1 + RGDP2 + RGDP3 + RGDP4)/4,
                                 QUARTER == 3 ~ (RGDP0 + RGDP1 + RGDP2 + RGDP3)/4,
                                 QUARTER == 4 ~ (RGDPm1 + RGDP0 + RGDP1 + RGDP2)/4),
         RGDPAm1 = case_when(QUARTER == 1 ~ (RGDPm2 + RGDPm1 + RGDP0 + RGDP1_bis)/4,
                             QUARTER == 2 ~ (RGDPm3 + RGDPm2 + RGDPm1 + RGDP0)/4,
                             QUARTER == 3 ~ (RGDPm4 + RGDPm3 + RGDPm2 + RGDPm1)/4,
                             QUARTER == 4 ~ (RGDPm5 + RGDPm4 + RGDPm3 + RGDPm2)/4),
         PE_1y = log(RGDPA/RGDPAm1)*100,
         PE_2y = log(RGDPB/RGDPA)*100,
         PE_3y = log(RGDPC/RGDPB)*100,
         PE_4y = log(RGDPD/RGDPC)*100,
         diff_A = abs(RGDPA - RGDPA_check)) %>%
  mutate(PE_1y = replace(PE_1y, which(diff_A>20), NA),
         PE_2y = replace(PE_2y, which(diff_A>20), NA),
         PE_3y = replace(PE_3y, which(diff_A>20), NA),
         PE_4y = replace(PE_4y, which(diff_A>20), NA)) # Problem with 1985-01-15, 1990-01-15, 1986-01-15 => remove them

SPF_US_G_PE_aggregate_reshape <- SPF_US_G_PE_aggregate %>%
  dplyr::select(date,PE_1y, PE_2y, PE_3y, PE_4y) %>%
  filter(date > "1981-04-15") %>%
  tidyr::gather(key=Forecast_Horizon, value =value, PE_1y:PE_4y) %>%
  mutate(Forecast_Horizon = case_when(Forecast_Horizon == "PE_1y" ~ 1,
                                      Forecast_Horizon == "PE_2y" ~ 2,
                                      Forecast_Horizon == "PE_3y" ~ 3,
                                      Forecast_Horizon == "PE_4y" ~ 4)) 


### DEFINE THE NAME OF THE SAMPLE CONSIDERED
US.G.SPF.names <- c("SPF_US_G_individual_1981_1992", "SPF_US_G_individual_1992_2009",	
                    "SPF_US_G_individual_2009_2020",	"SPF_US_G_individual_post_2020")

# Dates to remove
target_dates <- as.Date(c("1985-01-15","1990-01-15", "1986-01-15"),"%Y-%m-%d")

#col.names.US.G.SPF<- c("date", "F6_0",	"F5_0T5_9",	"F4_0T4_9", "F3_0T3_9",	"F2_0T2_9",	
                       #"F1_0T1_9",	"F0_0T0_9",	"FN1_0TN0_1", "FN2_0TN1_1", "TN2_0")


# Contains Bins but no point estimate
SPF_US_G_individual <-read_excel("./data/US/raw/Individual_PRGDP.xlsx", 
                               progress = readxl_progress(), .name_repair = "unique",
                               na = "#N/A", col_type="numeric") %>%
  mutate(date = as.Date(
    paste("15-",
          QUARTER*3-2,
          "-",
          YEAR,sep=""),
    "%d-%m-%Y"))

SPF_US_G_individual_1981_1992 <- SPF_US_G_individual %>%
  filter(date > "1981-04-15" & date <= "1991-10-15") %>%
  rename(Year = YEAR,
         Quarter = QUARTER,
         Forecaster_ID = ID) %>%
  dplyr::select(-c(INDUSTRY, PRGDP13:PRGDP44)) %>%
  tidyr::gather(key=Survey_Bin, value =value, PRGDP1:PRGDP12) %>%
  mutate(Forecast_Horizon = case_when(as.numeric(substr(Survey_Bin,6,7)) < 7 ~ 1,
                                      TRUE ~ 2)) %>% # put 1 if 1-6, 2if 7-12
  mutate(Survey_Bin = case_when(Forecast_Horizon == 2 ~ paste0("PRGDP", as.numeric(substr(Survey_Bin,6,7))-6),
                                TRUE ~ Survey_Bin)) %>%
  spread(key=Survey_Bin, value) %>%
  mutate(F6_0_T7_9 = PRGDP1, #F6_0
         F4_0T5_9 = PRGDP2, 
         F2_0T3_9 = PRGDP3, 
         F0_0T1_9 = PRGDP4, 
         FN2_0TN0_1 = PRGDP5, 
         FN4_0TN2_1 = PRGDP6) %>% #TN2_0
  mutate(Variable_Forecasted = "PRGDP") %>%
  dplyr::select(-starts_with("PRGDP"))  %>%
  dplyr::select(Year,	Quarter,	Forecaster_ID,	Variable_Forecasted,	Forecast_Horizon,
                everything()) %>%
  mutate(across(7:(ncol(.)), ~if_else(date %in% target_dates, NA_real_, .)))

SPF_US_G_individual_1992_2009 <- SPF_US_G_individual %>%
  filter(date > "1991-10-15" & date < "2009-04-15") %>%
  rename(Year = YEAR,
         Quarter = QUARTER,
         Forecaster_ID = ID) %>%
  dplyr::select(-c(INDUSTRY, PRGDP21:PRGDP44)) %>%
  tidyr::gather(key=Survey_Bin, value =value, PRGDP1:PRGDP20) %>%
  mutate(Forecast_Horizon = case_when(as.numeric(substr(Survey_Bin,6,7)) < 11 ~ 1,
                                      TRUE ~ 2)) %>% # put 1 if 1-10, 2if 11-20
  mutate(Survey_Bin = case_when(Forecast_Horizon == 2 ~ paste0("PRGDP", as.numeric(substr(Survey_Bin,6,7))-10),
                                TRUE ~ Survey_Bin)) %>%
  spread(key=Survey_Bin, value) %>%
  mutate(F6_0T6_9 = PRGDP1, #F6_0
         F5_0T5_9 = PRGDP2, 
         F4_0T4_9 = PRGDP3, 
         F3_0T3_9 = PRGDP4, 
         F2_0T2_9 = PRGDP5, 
         F1_0T1_9 = PRGDP6, 
         F0_0T0_9 = PRGDP7, 
         FN1_0TN0_1 = PRGDP8, 
         FN2_0TN1_1 = PRGDP9, 
         FN3_0TN2_1 = PRGDP10) %>% #TN2_0
  mutate(Variable_Forecasted = "PRGDP") %>%
  dplyr::select(-starts_with("PRGDP"))  %>%
  dplyr::select(Year,	Quarter,	Forecaster_ID,	Variable_Forecasted,	Forecast_Horizon,
                everything()) %>%
  mutate(across(7:(ncol(.)), ~if_else(date %in% target_dates, NA_real_, .)))

SPF_US_G_individual_2009_2020 <- SPF_US_G_individual %>%
  filter(date > "2009-01-15" & date < "2020-04-15") %>%
  rename(Year = YEAR,
         Quarter = QUARTER,
         Forecaster_ID = ID) %>%
  dplyr::select(-INDUSTRY) %>%
  tidyr::gather(key=Survey_Bin, value =value, PRGDP1:PRGDP44) %>%
  mutate(Forecast_Horizon = case_when(as.numeric(substr(Survey_Bin,6,8)) < 12 ~ 1,
                                      as.numeric(substr(Survey_Bin,6,8)) > 11 & as.numeric(substr(Survey_Bin,6,8)) < 23 ~ 2,
                                      as.numeric(substr(Survey_Bin,6,8)) > 22 & as.numeric(substr(Survey_Bin,6,8)) < 34 ~ 3,
                                      TRUE ~ 4)) %>%
  mutate(Survey_Bin = case_when(Forecast_Horizon == 4 ~ paste0("PRGDP", as.numeric(substr(Survey_Bin,6,8))-33),
                                Forecast_Horizon == 3 ~ paste0("PRGDP", as.numeric(substr(Survey_Bin,6,8))-22),
                                Forecast_Horizon == 2 ~ paste0("PRGDP", as.numeric(substr(Survey_Bin,6,8))-11),
                                TRUE ~ Survey_Bin)) %>%
  spread(key=Survey_Bin, value) %>%
  mutate(F6_0T6_9 = PRGDP1, #F6_0
         F5_0T5_9 = PRGDP2, 
         F4_0T4_9 = PRGDP3, 
         F3_0T3_9 = PRGDP4, 
         F2_0T2_9 = PRGDP5, 
         F1_0T1_9 = PRGDP6, 
         F0_0T0_9 = PRGDP7, 
         FN1_0TN0_1 = PRGDP8, 
         FN2_0TN1_1 = PRGDP9, 
         FN3_0TN2_1 = PRGDP10, 
         FN4_0TN3_1  = PRGDP11) %>% #TN3_0
  mutate(Variable_Forecasted = "PRGDP") %>%
  dplyr::select(-starts_with("PRGDP")) %>%
  dplyr::select(Year,	Quarter,	Forecaster_ID,	Variable_Forecasted,	Forecast_Horizon,
                everything()) %>%
  mutate(across(7:(ncol(.)), ~if_else(date %in% target_dates, NA_real_, .)))

SPF_US_G_individual_post_2020 <- SPF_US_G_individual %>%
  filter(date > "2020-01-15") %>%
  rename(Year = YEAR,
         Quarter = QUARTER,
         Forecaster_ID = ID) %>%
  dplyr::select(-INDUSTRY) %>%
  tidyr::gather(key=Survey_Bin, value =value, PRGDP1:PRGDP44) %>%
  mutate(Forecast_Horizon = case_when(as.numeric(substr(Survey_Bin,6,8)) < 12 ~ 1,
                                      as.numeric(substr(Survey_Bin,6,8)) > 11 & as.numeric(substr(Survey_Bin,6,8)) < 23 ~ 2,
                                      as.numeric(substr(Survey_Bin,6,8)) > 22 & as.numeric(substr(Survey_Bin,6,8)) < 34 ~ 3,
                                      TRUE ~ 4)) %>%
  mutate(Survey_Bin = case_when(Forecast_Horizon == 4 ~ paste0("PRGDP", as.numeric(substr(Survey_Bin,6,8))-33),
                                Forecast_Horizon == 3 ~ paste0("PRGDP", as.numeric(substr(Survey_Bin,6,8))-22),
                                Forecast_Horizon == 2 ~ paste0("PRGDP", as.numeric(substr(Survey_Bin,6,8))-11),
                                TRUE ~ Survey_Bin)) %>%
  spread(key=Survey_Bin, value) %>%
  mutate(F16_0T18_9 = PRGDP1, #F16_0
         F10_0T15_9 = PRGDP2, 
         F7_0T9_9 = PRGDP3,
         F4_0T6_9 = PRGDP4, 
         F2_5T3_9 = PRGDP5, 
         F1_5T2_4 = PRGDP6, 
         F0_0T1_4 = PRGDP7, 
         FN3_0TN0_1 = PRGDP8, 
         FN6_0TN3_1 = PRGDP9,
         FN12_0TN6_1 = PRGDP10,
         FN15_0TN12_1  = PRGDP11) %>% #TN12_0
  mutate(Variable_Forecasted = "PRGDP") %>%
  dplyr::select(-starts_with("PRGDP")) %>%
  dplyr::select(Year,	Quarter,	Forecaster_ID,	Variable_Forecasted,	Forecast_Horizon,
                everything()) %>%
  mutate(across(7:(ncol(.)), ~if_else(date %in% target_dates, NA_real_, .)))


#data.SPF.US <- rbind(SPF_pre_2014,SPF_post_2014)

US.G.SPF.dates <- SPF_US_G_individual %>% dplyr::select(date) %>% filter(date > "1981-04-15", !duplicated(date))

US.G.SPF.data <- data.frame(date=US.G.SPF.dates,
                          SPF.US.G.pe.1Q=NaN,SPF.US.G.pe.2Q=NaN,SPF.US.G.pe.3Q=NaN,SPF.US.G.pe.4Q=NaN,
                          SPF.US.G.pe.5Q=NaN,SPF.US.G.pe.6Q=NaN,SPF.US.G.pe.7Q=NaN,SPF.US.G.pe.8Q=NaN,
                          SPF.US.G.pe.9Q=NaN,SPF.US.G.pe.10Q=NaN,SPF.US.G.pe.11Q=NaN,SPF.US.G.pe.12Q=NaN,
                          SPF.US.G.pe.13Q=NaN,SPF.US.G.pe.14Q=NaN,SPF.US.G.pe.15Q=NaN,SPF.US.G.pe.16Q=NaN,
                          SPF.US.G.stdev.1Q=NaN,SPF.US.G.stdev.2Q=NaN,SPF.US.G.stdev.3Q=NaN,SPF.US.G.stdev.4Q=NaN,
                          SPF.US.G.stdev.5Q=NaN,SPF.US.G.stdev.6Q=NaN,SPF.US.G.stdev.7Q=NaN,SPF.US.G.stdev.8Q=NaN,
                          SPF.US.G.stdev.9Q=NaN,SPF.US.G.stdev.10Q=NaN,SPF.US.G.stdev.11Q=NaN,SPF.US.G.stdev.12Q=NaN,
                          SPF.US.G.stdev.13Q=NaN,SPF.US.G.stdev.14Q=NaN,SPF.US.G.stdev.15Q=NaN,SPF.US.G.stdev.16Q=NaN,
                          SPF.US.G.disagreement.1Q=NaN,SPF.US.G.disagreement.2Q=NaN,SPF.US.G.disagreement.3Q=NaN,SPF.US.G.disagreement.4Q=NaN,
                          SPF.US.G.disagreement.5Q=NaN,SPF.US.G.disagreement.6Q=NaN,SPF.US.G.disagreement.7Q=NaN,SPF.US.G.disagreement.8Q=NaN,
                          SPF.US.G.disagreement.9Q=NaN,SPF.US.G.disagreement.10Q=NaN,SPF.US.G.disagreement.11Q=NaN,SPF.US.G.disagreement.12Q=NaN,
                          SPF.US.G.disagreement.13Q=NaN,SPF.US.G.disagreement.14Q=NaN,SPF.US.G.disagreement.15Q=NaN,SPF.US.G.disagreement.16Q=NaN,
                          SPF.US.G.thirdcumulant.1Q=NaN,SPF.US.G.thirdcumulant.2Q=NaN,SPF.US.G.thirdcumulant.3Q=NaN,SPF.US.G.thirdcumulant.4Q=NaN,
                          SPF.US.G.thirdcumulant.5Q=NaN,SPF.US.G.thirdcumulant.6Q=NaN,SPF.US.G.thirdcumulant.7Q=NaN,SPF.US.G.thirdcumulant.8Q=NaN,
                          SPF.US.G.thirdcumulant.9Q=NaN,SPF.US.G.thirdcumulant.10Q=NaN,SPF.US.G.thirdcumulant.11Q=NaN,SPF.US.G.thirdcumulant.12Q=NaN,
                          SPF.US.G.thirdcumulant.13Q=NaN,SPF.US.G.thirdcumulant.14Q=NaN,SPF.US.G.thirdcumulant.15Q=NaN,SPF.US.G.thirdcumulant.16Q=NaN
)

count <- 0
for(spf.sample in US.G.SPF.names){
  
  # Create data.SPF.US that contains the considered sample
  eval(parse(text = gsub(" ","",paste("data.SPF.US.G <- ",
                                      spf.sample, sep="")))) 
  
  data.SPF.US.G$check.na <- !is.na(data.SPF.US.G[,7])
  data.SPF.US.G <- data.SPF.US.G %>%
    dplyr::select(Year,	Quarter,	Forecaster_ID,	Variable_Forecasted,	Forecast_Horizon, date, check.na,
                  everything())
  
  lower.bound.interval.initial <- round(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",colnames(data.SPF.US.G)[8:dim(data.SPF.US.G)[2]]))))*2)/2
  upper.bound.interval.initial <- round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(data.SPF.US.G)[8:dim(data.SPF.US.G)[2]]))))*2)/2 
  all.intervals.initial <- upper.bound.interval.initial - lower.bound.interval.initial

  med.classes.US.G.initial <- lower.bound.interval.initial + all.intervals.initial/2
  #med.classes.US <- seq(8.25,-.25,by=-.5)# note that quantiles in reverse order
  
  US.G.SPF.dates.indic <- data.SPF.US.G %>% dplyr::select(date) %>% filter(!duplicated(date))
  
  nbr.bins.pre <-length(colnames(data.SPF.US.G))-7
  
  # The distributions will be saved in dataframes US.SPF.DISTRI.XQ
  US.G.SPF.DISTRI.1Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.1Q <- cbind(US.G.SPF.DISTRI.1Q,matrix(NaN,dim(US.G.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.1Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.2Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.2Q <- cbind(US.G.SPF.DISTRI.2Q,matrix(NaN,dim(US.G.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.2Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.3Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.3Q <- cbind(US.G.SPF.DISTRI.3Q,matrix(NaN,dim(US.G.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.3Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.4Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.4Q <- cbind(US.G.SPF.DISTRI.4Q,matrix(NaN,dim(US.G.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.4Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.5Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.5Q <- cbind(US.G.SPF.DISTRI.5Q,matrix(NaN,dim(US.G.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.5Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.6Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.6Q <- cbind(US.G.SPF.DISTRI.6Q,matrix(NaN,dim(US.G.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.6Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.7Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.7Q <- cbind(US.G.SPF.DISTRI.7Q,matrix(NaN,dim(US.G.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.7Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.8Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.8Q <- cbind(US.G.SPF.DISTRI.8Q,matrix(NaN,dim(US.G.SPF.DISTRI.1Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.8Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.9Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.9Q <- cbind(US.G.SPF.DISTRI.9Q,matrix(NaN,dim(US.G.SPF.DISTRI.9Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.9Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.10Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.10Q <- cbind(US.G.SPF.DISTRI.10Q,matrix(NaN,dim(US.G.SPF.DISTRI.10Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.10Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.11Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.11Q <- cbind(US.G.SPF.DISTRI.11Q,matrix(NaN,dim(US.G.SPF.DISTRI.11Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.11Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.12Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.12Q <- cbind(US.G.SPF.DISTRI.12Q,matrix(NaN,dim(US.G.SPF.DISTRI.12Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.12Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.13Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.13Q <- cbind(US.G.SPF.DISTRI.13Q,matrix(NaN,dim(US.G.SPF.DISTRI.13Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.13Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.14Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.14Q <- cbind(US.G.SPF.DISTRI.14Q,matrix(NaN,dim(US.G.SPF.DISTRI.14Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.14Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.15Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.15Q <- cbind(US.G.SPF.DISTRI.15Q,matrix(NaN,dim(US.G.SPF.DISTRI.15Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.15Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  US.G.SPF.DISTRI.16Q <- data.frame(date=US.G.SPF.dates.indic$date)
  US.G.SPF.DISTRI.16Q <- cbind(US.G.SPF.DISTRI.16Q,matrix(NaN,dim(US.G.SPF.DISTRI.16Q)[1],nbr.bins.pre))
  colnames(US.G.SPF.DISTRI.16Q) <- colnames(data.SPF.US.G)[c(6,8:dim(data.SPF.US.G)[2])]
  
  for(y in US.G.SPF.dates.indic$date){
    indic.month <- which(US.G.SPF.dates.indic$date==c(y))
    count <- count + 1
    print(as.Date(y))
    for(horizon in 1:2){
      data.aux <- subset(data.SPF.US.G,(date==c(y))&(Forecast_Horizon==horizon)&check.na)
      avg.distri <- apply(data.aux[,8:dim(data.aux)[2]],2,mean,na.rm=T)/100
      
      # Modify first or last bin if one of the two if > 5% (take the larger of the two)
      if(max(avg.distri[1], avg.distri[length(avg.distri)])>0.1 & !is.na(max(avg.distri[1], avg.distri[length(avg.distri)]))){
        
        # Extract PE data
        data.PE <- SPF_US_G_PE_aggregate_reshape %>% filter(date==c(y)& (Forecast_Horizon==horizon))
        # Identify if we need to change the first or the last bins
        max.prob <- max(avg.distri[1], avg.distri[length(avg.distri)])
        indic.pos.max.prob <- match(max.prob, avg.distri)
        
        # Compute the interval that allows to fit perfectly the PE
        new.min.max.interval <- (data.PE$value - sum(((avg.distri)*(med.classes.US.G.initial))[-indic.pos.max.prob]))/max.prob
        
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
          med.classes.US.G <- lower.bound.interval + all.intervals/2
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
          med.classes.US.G <- lower.bound.interval.initial + all.intervals/2
        }
        
      } else{ #Condition to consider initial value
        med.classes.US.G <- med.classes.US.G.initial
      }
      
      mean.per.forecaster <- apply((data.aux[,8:dim(data.aux)[2]]/100) * (matrix(1,dim(data.aux)[1],1) %*% med.classes.US.G),1,sum)
      disagreement.aux <- sd(mean.per.forecaster)
      mean.aux <- sum((avg.distri)*(med.classes.US.G))
      stdev.aux <- sqrt(sum((avg.distri)*(med.classes.US.G)^2) - mean.aux^2)
      thirdcumulant.aux <- sum((avg.distri)*(med.classes.US.G - mean.aux)^3)
      # Simple average of forecasts??... potential pbm: not the same forecasters for horizons 1 and 2
      if(horizon==1){
        if(data.aux$Quarter[1]==1 & !is.na(data.aux$Quarter[1])){
          US.G.SPF.data$SPF.US.G.pe.4Q[count]    <- mean.aux
          US.G.SPF.data$SPF.US.G.stdev.4Q[count] <- stdev.aux
          US.G.SPF.data$SPF.US.G.thirdcumulant.4Q[count] <- thirdcumulant.aux
          US.G.SPF.data$SPF.US.G.disagreement.4Q[count] <- disagreement.aux
          US.G.SPF.DISTRI.4Q[indic.month,2:dim(US.G.SPF.DISTRI.4Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==2 & !is.na(data.aux$Quarter[1])){
          US.G.SPF.data$SPF.US.G.pe.3Q[count]    <- mean.aux
          US.G.SPF.data$SPF.US.G.stdev.3Q[count] <- stdev.aux
          US.G.SPF.data$SPF.US.G.thirdcumulant.3Q[count] <- thirdcumulant.aux
          US.G.SPF.data$SPF.US.G.disagreement.3Q[count] <- disagreement.aux
          US.G.SPF.DISTRI.3Q[indic.month,2:dim(US.G.SPF.DISTRI.3Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==3 & !is.na(data.aux$Quarter[1])){
          US.G.SPF.data$SPF.US.G.pe.2Q[count]    <- mean.aux
          US.G.SPF.data$SPF.US.G.stdev.2Q[count] <- stdev.aux
          US.G.SPF.data$SPF.US.G.thirdcumulant.2Q[count] <- thirdcumulant.aux
          US.G.SPF.data$SPF.US.G.disagreement.2Q[count] <- disagreement.aux
          US.G.SPF.DISTRI.2Q[indic.month,2:dim(US.G.SPF.DISTRI.2Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==4 & !is.na(data.aux$Quarter[1])){
          US.G.SPF.data$SPF.US.G.pe.1Q[count]    <- mean.aux
          US.G.SPF.data$SPF.US.G.stdev.1Q[count] <- stdev.aux
          US.G.SPF.data$SPF.US.G.thirdcumulant.1Q[count] <- thirdcumulant.aux
          US.G.SPF.data$SPF.US.G.disagreement.1Q[count] <- disagreement.aux
          US.G.SPF.DISTRI.1Q[indic.month,2:dim(US.G.SPF.DISTRI.1Q)[2]] <- avg.distri
        }
      }else{#horizon == 2
        if(data.aux$Quarter[1]==1 & !is.na(data.aux$Quarter[1])){
          US.G.SPF.data$SPF.US.G.pe.8Q[count]    <- mean.aux
          US.G.SPF.data$SPF.US.G.stdev.8Q[count] <- stdev.aux
          US.G.SPF.data$SPF.US.G.disagreement.8Q[count] <- disagreement.aux
          US.G.SPF.data$SPF.US.G.thirdcumulant.8Q[count] <- thirdcumulant.aux
          US.G.SPF.DISTRI.8Q[indic.month,2:dim(US.G.SPF.DISTRI.8Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==2 & !is.na(data.aux$Quarter[1])){
          US.G.SPF.data$SPF.US.G.pe.7Q[count]    <- mean.aux
          US.G.SPF.data$SPF.US.G.stdev.7Q[count] <- stdev.aux
          US.G.SPF.data$SPF.US.G.thirdcumulant.7Q[count] <- thirdcumulant.aux
          US.G.SPF.data$SPF.US.G.disagreement.7Q[count] <- disagreement.aux
          US.G.SPF.DISTRI.7Q[indic.month,2:dim(US.G.SPF.DISTRI.7Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==3 & !is.na(data.aux$Quarter[1])){
          US.G.SPF.data$SPF.US.G.pe.6Q[count]    <- mean.aux
          US.G.SPF.data$SPF.US.G.stdev.6Q[count] <- stdev.aux
          US.G.SPF.data$SPF.US.G.thirdcumulant.6Q[count] <- thirdcumulant.aux
          US.G.SPF.data$SPF.US.G.disagreement.6Q[count] <- disagreement.aux
          US.G.SPF.DISTRI.6Q[indic.month,2:dim(US.G.SPF.DISTRI.6Q)[2]] <- avg.distri
        }else if(data.aux$Quarter[1]==4 & !is.na(data.aux$Quarter[1])){
          US.G.SPF.data$SPF.US.G.pe.5Q[count]    <- mean.aux
          US.G.SPF.data$SPF.US.G.stdev.5Q[count] <- stdev.aux
          US.G.SPF.data$SPF.US.G.thirdcumulant.5Q[count] <- thirdcumulant.aux
          US.G.SPF.data$SPF.US.G.disagreement.5Q[count] <- disagreement.aux
          US.G.SPF.DISTRI.5Q[indic.month,2:dim(US.G.SPF.DISTRI.5Q)[2]] <- avg.distri
        }
      }
    }
    
     if(y > as.Date("2009-01-15")){
      for(horizon in 3:4){
        data.aux <- subset(data.SPF.US.G,(date==c(y))&(Forecast_Horizon==horizon)&check.na)
        avg.distri <- apply(data.aux[,8:dim(data.aux)[2]],2,mean, na.rm=T)/100
        
        # Modify first or last bin if one of the two if > 5% (take the larger of the two)
        if(max(avg.distri[1], avg.distri[length(avg.distri)])>0.1 & !is.na(max(avg.distri[1], avg.distri[length(avg.distri)]))){
          
          # Extract PE data
          data.PE <- SPF_US_G_PE_aggregate_reshape %>% filter(date==c(y)& (Forecast_Horizon==horizon))
          # Identify if we need to change the first or the last bins
          max.prob <- max(avg.distri[1], avg.distri[length(avg.distri)])
          indic.pos.max.prob <- match(max.prob, avg.distri)
          
          # Compute the interval that allows to fit perfectly the PE
          new.min.max.interval <- (data.PE$value - sum(((avg.distri)*(med.classes.US.G.initial))[-indic.pos.max.prob]))/max.prob
          
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
            med.classes.US.G <- lower.bound.interval + all.intervals/2
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
            med.classes.US.G <- lower.bound.interval.initial + all.intervals/2
          }
          
        } else{ #Condition to consider initial value
          med.classes.US.G <- med.classes.US.G.initial
        }
        
        mean.per.forecaster <- apply((data.aux[,8:dim(data.aux)[2]]/100) * (matrix(1,dim(data.aux)[1],1) %*% med.classes.US.G),1,sum)
        disagreement.aux <- sd(mean.per.forecaster)
        mean.aux <- sum((avg.distri)*(med.classes.US.G))
        stdev.aux <- sqrt(sum((avg.distri)*(med.classes.US.G)^2) - mean.aux^2)
        thirdcumulant.aux <- sum((avg.distri)*(med.classes.US.G - mean.aux)^3)
        # Simple average of forecasts??... potential pbm: not the same forecasters for horizons 1 and 2
        if(horizon==3){
          if(data.aux$Quarter[1]==1 & !is.na(data.aux$Quarter[1])){
            US.G.SPF.data$SPF.US.G.pe.12Q[count]    <- mean.aux
            US.G.SPF.data$SPF.US.G.stdev.12Q[count] <- stdev.aux
            US.G.SPF.data$SPF.US.G.thirdcumulant.12Q[count] <- thirdcumulant.aux
            US.G.SPF.data$SPF.US.G.disagreement.12Q[count] <- disagreement.aux
            US.G.SPF.DISTRI.12Q[indic.month,2:dim(US.G.SPF.DISTRI.12Q)[2]] <- avg.distri
          }else if(data.aux$Quarter[1]==2 & !is.na(data.aux$Quarter[1])){
            US.G.SPF.data$SPF.US.G.pe.11Q[count]    <- mean.aux
            US.G.SPF.data$SPF.US.G.stdev.11Q[count] <- stdev.aux
            US.G.SPF.data$SPF.US.G.thirdcumulant.11Q[count] <- thirdcumulant.aux
            US.G.SPF.data$SPF.US.G.disagreement.11Q[count] <- disagreement.aux
            US.G.SPF.DISTRI.11Q[indic.month,2:dim(US.G.SPF.DISTRI.11Q)[2]] <- avg.distri
          }else if(data.aux$Quarter[1]==3 & !is.na(data.aux$Quarter[1])){
            US.G.SPF.data$SPF.US.G.pe.10Q[count]    <- mean.aux
            US.G.SPF.data$SPF.US.G.stdev.10Q[count] <- stdev.aux
            US.G.SPF.data$SPF.US.G.thirdcumulant.10Q[count] <- thirdcumulant.aux
            US.G.SPF.data$SPF.US.G.disagreement.10Q[count] <- disagreement.aux
            US.G.SPF.DISTRI.10Q[indic.month,2:dim(US.G.SPF.DISTRI.10Q)[2]] <- avg.distri
          }else if(data.aux$Quarter[1]==4 & !is.na(data.aux$Quarter[1])){
            US.G.SPF.data$SPF.US.G.pe.9Q[count]    <- mean.aux
            US.G.SPF.data$SPF.US.G.stdev.9Q[count] <- stdev.aux
            US.G.SPF.data$SPF.US.G.thirdcumulant.9Q[count] <- thirdcumulant.aux
            US.G.SPF.data$SPF.US.G.disagreement.9Q[count] <- disagreement.aux
            US.G.SPF.DISTRI.9Q[indic.month,2:dim(US.G.SPF.DISTRI.9Q)[2]] <- avg.distri
          }
        }else{#horizon == 4
          if(data.aux$Quarter[1]==1 & !is.na(data.aux$Quarter[1])){
            US.G.SPF.data$SPF.US.G.pe.16Q[count]    <- mean.aux
            US.G.SPF.data$SPF.US.G.stdev.16Q[count] <- stdev.aux
            US.G.SPF.data$SPF.US.G.disagreement.16Q[count] <- disagreement.aux
            US.G.SPF.data$SPF.US.G.thirdcumulant.16Q[count] <- thirdcumulant.aux
            US.G.SPF.DISTRI.16Q[indic.month,2:dim(US.G.SPF.DISTRI.16Q)[2]] <- avg.distri
          }else if(data.aux$Quarter[1]==2 & !is.na(data.aux$Quarter[1])){
            US.G.SPF.data$SPF.US.G.pe.15Q[count]    <- mean.aux
            US.G.SPF.data$SPF.US.G.stdev.15Q[count] <- stdev.aux
            US.G.SPF.data$SPF.US.G.thirdcumulant.15Q[count] <- thirdcumulant.aux
            US.G.SPF.data$SPF.US.G.disagreement.15Q[count] <- disagreement.aux
            US.G.SPF.DISTRI.15Q[indic.month,2:dim(US.G.SPF.DISTRI.15Q)[2]] <- avg.distri
          }else if(data.aux$Quarter[1]==3 & !is.na(data.aux$Quarter[1])){
            US.G.SPF.data$SPF.US.G.pe.14Q[count]    <- mean.aux
            US.G.SPF.data$SPF.US.G.stdev.14Q[count] <- stdev.aux
            US.G.SPF.data$SPF.US.G.thirdcumulant.14Q[count] <- thirdcumulant.aux
            US.G.SPF.data$SPF.US.G.disagreement.14Q[count] <- disagreement.aux
            US.G.SPF.DISTRI.14Q[indic.month,2:dim(US.G.SPF.DISTRI.14Q)[2]] <- avg.distri
          }else if(data.aux$Quarter[1]==4 & !is.na(data.aux$Quarter[1])){
            US.G.SPF.data$SPF.US.G.pe.13Q[count]    <- mean.aux
            US.G.SPF.data$SPF.US.G.stdev.13Q[count] <- stdev.aux
            US.G.SPF.data$SPF.US.G.thirdcumulant.13Q[count] <- thirdcumulant.aux
            US.G.SPF.data$SPF.US.G.disagreement.13Q[count] <- disagreement.aux
            US.G.SPF.DISTRI.13Q[indic.month,2:dim(US.G.SPF.DISTRI.13Q)[2]] <- avg.distri
          }
        }
      }
    }

    print(sum(avg.distri))
  }
  
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.1Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.1Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.2Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.2Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.3Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.3Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.4Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.4Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.5Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.5Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.6Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.6Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.7Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.7Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.8Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.8Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.9Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.9Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.10Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.10Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.11Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.11Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.12Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.12Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.13Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.13Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.14Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.14Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.15Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.15Q", sep="")))) 
  eval(parse(text = gsub(" ","",paste("US.G.SPF.DISTRI.16Q.", str_sub(spf.sample,-9,-1), "<- US.G.SPF.DISTRI.16Q", sep="")))) 
  
}

# Add data to the big data.frame:
DATA.G.US <- merge(DATA.G.US, US.G.SPF.data,by="date",all=TRUE)




plot(US.G.SPF.data$date,US.G.SPF.data$SPF.US.G.stdev.1Q,ylim=c(0,2))
points(US.G.SPF.data$date,US.G.SPF.data$SPF.US.G.stdev.3Q,col='blue')
points(US.G.SPF.data$date,US.G.SPF.data$SPF.US.G.stdev.8Q,col='red')


plot(SPF_US_G_PE_aggregate$date,SPF_US_G_PE_aggregate$PE_1y, type="l", lwd=2)
points(US.G.SPF.data$date[!is.na(US.G.SPF.data$SPF.US.G.pe.1Q)],US.G.SPF.data$SPF.US.G.pe.1Q[!is.na(US.G.SPF.data$SPF.US.G.pe.1Q)], col="blue")
points(US.G.SPF.data$date[!is.na(US.G.SPF.data$SPF.US.G.pe.2Q)],US.G.SPF.data$SPF.US.G.pe.2Q[!is.na(US.G.SPF.data$SPF.US.G.pe.2Q)], col="red")
points(US.G.SPF.data$date[!is.na(US.G.SPF.data$SPF.US.G.pe.3Q)],US.G.SPF.data$SPF.US.G.pe.3Q[!is.na(US.G.SPF.data$SPF.US.G.pe.3Q)], col="green")
points(US.G.SPF.data$date[!is.na(US.G.SPF.data$SPF.US.G.pe.4Q)],US.G.SPF.data$SPF.US.G.pe.4Q[!is.na(US.G.SPF.data$SPF.US.G.pe.4Q)], col="purple")

Data.US.G.PE.1y <- US.G.SPF.data %>%
  dplyr::select(date, SPF.US.G.pe.1Q, SPF.US.G.pe.2Q, SPF.US.G.pe.3Q, SPF.US.G.pe.4Q) %>%
  mutate(PE.1y = coalesce(SPF.US.G.pe.1Q, SPF.US.G.pe.2Q, SPF.US.G.pe.3Q, SPF.US.G.pe.4Q))

diff <- Data.US.G.PE.1y$PE.1y - SPF_US_G_PE_aggregate$PE_1y[52:dim(SPF_US_G_PE_aggregate)[1]]
mean.diff <- mean(diff[-c(length(diff)-c(0,1,2,3))], na.rm=T)

plot(Data.US.G.PE.1y$date, diff, type="l")
abline(h=0, col="blue")

plot(SPF_US_G_PE_aggregate$date,SPF_US_G_PE_aggregate$PE_1y, type="l", lwd=2)
lines(Data.US.G.PE.1y$date, Data.US.G.PE.1y$PE.1y, col="red")




#===========================================================================
# load US PDSs
#===========================================================================


