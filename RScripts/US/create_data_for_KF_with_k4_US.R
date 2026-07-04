# ============================================
#        * * * * PAPER INFLATION - PHD * * * *
# ============================================
#           
# --------------------------------------------
#           Adrien Jean-Paul TSCHOPP 
# ============================================
#   CREATE DATA FOR KF 
# --------------------------------------------
#                 * ALL *
# ============================================

# 
# This program creates the matrix with all observables
#

#The vector of observables contain:
## 1. Date
## 2. PDS data ordered by maturities (1st definition of inflation)
##    - First PE, then Var and finally k3rd
## 3. US SPF data ordered by maturities (3rd definition). We exclude 1Q, 2Q, 3Q 
##    as the 3rd definition involves the last 4 quarter => use known info.
##    - First PE, then Var and finally k3rd.


# data.plot <- survey.DATA %>%
#   dplyr::select(date, PDS.5y.avg.beta.pe) %>%
#   na.omit()
# 
# data.plot1 <- DATA.US  %>%
#   dplyr::select(date, SPF.US.5y.pe) %>%
#   na.omit()
# 
# plot(data.plot1$date, data.plot1$SPF.US.5y.pe,  
#      col="red", ylim = c(1,4), type="l")
# lines(data.plot$date,data.plot$PDS.5y.avg.beta.pe)

observables.surveys.US <- DATA.US %>% 
  dplyr::select(date, US.infl.quarterly) %>% 
  rename(US.infl=US.infl.quarterly) %>%
  full_join(survey.DATA.US %>% 
              dplyr::select(date, ends_with(".pe")), by="date") %>% 
  mutate(PDS.5y.avg.beta.pe = coalesce(PDS.5y.avg.beta.pe, PDS.5y.beta.pe), #fills the NA from the first vector with values from the second vector at corresponding positions
         PDS.5y5y.avg.beta.pe = coalesce(PDS.5y5y.avg.beta.pe, PDS.5_10y.beta.pe)) %>% 
  dplyr::select(-c(PDS.5y.beta.pe, PDS.5_10y.beta.pe, contains(c("1Q", "2Q", "3Q"))))

observables.var.US <- survey.DATA.US %>% 
  dplyr::select(date, ends_with(".var")) %>% 
  mutate(PDS.5y.avg.beta.var = coalesce(PDS.5y.avg.beta.var, PDS.5y.beta.var), #fills the NA from the first vector with values from the second vector at corresponding positions
         PDS.5y5y.avg.beta.var = coalesce(PDS.5y5y.avg.beta.var, PDS.5_10y.beta.var)) %>% 
  dplyr::select(-c(PDS.5y.beta.var, PDS.5_10y.beta.var, contains(c("1Q", "2Q", "3Q"))))

observables.k3rd.US <- survey.DATA.US %>% 
  dplyr::select(date, ends_with(".k3rd")) %>% 
  mutate(PDS.5y.avg.beta.k3rd = coalesce(PDS.5y.avg.beta.k3rd, PDS.5y.beta.k3rd), #fills the NA from the first vector with values from the second vector at corresponding positions
         PDS.5y5y.avg.beta.k3rd = coalesce(PDS.5y5y.avg.beta.k3rd, PDS.5_10y.beta.k3rd)) %>% 
  dplyr::select(-c(PDS.5y.beta.k3rd, PDS.5_10y.beta.k3rd, contains(c("1Q", "2Q", "3Q"))))

observables.k4th.US <- survey.DATA.US %>% 
  dplyr::select(date, ends_with(".k4th")) %>% 
  mutate(PDS.5y.avg.beta.k4th = coalesce(PDS.5y.avg.beta.k4th, PDS.5y.beta.k4th), #fills the NA from the first vector with values from the second vector at corresponding positions
         PDS.5y5y.avg.beta.k4th = coalesce(PDS.5y5y.avg.beta.k4th, PDS.5_10y.beta.k4th)) %>% 
  mutate(across(-date, ~ if_else(. < 0, 0, .))) %>% 
  mutate(across(-date, ~ if_else(. > 50, 50, .))) %>% 
  dplyr::select(-c(PDS.5y.beta.k4th, PDS.5_10y.beta.k4th, contains(c("1Q", "2Q", "3Q"))))


observables.with.dates.US <- observables.surveys.US %>% 
  dplyr::select(-c(contains("SPF"))) %>% 
  full_join(observables.var.US %>% dplyr::select(-c(contains("SPF"))), by="date") %>% 
  full_join(observables.k3rd.US %>% dplyr::select(-c(contains("SPF"))), by="date") %>%
  full_join(observables.k4th.US %>% dplyr::select(-c(contains("SPF"))), by="date") %>%
  full_join(observables.surveys.US %>% dplyr::select(date, contains("SPF")), by="date") %>%
  full_join(observables.var.US %>% dplyr::select(date, contains("SPF")), by="date") %>% 
  full_join(observables.k3rd.US %>% dplyr::select(date, contains("SPF")), by="date") %>%
  full_join(observables.k4th.US %>% dplyr::select(date, contains("SPF")), by="date") %>%
  dplyr::filter(date >= "1981-07-15") 
  #dplyr::filter(date >= "1997-01-15" & date <= "2022-01-15") 

observables.US <- observables.with.dates.US %>% 
  dplyr::select(-date)


# Select horizons:
H.US <- matrix(c(12, 15, 18, 21, 24, 60, 120),7,1) # number of columns = number of areas

# Create the vector "select.inflation.types".
## select.inflation.types is an array of dimension:
##  dim(1) = (number of horizs) 
##  dim(2) = (number of areas) x
##  dim(3) = (4 = number of inflation types). Technically there is only one type
##           of inflation in Switzerland. But with keep 4 in order to remain as
##           flexible as possible if one need to extend to other countries.
## "1" indicates that point estimate (only) is available for this horiz and type of inflation
## "11" indicates that point estimate AND variance are available for this horiz and type of inflation
## "111" indicates that point estimate AND variance AND 3rd order cumulants are available for this horiz and type of inflation

select.inflation.types.US <- array(0,c(length(c(H.US)),1,4))
select.inflation.types.US[6:length(c(H.US)),1,1] <- 1111 # We have uncertainty, skew and kurtosis measures 
select.inflation.types.US[1:5,1,3] <- 1111 # We have uncertainty, skew and kurtosis measures


# Make ".US" matrices:
observables.US <- as_tibble(observables.US)
select.inflation.types.US <- select.inflation.types.US
observables.with.dates.US <- as_tibble(observables.with.dates.US)
r <- 1 # number of countries
