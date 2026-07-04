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
## 2. All available SPFs. The data chosen depends on the number of available 
##    observations. For now, since the differences between the means from the 
##    distributions and the means from the provided data are relatively negligible, 
##    we use the means from the data since the series is longer.
## 3. All available variance data. 

#observables.surveys <- SPF.data %>% 
#  dplyr::select(date, SPF.EA.pe.1y, SPF.EA.pe.2y, SPF.EA.pe.5y) 

# Remove 1y. 2022Q4 - 2022-10-15
#survey.DATA[96,2:8] <- NA


#### Reshape Survey.DATA according to maturity 
stat_order <- c("pe", "stdev", "var", "skew", "k3rd", "kurtosis", "k4th")

survey.DATA.G.new <- survey.DATA.G %>%
  pivot_longer(cols = -date, names_to = "var", values_to = "value") %>%
  mutate(
    quarter = quarter(date),
    maturity_q = case_when(
      str_detect(var, "\\.cy\\.")  ~ 5 - quarter,             # 4Q to 1Q
      str_detect(var, "\\.ny\\.")  ~ 5 - quarter + 4,         # 8Q to 5Q
      str_detect(var, "\\.nny\\.") ~ 5 - quarter + 8,         # 12Q to 9Q
      str_detect(var, "\\.5y\\.") & quarter == 1 ~ 20,
      str_detect(var, "\\.5y\\.") & quarter == 2 ~ 19,
      str_detect(var, "\\.5y\\.") & quarter == 3 ~ 22,
      str_detect(var, "\\.5y\\.") & quarter == 4 ~ 21,
      TRUE ~ NA_real_
    ),
    maturity_label = paste0(maturity_q, "Q"),
    new_var = var %>%
      str_replace("\\.cy\\.", paste0(".", maturity_label, ".")) %>%
      str_replace("\\.ny\\.", paste0(".", maturity_label, ".")) %>%
      str_replace("\\.nny\\.", paste0(".", maturity_label, ".")) %>%
      str_replace("\\.5y\\.", paste0(".", maturity_label, "."))
  ) %>%
  dplyr::select(date, new_var, value) %>%
  pivot_wider(names_from = new_var, values_from = value) %>%
  {
    # Reorder columns based on maturity and stat
    cols <- colnames(.)[-1]  # exclude date
    col_info <- tibble(col = cols) %>%
      mutate(
        maturity = str_extract(col, "\\.\\d+Q\\.") %>% str_remove_all("\\.") %>% str_remove("Q") %>% as.integer(),
        stat = str_extract(col, "beta\\.(\\w+)$") %>% str_remove("beta\\.")
      ) %>%
      arrange(maturity, match(stat, stat_order)) %>%
      pull(col)
    
    dplyr::select(., date, all_of(col_info))
  } 



survey.DATA.G.new.with.param <- survey.DATA.G.with.param %>%
  pivot_longer(cols = -date, names_to = "var", values_to = "value") %>%
  mutate(
    quarter = quarter(date),
    maturity_q = case_when(
      str_detect(var, "\\.cy\\.")  ~ 5 - quarter,             # 4Q to 1Q
      str_detect(var, "\\.ny\\.")  ~ 5 - quarter + 4,         # 8Q to 5Q
      str_detect(var, "\\.nny\\.") ~ 5 - quarter + 8,         # 12Q to 9Q
      str_detect(var, "\\.5y\\.") & quarter == 1 ~ 20,
      str_detect(var, "\\.5y\\.") & quarter == 2 ~ 19,
      str_detect(var, "\\.5y\\.") & quarter == 3 ~ 22,
      str_detect(var, "\\.5y\\.") & quarter == 4 ~ 21,
      TRUE ~ NA_real_
    ),
    maturity_label = paste0(maturity_q, "Q"),
    new_var = var %>%
      str_replace("\\.cy\\.", paste0(".", maturity_label, ".")) %>%
      str_replace("\\.ny\\.", paste0(".", maturity_label, ".")) %>%
      str_replace("\\.nny\\.", paste0(".", maturity_label, ".")) %>%
      str_replace("\\.5y\\.", paste0(".", maturity_label, "."))
  ) %>%
  dplyr::select(date, new_var, value) %>%
  pivot_wider(names_from = new_var, values_from = value) %>%
  {
    # Reorder columns based on maturity and stat
    cols <- colnames(.)[-1]  # exclude date
    col_info <- tibble(col = cols) %>%
      mutate(
        maturity = str_extract(col, "\\.\\d+Q\\.") %>% str_remove_all("\\.") %>% str_remove("Q") %>% as.integer(),
        stat = str_extract(col, "beta\\.(\\w+)$") %>% str_remove("beta\\.")
      ) %>%
      arrange(maturity, match(stat, stat_order)) %>%
      pull(col)
    
    dplyr::select(., date, all_of(col_info))
  } 

observables.surveys.EA.G  <- DATA.G %>% 
  dplyr::select(date, EA.growth.quarterly.deseasonalized) %>% 
  rename(EA.growth=EA.growth.quarterly.deseasonalized) %>% 
  filter(date != as.Date("2014-12-15")) %>%
  mutate(date = if_else(year(date) < 2015, 
                        date + months(1),  # Subtract one month
                        date)) %>% 
  full_join(survey.DATA.G.new %>% 
              dplyr::select(date, ends_with(".pe")), by="date")

observables.var.EA.G <- survey.DATA.G.new %>% 
  dplyr::select(date, ends_with(".var")) 

observables.k3rd.EA.G <- survey.DATA.G.new %>% 
  dplyr::select(date, ends_with(".k3rd")) 

observables.k4th.EA.G  <- survey.DATA.G.new %>% 
  dplyr::select(date, ends_with(".k4th")) %>%
  mutate(across(-date, ~ if_else(. > 50, 50, .))) #%>% 
#mutate(across(-date, ~ if_else(. < 0, 0, .)))

observables.with.dates.EA.G <- observables.surveys.EA.G %>% 
  full_join(observables.var.EA.G, by="date") %>% 
  full_join(observables.k3rd.EA.G, by="date") %>%
  full_join(observables.k4th.EA.G, by="date") %>%
  dplyr::filter(date >= "1999-01-15") %>%
  dplyr::select(-matches("\\b((1|2|3)Q|([9]|1[0-9]|2[0-9])Q)\\b"))
#dplyr::filter(date >= "1997-01-15" & date <= "2022-01-15") 

observables.EA.G  <- observables.with.dates.EA.G %>% 
  dplyr::select(-date)


# Select horizons:
H <- matrix(c(12, 15, 18, 21, 24),5,1) # number of columns = number of areas

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

select.inflation.types.EA.G  <- array(0,c(length(c(H)),1,4))
#select.inflation.types[1:length(c(H)),1,2] <- 11 # We have uncertainty measure but no 3rd order moments
#select.inflation.types[1:length(c(H)),1,2] <- 111 # We have uncertainty measure
select.inflation.types.EA.G [1:length(c(H)),1,3] <- 1111 # We have uncertainty measure + 3 and 4 cum


# Make ".EA" matrices:
observables.EA.G <- as_tibble(observables.EA.G )
H.EA.G <- H
select.inflation.types.EA.G  <- select.inflation.types.EA.G 
observables.with.dates.EA.G  <- as_tibble(observables.with.dates.EA.G )
r <- 2 # number of countries


