library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# load epidemic datasets with outbreak indicators
FluNet_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/FluNet_epidemic.csv")
FluNet_epidemic$SDATE <- as.POSIXct(FluNet_epidemic$SDATE)
HM_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HM_epidemic_ROC.csv") %>%
  select(-"X")
HM_epidemic$date <- as.POSIXct(HM_epidemic$date)
EIOS_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_epidemic_ROC.csv") %>%
  select(-"X")
EIOS_epidemic$date <- as.POSIXct(EIOS_epidemic$date)

country_list <- levels(FluNet_epidemic$Country)

# revise startend column of FluNet
FluNet_epidemic$startend <- NA
for(i in 2:nrow(FluNet_epidemic)){
  if(FluNet_epidemic$epidemic[i] == TRUE & FluNet_epidemic$epidemic[i-1] == FALSE){
    FluNet_epidemic$startend[i] <- "start"
  }else if(FluNet_epidemic$epidemic[i] == FALSE & FluNet_epidemic$epidemic[i-1] == TRUE){
    FluNet_epidemic$startend[i] <- "end"
  }
}
FluNet_epidemic$startend[which(FluNet_epidemic$SDATE == min(FluNet_epidemic$SDATE))] <- NA
# for Saudi Arabia (data start later)
sa <- which(FluNet_epidemic$Country == "Saudi Arabia")
FluNet_epidemic$startend[sa[1]] <- NA

# calculate summary information on outbreaks
FluNet_outbreak_length <- FluNet_epidemic %>% filter(epidemic == TRUE) %>% group_by(Country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n(), start_date = min(SDATE), end_date = max(SDATE), height = max(ALL_INF)) %>% 
  mutate(outbreak_no = row_number()) %>% select(-grp)
FluNet_outbreak_length$outbreak_interval <- interval(FluNet_outbreak_length$start_date, FluNet_outbreak_length$end_date)


### sensitivity per outbreak
# HealthMap
HM_sens_outbreak_list <- vector(mode = "list")

for(i in seq_along(country_list)){
  detected_list <- vector(mode = "list")
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  FluNet_HM_temp <- FluNet_HM_temp %>% group_by(Country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
  
  # state vector for every outbreak if detected or not
  for(j in 1:max(FluNet_HM_temp$outbreak_period)){
    FluNet_HM_temp_ob <- filter(FluNet_HM_temp, outbreak_period == j)
    HM_temp <- filter(HM_epidemic, country == country_list[i] & 
                        date %within% interval(min(FluNet_HM_temp_ob$SDATE), max(FluNet_HM_temp_ob$SDATE)))
    
    detected_list[[j]] <- sapply(57:73, function(x){
      state <- FluNet_HM_temp_ob$epidemic
      alarm <- HM_temp[x]
      return(sum(state[alarm == TRUE], na.rm = TRUE) > 0)
    })
  } 
  
  # transform list to df, sum up columns and calculate sens per outbreak
  detected_df <- data.frame(matrix(unlist(detected_list), nrow = max(FluNet_HM_temp$outbreak_period), byrow = TRUE), 
                            stringsAsFactors = FALSE)
  HM_sens_outbreak_list[[i]] <- apply(detected_df, 2, sum)/max(FluNet_HM_temp$outbreak_period)
}

HM_sens_outbreak <- data.frame(matrix(unlist(HM_sens_outbreak_list), ncol = 17, byrow = TRUE))
names(HM_sens_outbreak) <- paste("bcp_", seq(0.1, 0.9, 0.05), sep = "")
HM_sens_outbreak$country <- country_list
HM_sens_outbreak <- select(HM_sens_outbreak, country, everything())


# EIOS
EIOS_sens_outbreak_list <- vector(mode = "list")

for(i in seq_along(country_list)){
  detected_list <- vector(mode = "list")
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_EIOS_temp <- FluNet_EIOS_temp %>% group_by(Country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
  
  # state vector for every outbreak if detected or not
  for(j in 1:max(FluNet_EIOS_temp$outbreak_period)){
    FluNet_EIOS_temp_ob <- filter(FluNet_EIOS_temp, outbreak_period == j)
    EIOS_temp <- filter(EIOS_epidemic, country == country_list[i] & 
                        date %within% interval(min(FluNet_EIOS_temp_ob$SDATE), max(FluNet_EIOS_temp_ob$SDATE)))
    
    detected_list[[j]] <- sapply(57:73, function(x){
      state <- FluNet_EIOS_temp_ob$epidemic
      alarm <- EIOS_temp[x]
      return(sum(state[alarm == TRUE], na.rm = TRUE) > 0)
    })
  } 
  
  # transform list to df, sum up columns and calculate sens per outbreak
  detected_df <- data.frame(matrix(unlist(detected_list), nrow = max(FluNet_EIOS_temp$outbreak_period), byrow = TRUE), 
                            stringsAsFactors = FALSE)
  EIOS_sens_outbreak_list[[i]] <- apply(detected_df, 2, sum)/max(FluNet_EIOS_temp$outbreak_period)
}

EIOS_sens_outbreak <- data.frame(matrix(unlist(EIOS_sens_outbreak_list), ncol = 17, byrow = TRUE))
names(EIOS_sens_outbreak) <- paste("bcp_", seq(0.1, 0.9, 0.05), sep = "")
EIOS_sens_outbreak$country <- country_list
EIOS_sens_outbreak <- select(EIOS_sens_outbreak, country, everything())


### specificity
