library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# load epidemic datasets with outbreak indicators
FluNet_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/FluNet_epidemic.csv")
FluNet_epidemic$SDATE <- as.POSIXct(FluNet_epidemic$SDATE)
HM_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HealthMap_epidemic.csv")
HM_epidemic$date <- as.POSIXct(HM_epidemic$date)
EIOS_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_epidemic.csv")
EIOS_epidemic$date <- as.POSIXct(EIOS_epidemic$date)

country_list <- levels(FluNet_epidemic$Country)

# revise startend columns of HM and EIOS datasets
HM_epidemic$startend <- NA
for(i in 2:nrow(HM_epidemic)){
  if(HM_epidemic$epidemic[i] == TRUE & HM_epidemic$epidemic[i-1] == FALSE){
    HM_epidemic$startend[i] <- "start"
  }else if(HM_epidemic$epidemic[i] == FALSE & HM_epidemic$epidemic[i-1] == TRUE){
    HM_epidemic$startend[i] <- "end"
  }
}

EIOS_epidemic$startend <- NA
for(i in 2:nrow(EIOS_epidemic)){
  if(EIOS_epidemic$epidemic[i] == TRUE & EIOS_epidemic$epidemic[i-1] == FALSE){
    EIOS_epidemic$startend[i] <- "start"
  }else if(EIOS_epidemic$epidemic[i] == FALSE & EIOS_epidemic$epidemic[i-1] == TRUE){
    EIOS_epidemic$startend[i] <- "end"
  }
}

# calculate summary information on outbreaks
FluNet_outbreak_length <- FluNet_epidemic %>% filter(epidemic == TRUE) %>% group_by(Country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n(), start_date = min(SDATE), end_date = max(SDATE), height = max(ALL_INF)) %>% 
  mutate(outbreak_no = row_number()) %>% select(-grp)
FluNet_outbreak_length$outbreak_interval <- interval(FluNet_outbreak_length$start_date, FluNet_outbreak_length$end_date)

HM_outbreak_length <- HM_epidemic %>% filter(epidemic == TRUE) %>% group_by(country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n(), start_date = min(date), end_date = max(date), height = max(counts)) %>% 
  mutate(outbreak_no = row_number()) %>% select(-grp)
HM_outbreak_length$outbreak_interval <- interval(HM_outbreak_length$start_date, HM_outbreak_length$end_date)

EIOS_outbreak_length <- EIOS_epidemic %>% filter(epidemic == TRUE) %>% group_by(country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n(), start_date = min(date), end_date = max(date), height = max(counts)) %>% 
  mutate(outbreak_no = row_number()) %>% select(-grp)
EIOS_outbreak_length$outbreak_interval <- interval(EIOS_outbreak_length$start_date, EIOS_outbreak_length$end_date)



##### calculate sensitivity per outbreak for HM and EIOS #####
# input truth: outbreak intervals from FluNet data for every country
# input source: start and end dates from EIOS or HM data for every country

outbreaks_detected <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                               start_date %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                             start_date %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_outbreak_length, country == country_list[i])
  EIOS_temp <- filter(EIOS_outbreak_length, country == country_list[i])
 
  sum_HM <- NA
  for(j in 1:nrow(HM_temp)){
    sum_HM[j] <- sum(int_overlaps(HM_temp$outbreak_interval[j], FluNet_HM_temp$outbreak_interval))
    no_detected_outbreaks_HM <- sum(sum_HM)
  }
  sum_EIOS <- NA
  for(k in 1:nrow(EIOS_temp)){
    sum_EIOS[k] <- sum(int_overlaps(EIOS_temp$outbreak_interval[k], FluNet_EIOS_temp$outbreak_interval))
    no_detected_outbreaks_EIOS <- sum(sum_EIOS)
  } 
  
  outbreaks_detected$country[i] <- country_list[i]
  outbreaks_detected$no_HM[i] <- no_detected_outbreaks_HM
  outbreaks_detected$no_EIOS[i] <- no_detected_outbreaks_EIOS
}

FluNet_outbreaks_for_HM <- filter(FluNet_outbreak_length, start_date %within% interval(min(HM_epidemic$date), max(HM_epidemic$date))) %>% 
  group_by(Country) %>% summarize(no_FluNet_for_HM = n())

FluNet_outbreaks_for_EIOS <- filter(FluNet_outbreak_length, start_date %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date))) %>% 
  group_by(Country) %>% summarize(no_FluNet_for_EIOS = n())


no_outbreaks <- data.frame(outbreaks_detected$country, FluNet_outbreaks_for_HM$no_FluNet_for_HM, outbreaks_detected$no_HM,
                           FluNet_outbreaks_for_EIOS$no_FluNet_for_EIOS, outbreaks_detected$no_EIOS)
names(no_outbreaks) <- c("country", "no_FluNet_for_HM", "no_HM", "no_FluNet_for_EIOS", "no_EIOS")

no_outbreaks <- no_outbreaks %>% mutate(sens_HM = no_HM/no_FluNet_for_HM*100, sens_EIOS = no_EIOS/no_FluNet_for_EIOS*100)

#write.csv(no_outbreaks, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/sensitivity.csv")


# prepare data for plotting and plot
source_sensitivity <- pivot_longer(no_outbreaks, cols = c(sens_HM, sens_EIOS), names_to = "source", values_to = "sensitivity")
source_sensitivity$source <- factor(source_sensitivity$source, levels = c("sens_HM", "sens_EIOS"), labels = c("HealthMap", "EIOS"))
source_sensitivity$sensitivity <- ifelse(source_sensitivity$sensitivity > 100, 100, source_sensitivity$sensitivity)

ggplot(source_sensitivity, aes(x = country, y = sensitivity, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Sensitivity per outbreak of HealthMap and EIOS systems", y = "sensitivity (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(source_sensitivity$country))) +
  scale_fill_brewer(palette = "Set1")
  

##### calculate false alarm rate #####