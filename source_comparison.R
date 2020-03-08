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

# revise startend columns of datasets
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

FluNet_epidemic$startend <- NA
for(i in 2:nrow(FluNet_epidemic)){
  if(FluNet_epidemic$epidemic[i] == TRUE & FluNet_epidemic$epidemic[i-1] == FALSE){
    FluNet_epidemic$startend[i] <- "start"
  }else if(FluNet_epidemic$epidemic[i] == FALSE & FluNet_epidemic$epidemic[i-1] == TRUE){
    FluNet_epidemic$startend[i] <- "end"
  }
}

# outbreak periods (from one start to next start)
HM_epidemic <- HM_epidemic %>% group_by(country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
EIOS_epidemic <- EIOS_epidemic %>% group_by(country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
FluNet_epidemic <- FluNet_epidemic %>% group_by(Country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)



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
  FluNet_EIOS_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                               start_date %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                             start_date %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp_ol <- filter(HM_outbreak_length, country == country_list[i])
  EIOS_temp_ol <- filter(EIOS_outbreak_length, country == country_list[i])
 
  sum_HM <- NA
  for(j in 1:nrow(HM_temp_ol)){
    sum_HM[j] <- sum(int_overlaps(HM_temp_ol$outbreak_interval[j], FluNet_HM_temp_ol$outbreak_interval))
    no_detected_outbreaks_HM <- sum(sum_HM)
  }
  sum_EIOS <- NA
  for(k in 1:nrow(EIOS_temp_ol)){
    sum_EIOS[k] <- sum(int_overlaps(EIOS_temp_ol$outbreak_interval[k], FluNet_EIOS_temp_ol$outbreak_interval))
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


sens_per_outbreak <- data.frame(outbreaks_detected$country, FluNet_outbreaks_for_HM$no_FluNet_for_HM, outbreaks_detected$no_HM,
                           FluNet_outbreaks_for_EIOS$no_FluNet_for_EIOS, outbreaks_detected$no_EIOS)
names(sens_per_outbreak) <- c("country", "no_FluNet_for_HM", "no_HM", "no_FluNet_for_EIOS", "no_EIOS")

sens_per_outbreak <- sens_per_outbreak %>% mutate(sens_HM = no_HM/no_FluNet_for_HM, sens_EIOS = no_EIOS/no_FluNet_for_EIOS)

#write.csv(no_outbreaks, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/sensitivity.csv")


# prepare data for plotting and plot
sens_per_outbreak_long <- pivot_longer(sens_per_outbreak, cols = c(sens_HM, sens_EIOS), names_to = "source", values_to = "sensitivity")
sens_per_outbreak_long$source <- factor(sens_per_outbreak_long$source, levels = c("sens_HM", "sens_EIOS"), labels = c("HealthMap", "EIOS"))
sens_per_outbreak_long$sensitivity <- ifelse(sens_per_outbreak_long$sensitivity > 1, 1, sens_per_outbreak_long$sensitivity)

ggplot(sens_per_outbreak_long, aes(x = country, y = sensitivity * 100, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Sensitivity per outbreak of HealthMap and EIOS systems", y = "Sensitivity per outbreak (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(sens_per_outbreak_long$country))) +
  scale_fill_brewer(palette = "Set1")



##### calculate sensitivity per week #####
sens_per_week_list <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  sens_HM <- sum(alarm_HM[state_HM == TRUE], na.rm = TRUE)
  sens_HM <- sens_HM / length(alarm_HM[state_HM == TRUE])
  
  state_EIOS <- FluNet_EIOS_temp$epidemic
  alarm_EIOS <- EIOS_temp$epidemic
  sens_EIOS <- sum(alarm_EIOS[state_EIOS == TRUE], na.rm = TRUE)
  sens_EIOS <- sens_EIOS / length(alarm_EIOS[state_EIOS == TRUE])
  
  sens_per_week_list$country[i] <- country_list[i]
  sens_per_week_list$sens_HM[i] <- sens_HM
  sens_per_week_list$sens_EIOS[i] <- sens_EIOS
}

sens_per_week <- data.frame(sens_per_week_list$country, sens_per_week_list$sens_HM, sens_per_week_list$sens_EIOS)
names(sens_per_week) <- c("country", "sens_HM", "sens_EIOS")

sens_per_week_long <- pivot_longer(sens_per_week, cols = c(sens_HM, sens_EIOS), names_to = "source", 
                                   values_to = "sensitivity_per_week")
sens_per_week_long$source <- factor(sens_per_week_long$source, levels = c("sens_HM", "sens_EIOS"), labels = c("HealthMap", "EIOS"))

ggplot(sens_per_week_long, aes(x = country, y = sensitivity_per_week * 100, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Sensitivity per week of HealthMap and EIOS systems", y = "Sensitivity per week (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(sens_per_week_long$country))) +
  scale_y_continuous(limits = c(0, 100)) + 
  scale_fill_brewer(palette = "Set1")



##### calculate false alarm rate #####
# calculate false alert rate during non-outbreak intervals

false_alarm_rate <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  FA_HM <- sum(alarm_HM[state_HM == FALSE], na.rm = TRUE)
  FAR_HM <- FA_HM / length(alarm_HM[state_HM == FALSE])
  
  state_EIOS <- FluNet_EIOS_temp$epidemic
  alarm_EIOS <- EIOS_temp$epidemic
  FA_EIOS <- sum(alarm_EIOS[state_EIOS == FALSE], na.rm = TRUE)
  FAR_EIOS <- FA_EIOS / length(alarm_EIOS[state_EIOS == FALSE])
  
  false_alarm_rate$country[i] <- country_list[i]
  false_alarm_rate$FAR_HM[i] <- FAR_HM
  false_alarm_rate$FAR_EIOS[i] <- FAR_EIOS
}

FAR <- data.frame(false_alarm_rate$country, false_alarm_rate$FAR_HM, false_alarm_rate$FAR_EIOS)
names(FAR) <- c("country", "FAR_HM", "FAR_EIOS")

FAR_long <- pivot_longer(FAR, cols = c(FAR_HM, FAR_EIOS), names_to = "source", values_to = "false_alarm_rate")
FAR_long$source <- factor(FAR_long$source, levels = c("FAR_HM", "FAR_EIOS"), labels = c("HealthMap", "EIOS"))
FAR_long$specificity <- 1 - FAR_long$false_alarm_rate

ggplot(FAR_long, aes(x = country, y = false_alarm_rate, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "False alarm rate of HealthMap and EIOS systems", y = "False alarm rate", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(FAR_long$country))) +
  scale_fill_brewer(palette = "Set1")

ggplot(FAR_long, aes(x = country, y = specificity * 100, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Specificity of HealthMap and EIOS systems", y = "Specificity (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(FAR_long$country))) +
  scale_fill_brewer(palette = "Set1")


##### calculate timeliness#####
#  input df is data frame with 3 col: date, baseline, baseline + outbreak
#  input truth is list with element: state, which is binary vector indicating presence of outbreak
#  output is proportion of outbreak "prevented" (from 1 for first day to 0 for last day of outbreak, or later)


timeliness_HM <- vector(mode = "list")
timeliness_EIOS <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  prevented <- 0
  for(j in 1:max(FluNet_HM_temp$outbreak_period)){
    FluNet_HM_temp_ob <- filter(FluNet_HM_temp, outbreak_period == j)
    HM_temp <- filter(HM_epidemic, country == country_list[i] & 
                        date %within% interval(min(FluNet_HM_temp_ob$SDATE), max(FluNet_HM_temp_ob$SDATE)))
    
    state <- FluNet_HM_temp_ob$epidemic
    alarm <- HM_temp$epidemic
    length <- sum(FluNet_HM_temp_ob$epidemic)
    
    detect <- ifelse(sum(alarm[state == TRUE], na.rm = TRUE) > 0, TRUE, FALSE)
    if (detect) {
      first.alarm <- min(which(alarm[state == TRUE] == TRUE))
      prevented[j] <- (length - first.alarm) / length
    }else{
      prevented[j] <- 0
    }
  }
  timeliness_HM$country[i] <- country_list[i]
  timeliness_HM$timeliness[i] <- mean(prevented, na.rm = TRUE)
}  

for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_EIOS_temp <- FluNet_EIOS_temp %>% group_by(Country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
  
  prevented <- 0
  for(j in 1:max(FluNet_EIOS_temp$outbreak_period)){
    FluNet_EIOS_temp_ob <- filter(FluNet_EIOS_temp, outbreak_period == j)
    EIOS_temp <- filter(EIOS_epidemic, country == country_list[i] & 
                        date %within% interval(min(FluNet_EIOS_temp_ob$SDATE), max(FluNet_EIOS_temp_ob$SDATE)))
    
    state <- FluNet_EIOS_temp_ob$epidemic
    alarm <- EIOS_temp$epidemic
    length <- sum(FluNet_EIOS_temp_ob$epidemic)
    
    detect <- ifelse(sum(alarm[state == TRUE], na.rm = TRUE) > 0, TRUE, FALSE)
    if (detect == TRUE){
      first.alarm <- min(which(alarm[state == TRUE] == TRUE))
      prevented[j] <- (length - first.alarm) / length
    }else{
      prevented[j] <- 0
    }
  }
  timeliness_EIOS$country[i] <- country_list[i]
  timeliness_EIOS$timeliness[i] <- mean(prevented, na.rm = TRUE)
} 

timeliness <- data.frame(timeliness_HM$country, timeliness_HM$timeliness, timeliness_EIOS$timeliness)
names(timeliness) <- c("country", "frac_prevented_HM", "frac_prevented_EIOS")

timeliness_long <- pivot_longer(timeliness, cols = c(frac_prevented_HM, frac_prevented_EIOS), 
                                names_to = "source", values_to = "frac_prevented")
timeliness_long$source <- factor(timeliness_long$source, levels = c("frac_prevented_HM", "frac_prevented_EIOS"), 
                                 labels = c("HealthMap", "EIOS"))


ggplot(timeliness_long, aes(x = country, y = frac_prevented, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Timeliness (prevented fraction) of HealthMap and EIOS systems", y = "prevented fraction", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(timeliness_long$country))) +
  scale_fill_brewer(palette = "Set1")


##### combine all metrics into one dataframe ##### 
metrics <- data.frame(sens_per_outbreak_long$country, sens_per_outbreak_long$source, sens_per_outbreak_long$sensitivity, 
                       sens_per_week_long$sensitivity_per_week, FAR_long$specificity, timeliness_long$frac_prevented)
names(metrics) <- c("country", "source", "sens_per_outbreak", "sens_per_week", "specificity", "frac_prevented")

metrics_long <- pivot_longer(metrics, cols = c(sens_per_outbreak, sens_per_week, specificity, frac_prevented), 
                             names_to = "metric", values_to = "values")


for(i in seq_along(country_list)){
  metrics_temp <- filter(metrics_long, country == country_list[i])

  plot <- ggplot(metrics_temp, aes(x = metric, y = values, fill = source)) +
    geom_col(position = "dodge") +
    labs(title = paste(country_list[i], ": Evaluation of HealthMap and EIOS systems", sep = ""), y = "", x = "", 
         caption = "WHO FluNet data were used as gold standard") +
    scale_x_discrete(labels = c('Prevented fraction \n (Timeliness)','Sensitivity per outbreak', 'Sensitivity per week', 'Specificity')) +
    scale_y_continuous(limits = c(0, 1)) + 
    scale_fill_brewer(palette = "Set1") 
  print(plot)
  ggsave(plot = plot, file = paste("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/Evaluation metrics", country_list[i], ".jpeg", sep=' '))
}
