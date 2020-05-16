library(plyr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(bcp)
library(readxl)
library(qcc)
library(surveillance)
library(cpm)
library(tidyr)

setwd("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet")
# load HealthMap data
dataHM1 <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/mcgill flu iris 20190711144158.csv", header = F, 
                    stringsAsFactors = F)
colnames(dataHM1) <- c("place_name", "country", "disease_name", "species_name", "alert_id", "Headline", "Link", 
                       "issue_date", "load_date", "Snippet", "Tag", "Source", "lon", "lat")
dataHM1 <- dataHM1 %>% select(-species_name)
dataHM1$load_date <- as.POSIXct(dataHM1$load_date, format="%Y-%m-%d %H:%M:%S")
dataHM1$issue_date <- as.POSIXct(dataHM1$issue_date, format="%Y-%m-%d %H:%M:%S")
dataHM1 <- dataHM1[!duplicated(dataHM1$alert_id), ]

dataHM2 <- read_excel("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/mcgill flu iris 20191212113549.xlsx")
dataHM2 <- dataHM2 %>% select(-Language)
colnames(dataHM2) <- c("place_name", "country", "disease_name", "alert_id", "Headline", "Link", 
                       "issue_date", "load_date", "Snippet", "Tag", "Source", "lon", "lat")
dataHM2 <- filter(dataHM2, disease_name != "Respiratory Illness")
dataHM2$load_date <- as.POSIXct(dataHM2$load_date, format="%Y-%m-%d %H:%M:%S", tz = "GMT")
dataHM2$issue_date <- as.POSIXct(dataHM2$issue_date, format="%Y-%m-%d %H:%M:%S", tz = "GMT")
dataHM2 <- dataHM2[!duplicated(dataHM2$alert_id), ]


dataHM <- rbind(dataHM1, dataHM2) 

dataHM$country <- as.factor(dataHM$country)
dataHM$Source <- as.factor(dataHM$Source)

overseas_territories <- c("Bermuda [UK]", "CollectivitÃ© d'outre-mer de Saint BarthÃ©lemy, France", "Cayman Islands [UK]", 
                          "Pays d'outre-mer de French Polynesia, France", "RÃ©gion d'outre-mer de Mayotte, France", 
                          "RÃ©gion d'outre-mer de French Guiana, France", "RÃ©gion d'outre-mer de RÃ©union, France", 
                          "RÃ©gion d'outre-mer de Guadeloupe, France", "RÃ©gion d'outre-mer de Martinique, France", 
                          "American Samoa [USA]", "Northern Mariana Islands [United States]", "Guam [USA]")

dataHM <- filter(dataHM, !(country %in% overseas_territories | place_name %in% overseas_territories))
dataHM$country <- factor(dataHM$country)


HM_byweek <- dataHM %>% group_by(date = floor_date(load_date, "week", week_start = 1), country = country, .drop = FALSE) %>%
  summarize(counts=n()) %>% as.data.frame()
HM_byweek <- HM_byweek[order(HM_byweek$country), ]
country_list <- levels(HM_byweek$country)

# smooth counts
HM_byweek$count_smooth <- NA
for (i in seq_along(country_list)) { 
  HM_temp <- filter(HM_byweek, HM_byweek$country==country_list[i])
  count_smooth_temp <- loess(HM_temp$counts ~ as.numeric(HM_temp$date), span = 0.07, degree = 2)
  HM_byweek$count_smooth[HM_byweek$country==country_list[i]] <- count_smooth_temp$fitted
}

# bcp
for (i in seq_along(country_list)) { 
  HM_temp <- filter(HM_byweek, country == country_list[i])
  
  bcp_temp <- bcp(HM_temp$counts, burnin = 100, mcmc = 500)
  HM_byweek$bcp.postprob[HM_byweek$country==country_list[i]] <- bcp_temp$posterior.prob
}

# criteria for start and end of epidemics
HM_byweek$bcp_start <- NA
for (i in seq_along(country_list)) {
  HM_temp <- filter(HM_byweek, country == country_list[i] & is.na(HM_byweek$bcp.postprob) == FALSE)
  
  for(j in 2:13){ # first few rows without looking back for previous outbreak
    if(HM_temp$bcp.postprob[j] >= 0.5 & HM_temp$bcp.postprob[j-1] < 0.5 &
       HM_temp$count_smooth[j] > HM_temp$count_smooth[j-1] &
       sum(HM_temp$counts[(j-1):(j+6)] > 4) != 0){
      HM_temp$bcp_start[j] <- HM_temp$date[j]
    } else {
      HM_temp$bcp_start[j] <- NA
    }
  }
  for(j in 12:(nrow(HM_temp)-6)){
    if(HM_temp$bcp.postprob[j] >= 0.5 & HM_temp$bcp.postprob[(j-1)] < 0.5 & # transition from non-epidemic to epidemic
       HM_temp$count_smooth[j] > HM_temp$count_smooth[j-1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(HM_temp$bcp_start[(j-12):(j-1)])) == 0 & # No outbreak flagged during the previous 10 weeks
       sum(HM_temp$counts[(j-1):(j+6)] > 4) != 0){  
      HM_temp$bcp_start[j] <- HM_temp$date[j]
    } else {
      HM_temp$bcp_start[j] <- NA
    }
  }
  HM_byweek$bcp_start[HM_byweek$country==country_list[i]] <- HM_temp$bcp_start
}


#  apply criteria for end of epidemics
HM_byweek$bcp_end <- NA
for (i in seq_along(country_list)) { 
  HM_temp <- filter(HM_byweek, country == country_list[i] & is.na(HM_byweek$bcp.postprob) == FALSE)
  
  for(j in (nrow(HM_temp)-2):2){# run in reverse because otherwise, third criterion cannot be recognized
    if(HM_temp$bcp.postprob[j-1] >= 0.5 & HM_temp$bcp.postprob[j] < 0.5 & # transition from non-epidemic to epidemic
       HM_temp$count_smooth[j] > HM_temp$count_smooth[j+1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(HM_temp$bcp_end[(j+1):(j+15)])) == 0){ # No outbreak flagged during the previous 15 weeks
      HM_temp$bcp_end[j] <- HM_temp$date[j]
    } else {
      HM_temp$bcp_end[j] <- NA
    }
  }
  HM_byweek$bcp_end[HM_byweek$country==country_list[i]] <- HM_temp$bcp_end
}

# start and end indicators of epidemics
HM_byweek$startend_bcp <- ifelse(is.na(HM_byweek$bcp_start) == FALSE, "start", 
                                 ifelse(is.na(HM_byweek$bcp_end) == FALSE, "end", NA))

#for spikes: if 'start' is not followed by an 'end' within 30 weeks, add an 'end' directly after - this period should be optimized
for(i in 1:nrow(HM_byweek)){
  if(HM_byweek$date[i] < "2019-09-01"){
    if(is.na(HM_byweek$startend_bcp[i]) == FALSE){
      if(HM_byweek$startend_bcp[i] == "start" & !sum(HM_byweek$startend_bcp[(i+1):(i+40)] == "end", na.rm = TRUE) > 0){
        HM_byweek$startend_bcp[i+1] <- "end"
      }
    }
  }
}

# manually insert one 'end' after a spike for Russia
which(HM_byweek$country == "Russia" & HM_byweek$startend_bcp == "start")
HM_byweek$startend_bcp[5359] <- "end"
# manually insert one 'end' after a spike for Mexico
which(HM_byweek$country == "Mexico" & HM_byweek$startend_bcp == "start")
HM_byweek$startend_bcp[4690] <- "end"
# manually remove one 'end' in an epidemic for Brazil
which(HM_byweek$country == "Brazil" & HM_byweek$startend_bcp == "end")
HM_byweek$startend_bcp[860] <- NA
# manually remove one 'end' in an epidemic for USA
which(HM_byweek$country == "United States" & HM_byweek$startend_bcp == "end")
HM_byweek$startend_bcp[7254] <- NA
#manually remove one 'end' in an epidemic for Germany
which(HM_byweek$country == "Germany" & HM_byweek$startend_bcp == "end")
HM_byweek$startend_bcp[3272] <- NA

# repeat date copying into bcp_start and bcp_end columns to correct for inserted 'ends' (and manually curated values)
HM_byweek$bcp_end <- ifelse(HM_byweek$startend_bcp == "end", HM_byweek$date, NA)
HM_byweek$bcp_start <- ifelse(HM_byweek$startend_bcp == "start", HM_byweek$date, NA)


HM_epidemic_sens <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(startend_bcp))) %>% 
  mutate(epidemic = replace(startend_bcp, first(startend_bcp) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp) 
HM_epidemic_sens$epidemic[which(HM_epidemic_sens$epidemic == "end")] <- TRUE
HM_epidemic_sens$epidemic[which(is.na(HM_epidemic_sens$epidemic) == TRUE)] <- FALSE
HM_epidemic_sens$epidemic <- as.logical(HM_epidemic_sens$epidemic)


for (i in seq_along(country_list)) {
  HM_temp <- filter(HM_epidemic_sens, country == country_list[i])
  
  plot <- ggplot(HM_temp, aes(x = date, y = counts, col = epidemic)) + 
    geom_line(aes(group = 1), size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "HealthMap event counts", 
         title = paste("HealthMap data for", country_list[i], "with start and end of epidemics", sep = " ")) + 
    geom_vline(xintercept = na.omit(HM_temp$bcp_start[HM_temp$startend_bcp == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(HM_temp$bcp_end[HM_temp$startend_bcp == "end"]), lty = 2, col = "darkgreen")  + 
    scale_color_manual(values = c("#6e6868", "#e64040"))
  
  print(plot)
  ggsave(plot = plot, file = paste("HealthMap outbreak min 5 counts", country_list[i], ".jpeg", sep=' '))
}

write.csv(HM_epidemic_sens, 
          file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HealthMap_epidemic_min5counts.csv",
          row.names = FALSE)


##### calculate evaluation metrics #####
HM_epidemic_sens$startend <- NA
for(i in 2:nrow(HM_epidemic_sens)){
  if(HM_epidemic_sens$epidemic[i] == TRUE & HM_epidemic_sens$epidemic[i-1] == FALSE){
    HM_epidemic_sens$startend[i] <- "start"
  }else if(HM_epidemic_sens$epidemic[i] == TRUE & HM_epidemic_sens$epidemic[i+1] == FALSE){
    HM_epidemic_sens$startend[i] <- "end"
  }
}

FluNet_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/FluNet_epidemic.csv")
FluNet_epidemic$SDATE <- as.POSIXct(FluNet_epidemic$SDATE)

FluNet_epidemic$startend <- NA
for(i in 2:nrow(FluNet_epidemic)){
  if(FluNet_epidemic$epidemic[i] == TRUE & FluNet_epidemic$epidemic[i-1] == FALSE){
    FluNet_epidemic$startend[i] <- "start"
  }else if(FluNet_epidemic$epidemic[i] == TRUE & FluNet_epidemic$epidemic[i+1] == FALSE){
    FluNet_epidemic$startend[i] <- "end"
  }
}
FluNet_epidemic$startend[which(FluNet_epidemic$SDATE == min(FluNet_epidemic$SDATE))] <- NA
# for Saudi Arabia (data start later)
sa <- which(FluNet_epidemic$Country == "Saudi Arabia")
FluNet_epidemic$startend[sa[1]] <- NA


HM_epidemic_sens <- HM_epidemic_sens %>% group_by(country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
FluNet_epidemic <- FluNet_epidemic %>% group_by(Country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)


# calculate summary information on outbreaks
FluNet_outbreak_length <- FluNet_epidemic %>% filter(epidemic == TRUE) %>% group_by(Country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n(), start_date = min(SDATE), end_date = max(SDATE), height = max(ALL_INF)) %>% 
  mutate(outbreak_no = row_number()) %>% select(-grp)
FluNet_outbreak_length$outbreak_interval <- interval(FluNet_outbreak_length$start_date, FluNet_outbreak_length$end_date)

HM_outbreak_length <- HM_epidemic_sens %>% filter(epidemic == TRUE) %>% group_by(country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n(), start_date = min(date), end_date = max(date), height = max(counts)) %>% 
  mutate(outbreak_no = row_number()) %>% select(-grp)
HM_outbreak_length$outbreak_interval <- interval(HM_outbreak_length$start_date, HM_outbreak_length$end_date)

FluNet_outbreaks_for_HM <- filter(FluNet_outbreak_length, start_date %within% interval(min(HM_epidemic_sens$date), max(HM_epidemic_sens$date))) %>% 
  group_by(Country) %>% summarize(no_FluNet_for_HM = n())


# sens per ourbtreak
outbreaks_detected <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_HM_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                                start_date %within% interval(min(HM_epidemic_sens$date), max(HM_epidemic_sens$date)))
  HM_temp_ol <- filter(HM_outbreak_length, country == country_list[i])
  
  sum_HM <- NA
  for(j in 1:nrow(HM_temp_ol)){
    sum_HM[j] <- sum(int_overlaps(HM_temp_ol$outbreak_interval[j], FluNet_HM_temp_ol$outbreak_interval))
    no_detected_outbreaks_HM <- sum(sum_HM)
  }

  outbreaks_detected$country[i] <- country_list[i]
  outbreaks_detected$no_HM[i] <- no_detected_outbreaks_HM
}
sens_per_outbreak <- data.frame(outbreaks_detected$country, FluNet_outbreaks_for_HM$no_FluNet_for_HM, outbreaks_detected$no_HM)
names(sens_per_outbreak) <- c("country", "no_FluNet_for_HM", "no_HM")
sens_per_outbreak <- sens_per_outbreak %>% mutate(sens_per_outbreak = no_HM/no_FluNet_for_HM)

# sens per week
sens_per_week_list <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic_sens$date), max(HM_epidemic_sens$date)))
  HM_temp <- filter(HM_epidemic_sens, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  sens_HM <- sum(alarm_HM[state_HM == TRUE], na.rm = TRUE)
  sens_HM <- sens_HM / length(alarm_HM[state_HM == TRUE])
 
  sens_per_week_list$country[i] <- country_list[i]
  sens_per_week_list$sens_HM[i] <- sens_HM
}

sens_per_week <- data.frame(sens_per_week_list$country, sens_per_week_list$sens_HM)
names(sens_per_week) <- c("country", "sens_per_week")

# timely sens
detected_exact <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_HM_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                                start_date %within% interval(min(HM_epidemic_sens$date), max(HM_epidemic_sens$date)))
  HM_temp_ol <- filter(HM_outbreak_length, country == country_list[i])
  
  sum_HM <- NA
  for(j in 1:nrow(HM_temp_ol)){
    sum_HM[j] <- sum(HM_temp_ol$start_date[j] %within%
                       interval(FluNet_HM_temp_ol$start_date[j] - weeks(2), FluNet_HM_temp_ol$start_date[j] + weeks(2)),
                     na.rm = TRUE)
    no_detected_outbreaks_HM <- sum(sum_HM)
  }
  
  detected_exact$country[i] <- country_list[i]
  detected_exact$no_HM[i] <- no_detected_outbreaks_HM
}

sens_exact <- data.frame(detected_exact$country, FluNet_outbreaks_for_HM$no_FluNet_for_HM, detected_exact$no_HM)
names(sens_exact) <- c("country", "no_FluNet_for_HM", "no_HM")
sens_exact <- sens_exact %>% mutate(sens_exact = no_HM/no_FluNet_for_HM)


# PPV
PPV_list <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic_sens$date), max(HM_epidemic_sens$date)))
  HM_temp <- filter(HM_epidemic_sens, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  TP_HM <- sum(alarm_HM[state_HM == TRUE], na.rm = TRUE)
  PPV_HM <- TP_HM / sum(alarm_HM, na.rm = TRUE)
   
  PPV_list$country[i] <- country_list[i]
  PPV_list$PPV_HM[i] <- PPV_HM
}

PPV <- data.frame(PPV_list$country, PPV_list$PPV_HM)
names(PPV) <- c("country", "PPV")


# false alarm rate
false_alarm_rate <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic_sens$date), max(HM_epidemic_sens$date)))
  HM_temp <- filter(HM_epidemic_sens, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  FA_HM <- sum(alarm_HM[state_HM == FALSE], na.rm = TRUE)
  FAR_HM <- FA_HM / length(alarm_HM[state_HM == FALSE])
  
  false_alarm_rate$country[i] <- country_list[i]
  false_alarm_rate$FAR_HM[i] <- FAR_HM
}

specificity <- data.frame(false_alarm_rate$country, 1- false_alarm_rate$FAR_HM)
names(specificity) <- c("country", "specificity")


# prevented fraction
timeliness_HM <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic_sens$date), max(HM_epidemic_sens$date)))
  prevented <- 0
  for(j in 1:max(FluNet_HM_temp$outbreak_period)){
    FluNet_HM_temp_ob <- filter(FluNet_HM_temp, outbreak_period == j)
    HM_temp <- filter(HM_epidemic_sens, country == country_list[i] & 
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

timeliness <- data.frame(timeliness_HM$country, timeliness_HM$timeliness)
names(timeliness) <- c("country", "frac_prevented")


# accuracy
accuracy_list <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic_sens$date), max(HM_epidemic_sens$date)))
  HM_temp <- filter(HM_epidemic_sens, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  TP_HM <- sum(alarm_HM[state_HM == TRUE], na.rm = TRUE)
  TN_HM <- sum(alarm_HM[state_HM == FALSE] == FALSE, na.rm = TRUE)
  accuracy_HM <- (TP_HM + TN_HM)/nrow(HM_temp)
   
  accuracy_list$country[i] <- country_list[i]
  accuracy_list$accuracy_HM[i] <- accuracy_HM
}

accuracy <- data.frame(accuracy_list$country, accuracy_list$accuracy_HM)
names(accuracy) <- c("country", "accuracy")


# combine all metrics into one df
metrics_sens <- data.frame(sens_per_outbreak$country, sens_per_outbreak$sens_per_outbreak, sens_exact$sens_exact, 
                      sens_per_week$sens_per_week, PPV$PPV, specificity$specificity, timeliness$frac_prevented)
names(metrics_sens) <- c("country", "sens_per_outbreak", "sens_exact", "sens_per_week", "PPV", "specificity", "frac_prevented")

write.csv(metrics_sens, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics_min5counts.csv")



##### compare with previous results for HM #####
metrics <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics.csv", stringsAsFactors = TRUE)
metrics_HM <- metrics %>% filter(source == "HealthMap") %>% select(-c(X, source))

comp_df <- rbind(metrics_HM, metrics_sens)
comp_df$count_limit <- rep(c("no limit", "min 5"), each = 24)

mean_metrics <- comp_df %>% group_by(count_limit) %>% summarize_if(is.numeric, mean)
median_metrics <- comp_df %>% group_by(count_limit) %>% summarize_if(is.numeric, median)

mean_metrics_long <- pivot_longer(mean_metrics, cols = 2:7, names_to = "metric", values_to = "value")
median_metrics_long <- pivot_longer(median_metrics, cols = 2:7, names_to = "metric", values_to = "value")


ggplot(mean_metrics_long, aes(y = value, x = metric, fill = count_limit)) + 
  geom_col(position = "dodge") + 
  labs(title = "Mean evaluation metrics", y ="", x = "", fill = "count limit") + 
  scale_x_discrete(labels = c("prevented fraction", "positive \npredictive value", "exact sensitivity", "sensitivity \nper outbreak",
                              "sensitivity \nper week", "specificity")) + 
  scale_fill_brewer(palette = "Paired") + 
  scale_y_continuous(limits = c(0, 0.99))
ggsave("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/limit_comp_mean.jpeg")

ggplot(median_metrics_long, aes(y = value, x = metric, fill = count_limit)) + 
  geom_col(position = "dodge") + 
  labs(title = "Median evaluation metrics", y ="", x = "", fill = "count limit") + 
  scale_x_discrete(labels = c("prevented fraction", "positive \npredictive value", "exact sensitivity", "sensitivity \nper outbreak",
                              "sensitivity \nper week", "specificity")) + 
  scale_fill_brewer(palette = "Paired") + 
  scale_y_continuous(limits = c(0, 0.99))
ggsave("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/limit_comp_median.jpeg")


comp_df_long <- pivot_longer(comp_df, cols = 2:7, names_to = "metric", values_to = "value")
metrics.labs <- c("sensitivity per outbreak", "timely sensitivity", "sensitivity per week", "positive predictive value", 
                   "specificity", "prevented fraction")
names(metrics.labs) <- names(comp_df)[2:7]

ggplot(comp_df_long, aes(x = country, y = value, fill = count_limit)) + 
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Comparison of count limit for outbreak detection in HealthMap", y = "", x = "", fill = "count limit") +
  scale_x_discrete(limits = rev(levels(comp_df_long$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Paired") + 
  facet_wrap(facets = ~metric, ncol = 2, labeller = labeller(metric = metrics.labs))
ggsave("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/limit_comp_all_metrics.jpeg", 
       width = 13, height = 10, unit = "in")
