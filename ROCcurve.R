library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(ROCR)

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
      alarm <- HM_temp[, x]
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
      alarm <- EIOS_temp[, x]
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


### sensitivity per week ####
HM_sens_week_list <- vector(mode = "list")
EIOS_sens_week_list <- vector(mode = "list")

for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  HM_sens_list <- sapply(57:73, function(x){
    state_HM <- FluNet_HM_temp$epidemic
    alarm_HM <- HM_temp[, x]
    sens_HM <- sum(alarm_HM[state_HM == TRUE], na.rm = TRUE) / length(alarm_HM[state_HM == TRUE])
    return(sens_HM)
  })
  
  EIOS_sens_list <- sapply(57:73, function(x){
    state_EIOS <- FluNet_EIOS_temp$epidemic
    alarm_EIOS <- EIOS_temp[, x]
    sens_EIOS <- sum(alarm_EIOS[state_EIOS == TRUE], na.rm = TRUE) / length(alarm_EIOS[state_EIOS == TRUE])
    return(sens_EIOS)
  })
  
  HM_sens_week_list[[i]] <- HM_sens_list
  EIOS_sens_week_list[[i]] <- EIOS_sens_list
}

HM_sens_week <- data.frame(matrix(unlist(HM_sens_week_list), ncol = 17, byrow = TRUE))
names(HM_sens_week) <- paste("bcp_", seq(0.1, 0.9, 0.05), sep = "")
HM_sens_week$country <- country_list
HM_sens_week <- select(HM_sens_week, country, everything())

EIOS_sens_week <- data.frame(matrix(unlist(EIOS_sens_week_list), ncol = 17, byrow = TRUE))
names(EIOS_sens_week) <- paste("bcp_", seq(0.1, 0.9, 0.05), sep = "")
EIOS_sens_week$country <- country_list
EIOS_sens_week <- select(EIOS_sens_week, country, everything())



### false alarm rate
HM_far_list <- vector(mode = "list")
EIOS_far_list <- vector(mode = "list")

for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  HM_far <- sapply(57:73, function(x){
    state_HM <- FluNet_HM_temp$epidemic
    alarm_HM <- HM_temp[, x]
    far_HM <- sum(alarm_HM[state_HM == FALSE], na.rm = TRUE) / length(alarm_HM[state_HM == FALSE])
    return(far_HM)
  })
  
  EIOS_far <- sapply(57:73, function(x){
    state_EIOS <- FluNet_EIOS_temp$epidemic
    alarm_EIOS <- EIOS_temp[, x]
    far_EIOS <- sum(alarm_EIOS[state_EIOS == FALSE], na.rm = TRUE) / length(alarm_EIOS[state_EIOS == FALSE])
    return(far_EIOS)
  })
  
  HM_far_list[[i]] <- HM_far
  EIOS_far_list[[i]] <- EIOS_far
}

HM_false_alarm_rate <- data.frame(matrix(unlist(HM_far_list), ncol = 17, byrow = TRUE))
names(HM_false_alarm_rate) <- paste("bcp_", seq(0.1, 0.9, 0.05), sep = "")
HM_false_alarm_rate$country <- country_list
HM_false_alarm_rate <- select(HM_false_alarm_rate, country, everything())

EIOS_false_alarm_rate <- data.frame(matrix(unlist(EIOS_far_list), ncol = 17, byrow = TRUE))
names(EIOS_false_alarm_rate) <- paste("bcp_", seq(0.1, 0.9, 0.05), sep = "")
EIOS_false_alarm_rate$country <- country_list
EIOS_false_alarm_rate <- select(EIOS_false_alarm_rate, country, everything())


### PPV ####
HM_PPV_list <- vector(mode = "list")
EIOS_PPV_list <- vector(mode = "list")

for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  HM_ppv <- sapply(57:73, function(x){
    state_HM <- FluNet_HM_temp$epidemic
    alarm_HM <- HM_temp[, x]
    ppv_HM <- sum(alarm_HM[state_HM == FALSE], na.rm = TRUE) / sum(alarm_HM, na.rm = TRUE)
    return(ppv_HM)
  })
  
  EIOS_ppv <- sapply(57:73, function(x){
    state_EIOS <- FluNet_EIOS_temp$epidemic
    alarm_EIOS <- EIOS_temp[, x]
    ppv_EIOS <- sum(alarm_EIOS[state_EIOS == FALSE], na.rm = TRUE) / sum(alarm_EIOS, na.rm = TRUE)
    return(ppv_EIOS)
  })
  
  HM_PPV_list[[i]] <- HM_ppv
  EIOS_PPV_list[[i]] <- EIOS_ppv
}

HM_PPV <- data.frame(matrix(unlist(HM_PPV_list), ncol = 17, byrow = TRUE))
names(HM_PPV) <- paste("bcp_", seq(0.1, 0.9, 0.05), sep = "")
HM_PPV$country <- country_list
HM_PPV <- select(HM_PPV, country, everything())

EIOS_PPV <- data.frame(matrix(unlist(EIOS_PPV_list), ncol = 17, byrow = TRUE))
names(EIOS_PPV) <- paste("bcp_", seq(0.1, 0.9, 0.05), sep = "")
EIOS_PPV$country <- country_list
EIOS_PPV <- select(EIOS_PPV, country, everything())



### timeliness #### 
HM_timeliness_list <- vector(mode = "list")
for(i in seq_along(country_list)){
  prevented_list <- vector(mode = "list")
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  FluNet_HM_temp <- FluNet_HM_temp %>% group_by(Country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
  
    
  for(j in 1:max(FluNet_HM_temp$outbreak_period)){
    FluNet_HM_temp_ob <- filter(FluNet_HM_temp, outbreak_period == j)
    HM_temp <- filter(HM_epidemic, country == country_list[i] & 
                        date %within% interval(min(FluNet_HM_temp_ob$SDATE), max(FluNet_HM_temp_ob$SDATE)))
    prevented <- 0
    
    prevented_list[[j]] <- sapply(57:73, function(x){
      state <- FluNet_HM_temp_ob$epidemic
      alarm <- HM_temp[, x]
      length <- sum(FluNet_HM_temp_ob$epidemic)
      
      detect <- ifelse(sum(alarm[state == TRUE], na.rm = TRUE) > 0, TRUE, FALSE)
      if (detect) {
        first.alarm <- min(which(alarm[state == TRUE] == TRUE))
        prevented <- (length - first.alarm) / length
      }else{
        prevented <- 0
      }
      return(prevented)
    })

  }
  
  timeliness_df <- data.frame(matrix(unlist(prevented_list), nrow = max(FluNet_HM_temp$outbreak_period), byrow = TRUE), 
                            stringsAsFactors = FALSE)
  HM_timeliness_list[[i]] <- apply(timeliness_df, 2, mean)
}  

HM_timeliness <- data.frame(matrix(unlist(HM_timeliness_list), ncol = 17, byrow = TRUE))
names(HM_timeliness) <- paste("bcp_", seq(0.1, 0.9, 0.05), sep = "")
HM_timeliness$country <- country_list
HM_timeliness <- select(HM_timeliness, country, everything())


EIOS_timeliness_list <- vector(mode = "list")
for(i in seq_along(country_list)){
  prevented_list <- vector(mode = "list")
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_EIOS_temp <- FluNet_EIOS_temp %>% group_by(Country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
  
  
  for(j in 1:max(FluNet_EIOS_temp$outbreak_period)){
    FluNet_EIOS_temp_ob <- filter(FluNet_EIOS_temp, outbreak_period == j)
    EIOS_temp <- filter(EIOS_epidemic, country == country_list[i] & 
                        date %within% interval(min(FluNet_EIOS_temp_ob$SDATE), max(FluNet_EIOS_temp_ob$SDATE)))
    prevented <- 0
    
    prevented_list[[j]] <- sapply(57:73, function(x){
      state <- FluNet_EIOS_temp_ob$epidemic
      alarm <- EIOS_temp[, x]
      length <- sum(FluNet_EIOS_temp_ob$epidemic)
      
      detect <- ifelse(sum(alarm[state == TRUE], na.rm = TRUE) > 0, TRUE, FALSE)
      if (detect) {
        first.alarm <- min(which(alarm[state == TRUE] == TRUE))
        prevented <- (length - first.alarm) / length
      }else{
        prevented <- 0
      }
      return(prevented)
    })
    
  }
  
  timeliness_df <- data.frame(matrix(unlist(prevented_list), nrow = max(FluNet_EIOS_temp$outbreak_period), byrow = TRUE), 
                              stringsAsFactors = FALSE)
  EIOS_timeliness_list[[i]] <- apply(timeliness_df, 2, mean)
}  

EIOS_timeliness <- data.frame(matrix(unlist(EIOS_timeliness_list), ncol = 17, byrow = TRUE))
names(EIOS_timeliness) <- paste("bcp_", seq(0.1, 0.9, 0.05), sep = "")
EIOS_timeliness$country <- country_list
EIOS_timeliness <- select(EIOS_timeliness, country, everything())


#### combine all metrics into one df ####
metrics_df_HM <- rbind(HM_false_alarm_rate, HM_PPV, HM_sens_outbreak, HM_sens_week, HM_timeliness)
metrics_df_HM$metric <- rep(c("FAR", "PPV", "sens_outbreak", "sens_week", "frac_prevented"), each = 24)
metrics_df_HM_long <- pivot_longer(metrics_df_HM, cols = 2:18, names_to = "cutoff", values_to = "value")
metrics_df_HM_wide <- pivot_wider(metrics_df_HM_long, names_from = "metric", values_from = "value")
metrics_df_HM_wide <- metrics_df_HM_wide %>% mutate(cutoff = gsub("bcp_", "", cutoff))

metrics_df_EIOS <- rbind(EIOS_false_alarm_rate, EIOS_PPV, EIOS_sens_outbreak, EIOS_sens_week, EIOS_timeliness)
metrics_df_EIOS$metric <- rep(c("FAR", "PPV", "sens_outbreak", "sens_week", "frac_prevented"), each = 24)
metrics_df_EIOS_long <- pivot_longer(metrics_df_EIOS, cols = 2:18, names_to = "cutoff", values_to = "value")
metrics_df_EIOS_wide <- pivot_wider(metrics_df_EIOS_long, names_from = "metric", values_from = "value")
metrics_df_EIOS_wide <- metrics_df_EIOS_wide %>% mutate(cutoff = gsub("bcp_", "", cutoff))


for(i in seq_along(country_list)){
  plot <- ggplot(filter(metrics_df_EIOS_wide, country == country_list[i]), aes(x = FAR, y = sens_week)) + 
    geom_point(size = 3, shape = 19) + 
    labs(title = paste("EIOS ROC curve for", country_list[i]), x = "false alarm rate", y = "sensitivity per week") + 
    scale_y_continuous(limits = c(-0.01, 1)) + 
    scale_x_continuous(limits = c(-0.01, 1)) + 
    geom_text_repel(aes(label = cutoff))
  print(plot)
  ggsave(plot = plot, file = paste("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS ROC", 
                                   country_list[i], ".jpeg", sep=' '))
}

for(i in seq_along(country_list)){
  plot <- ggplot(filter(metrics_df_HM_wide, country == country_list[i]), aes(x = FAR, y = sens_week)) + 
    geom_point(size = 3, shape = 19) + 
    labs(title = paste("HealthMap ROC curve for", country_list[i]), x = "false alarm rate", y = "sensitivity per week") + 
    scale_y_continuous(limits = c(-0.01, 1)) + 
    scale_x_continuous(limits = c(-0.01, 1)) + 
    geom_text_repel(aes(label = cutoff))
  print(plot)
  ggsave(plot = plot, file = paste("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HM ROC", 
                                   country_list[i], ".jpeg", sep=' '))
}

# AMOC 
for(i in seq_along(country_list)){
  plot <- ggplot(filter(metrics_df_EIOS_wide, country == country_list[i]), aes(x = FAR, y = frac_prevented)) + 
    geom_point(size = 3, shape = 19) + 
    labs(title = paste("EIOS AMOC curve for", country_list[i]), x = "false alarm rate", y = "prevented fraction") + 
    scale_y_continuous(limits = c(0, 1)) + 
    scale_x_continuous(limits = c(0, 1)) + 
    geom_text_repel(aes(label = cutoff))
  print(plot)
  ggsave(plot = plot, file = paste("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS AMOC",
                                   country_list[i], ".jpeg", sep=' '))
}

for(i in seq_along(country_list)){
  plot <- ggplot(filter(metrics_df_HM_wide, country == country_list[i]), aes(x = FAR, y = frac_prevented)) + 
    geom_point(size = 3, shape = 19) + 
    labs(title = paste("HealthMap ROC curve for", country_list[i]), x = "false alarm rate", y = "prevented fraction") + 
    scale_y_continuous(limits = c(-0, 1)) + 
    scale_x_continuous(limits = c(-0, 1)) + 
    geom_text_repel(aes(label = cutoff))
  print(plot)
  ggsave(plot = plot, file = paste("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HM AMOC",
                                   country_list[i], ".jpeg", sep=' '))
}



