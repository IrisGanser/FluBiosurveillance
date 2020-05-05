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

HM_byweek$count_smooth <- NA
for (i in seq_along(country_list)) { 
  HM_temp <- filter(HM_byweek, HM_byweek$country==country_list[i])
  count_smooth_temp <- loess(HM_temp$counts ~ as.numeric(HM_temp$date), span = 0.07, degree = 2)
  HM_byweek$count_smooth[HM_byweek$country==country_list[i]] <- count_smooth_temp$fitted
}

country_list <- levels(HM_byweek$country)


##### functions #####

epi_start <- function(df, col, cutoff = 0.5){# col in number
  for(j in 2:13){ # first few rows without looking back for previous outbreak
    if(df[j, col] >= cutoff & df[(j-1), col] < cutoff &
       df$count_smooth[j] > df$count_smooth[j-1]){
      bcp_start[j] <- df$date[j]
    } else {
      bcp_start[j] <- NA
    }
  }
  for(j in 10:nrow(df)){
    if(df[j, col] >= cutoff & df[(j-1), col] < cutoff & # transition from non-epidemic to epidemic
       df$count_smooth[j] > df$count_smooth[j-1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(bcp_start[(j-10):(j-1)])) == 0){ # No outbreak flagged during the previous 10 weeks
      bcp_start[j] <- df$date[j]
    } else {
      bcp_start[j] <- NA
    }
  }
  return(bcp_start)
  
}

epi_end <- function(df, col, cutoff = 0.5){ 
  for(j in (nrow(df)-2):2){# run in reverse because otherwise, third criterion cannot be recognized
    if(df[j-1, col] >= cutoff & df[j, col] < cutoff & # transition from non-epidemic to epidemic
       df$count_smooth[j] > df$count_smooth[j+1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(bcp_end[(j+1):(j+15)])) == 0){ # No outbreak flagged during the previous 15 weeks
      bcp_end[j] <- df$date[j]
    } else {
      bcp_end[j] <- NA
    }
  }
  return(bcp_end)
}

##### varying bcp p0 and w0 values ##### 
### start and end of epidemics in HealthMap

p_prior <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)
w_prior <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)
p_prior_col <- paste("p_prior", p_prior, sep = "")
w_prior_col <- paste("w_prior", w_prior, sep = "")
HM_byweek <- cbind(HM_byweek, setNames(lapply(p_prior_col, function(x) x=NA), p_prior_col))
HM_byweek <- cbind(HM_byweek, setNames(lapply(w_prior_col, function(x) x=NA), w_prior_col))

bcp_list_p <- vector(mode = "list")
bcp_list_w <- vector(mode = "list")
bcp_postprob_list_p <- vector(mode = "list")
bcp_postprob_list_w <- vector(mode = "list")

for (i in seq_along(country_list)) { 
  HM_temp <- filter(HM_byweek, country == country_list[i])
  
  bcp_list_p <- lapply(p_prior, function(x) bcp(HM_temp$counts, p0 = x, burnin = 100, mcmc = 500))
  bcp_postprob_list_p <- lapply(bcp_list_p, '[[', "posterior.prob")
  bcp_postprob_df_p <- data.frame(matrix(unlist(bcp_postprob_list_p), nrow = 341, byrow = FALSE))
  HM_byweek[HM_byweek$country == country_list[i], 5:10] <- bcp_postprob_df_p
  
  bcp_list_w <- lapply(w_prior, function(x) bcp(HM_temp$counts, w0 = x, burnin = 100, mcmc = 500))
  bcp_postprob_list_w <- lapply(bcp_list_w, '[[', "posterior.prob")
  bcp_postprob_df_w <- data.frame(matrix(unlist(bcp_postprob_list_w), nrow = 341, byrow = FALSE))
  HM_byweek[HM_byweek$country == country_list[i], 11:16] <- bcp_postprob_df_w
}


start_col <- paste("bcp_start", rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)), sep = "_")
end_col <- paste("bcp_end", rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)), sep = "_")
HM_byweek <- cbind(HM_byweek, setNames(lapply(start_col, function(x) x=NA), start_col))
HM_byweek <- cbind(HM_byweek, setNames(lapply(end_col, function(x) x=NA), end_col))

bcp_start_list <- vector(mode = "list")
bcp_end_list <- vector(mode = "list")

for (i in seq_along(country_list)) {
  HM_temp <- filter(HM_byweek, country == country_list[i] & is.na(HM_byweek$p_prior0.05) == FALSE)
  
  # start of epidemics
  bcp_start_list <- lapply(5:16, function(x) epi_start(HM_temp, col = x))
  bcp_start_df <- data.frame(matrix(unlist(bcp_start_list), nrow = 341, byrow = FALSE))
  
  HM_byweek[HM_byweek$country == country_list[i], 17:28] <- bcp_start_df
  
  # end of epidemics
  bcp_end_list <- lapply(5:16, function(x) epi_end(HM_temp, col = x))
  bcp_end_df <- data.frame(matrix(unlist(bcp_end_list), nrow = 341, byrow = FALSE))
  
  HM_byweek[HM_byweek$country == country_list[i], 29:40] <- bcp_end_df
}

# replace all 0 with NA
HM_byweek[17:40] <- sapply(HM_byweek[17:40], na_if, 0)

# start end indicator column
patterns <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")

start_end_list <- vector(mode = "list")
for(pattern in patterns){
  HM_temp <- HM_byweek[, grep(names(HM_byweek), pattern = pattern)]
  start_end_list[[pattern]] <- ifelse(is.na(HM_temp[1]) == FALSE, "start", 
                                      ifelse(is.na(HM_temp[2]) == FALSE, "end", NA))
}

start_end_col <- paste("start_end", patterns, sep = "_")
HM_byweek <- cbind(HM_byweek, setNames(lapply(start_end_col, function(x) x=NA), start_end_col))
HM_byweek[41:52] <- data.frame(matrix(unlist(start_end_list), nrow = 8184, byrow = FALSE), stringsAsFactors = FALSE)

# for spikes: if 'start' is not followed by an 'end' within 30 weeks, add an 'end' directly after - this period should be optimized
# for(i in 1:nrow(HM_byweek)){
#   if(HM_byweek$date[i] < "2019-09-01"){
#     if(is.na(HM_byweek[i, 40:56]) == FALSE){
#       if(HM_byweek[i, 40:56] == "start" & sum(HM_byweek[((i+1):(i+30)), 40:56] == "end", na.rm = TRUE) == 0){
#         HM_byweek[(i+1), 40:56] <- "end"
#       }
#     }
#   }
# }


# epidemic indicator

HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.05))) %>%
  mutate(epidemic_p_0.05 = replace(start_end_p_0.05, first(start_end_p_0.05) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.1))) %>%
  mutate(epidemic_p_0.1 = replace(start_end_p_0.1, first(start_end_p_0.1) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.2))) %>%
  mutate(epidemic_p_0.2 = replace(start_end_p_0.2, first(start_end_p_0.2) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.3))) %>%
  mutate(epidemic_p_0.3 = replace(start_end_p_0.3, first(start_end_p_0.3) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.5))) %>%
  mutate(epidemic_p_0.5 = replace(start_end_p_0.5, first(start_end_p_0.5) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.8))) %>%
  mutate(epidemic_p_0.8 = replace(start_end_p_0.8, first(start_end_p_0.8) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)

HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.05))) %>%
  mutate(epidemic_w_0.05 = replace(start_end_w_0.05, first(start_end_w_0.05) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.1))) %>%
  mutate(epidemic_w_0.1 = replace(start_end_w_0.1, first(start_end_w_0.1) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.2))) %>%
  mutate(epidemic_w_0.2 = replace(start_end_w_0.2, first(start_end_w_0.2) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.3))) %>%
  mutate(epidemic_w_0.3 = replace(start_end_w_0.3, first(start_end_w_0.3) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.5))) %>%
  mutate(epidemic_w_0.5 = replace(start_end_w_0.5, first(start_end_w_0.5) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.8))) %>%
  mutate(epidemic_w_0.8 = replace(start_end_w_0.8, first(start_end_w_0.8) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)

HM_byweek<- HM_byweek %>% mutate_at(vars(53:64), ~replace(., is.na(.), FALSE))
HM_byweek<- HM_byweek %>% mutate_at(vars(53:64), ~replace(., . == "end", TRUE))


write.csv(HM_byweek, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HM_epidemic_bcp_p0_w0.csv")


#### calculate evaluation metrics #####
# load epidemic datasets with outbreak indicators
FluNet_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/FluNet_epidemic.csv")
FluNet_epidemic$SDATE <- as.POSIXct(FluNet_epidemic$SDATE)
HM_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HM_epidemic_bcp_p0_w0.csv") %>%
  select(-"X")
HM_epidemic$date <- as.POSIXct(HM_epidemic$date)

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
    
    detected_list[[j]] <- sapply(53:64, function(x){
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

HM_sens_outbreak <- data.frame(matrix(unlist(HM_sens_outbreak_list), ncol = 12, byrow = TRUE))
names(HM_sens_outbreak) <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")
HM_sens_outbreak$country <- country_list
HM_sens_outbreak <- select(HM_sens_outbreak, country, everything())


### sensitivity per week ###
HM_sens_week_list <- vector(mode = "list")

for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  
  HM_sens_list <- sapply(53:64, function(x){
    state_HM <- FluNet_HM_temp$epidemic
    alarm_HM <- HM_temp[, x]
    sens_HM <- sum(alarm_HM[state_HM == TRUE], na.rm = TRUE) / length(alarm_HM[state_HM == TRUE])
    return(sens_HM)
  })
  
  HM_sens_week_list[[i]] <- HM_sens_list
}

HM_sens_week <- data.frame(matrix(unlist(HM_sens_week_list), ncol = 12, byrow = TRUE))
names(HM_sens_week) <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")
HM_sens_week$country <- country_list
HM_sens_week <- select(HM_sens_week, country, everything())


### false alarm rate
HM_far_list <- vector(mode = "list")

for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  
  HM_far <- sapply(53:64, function(x){
    state_HM <- FluNet_HM_temp$epidemic
    alarm_HM <- HM_temp[, x]
    far_HM <- sum(alarm_HM[state_HM == FALSE], na.rm = TRUE) / length(alarm_HM[state_HM == FALSE])
    return(far_HM)
  })
  
  HM_far_list[[i]] <- HM_far
}

HM_false_alarm_rate <- data.frame(matrix(unlist(HM_far_list), ncol = 12, byrow = TRUE))
names(HM_false_alarm_rate) <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")
HM_false_alarm_rate$country <- country_list
HM_false_alarm_rate <- select(HM_false_alarm_rate, country, everything())


### PPV ####
HM_PPV_list <- vector(mode = "list")

for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  
  HM_ppv <- sapply(53:64, function(x){
    state_HM <- FluNet_HM_temp$epidemic
    alarm_HM <- HM_temp[, x]
    ppv_HM <- sum(alarm_HM[state_HM == FALSE], na.rm = TRUE) / sum(alarm_HM, na.rm = TRUE)
    return(ppv_HM)
  })
  
  HM_PPV_list[[i]] <- HM_ppv
}

HM_PPV <- data.frame(matrix(unlist(HM_PPV_list), ncol = 12, byrow = TRUE))
names(HM_PPV) <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")
HM_PPV$country <- country_list
HM_PPV <- select(HM_PPV, country, everything())


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
    
    prevented_list[[j]] <- sapply(53:64, function(x){
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

HM_timeliness <- data.frame(matrix(unlist(HM_timeliness_list), ncol = 12, byrow = TRUE))
names(HM_timeliness) <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")
HM_timeliness$country <- country_list
HM_timeliness <- select(HM_timeliness, country, everything())


#### combine all metrics into one df ####
metrics_df_HM <- rbind(HM_false_alarm_rate, HM_PPV, HM_sens_outbreak, HM_sens_week, HM_timeliness)
metrics_df_HM$metric <- rep(c("FAR", "PPV", "sens_outbreak", "sens_week", "frac_prevented"), each = 24)
metrics_df_HM_long <- pivot_longer(metrics_df_HM, cols = 2:13, names_to = "prior", values_to = "value")
metrics_df_HM_wide <- pivot_wider(metrics_df_HM_long, names_from = "metric", values_from = "value")
metrics_df_HM_wide$country <- as.factor(metrics_df_HM_wide$country)
metrics_df_HM_wide$prior <- factor(metrics_df_HM_wide$prior, 
                                      levels = paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_"),
                                      labels = paste(rep(c("p0", "w0"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = ": "))

### make plots ###
ggplot(metrics_df_HM_wide, aes(x = country, y = sens_outbreak, fill = prior)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Sensitivity per outbreak of HealthMap according to prior choices", y = "Sensitivity per outbreak", x = "") +
  scale_x_discrete(limits = rev(levels(metrics_df_HM_wide$country))) +
  scale_y_continuous(limits = c(0, 1)) 
ggsave(filename = "sens_per_outbreak_p0_w0.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

ggplot(metrics_df_HM_wide, aes(x = country, y = sens_week, fill = prior)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Sensitivity per week of HealthMap according to prior choices", y = "Sensitivity per week", x = "") +
  scale_x_discrete(limits = rev(levels(metrics_df_HM_wide$country))) +
  scale_y_continuous(limits = c(0, 1)) 
ggsave(filename = "sens_per_week_p0_w0.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

ggplot(metrics_df_HM_wide, aes(x = country, y = 1-FAR, fill = prior)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Specificity of HealthMap according to prior choices", y = "Specificity", x = "") +
  scale_x_discrete(limits = rev(levels(metrics_df_HM_wide$country))) +
  scale_y_continuous(limits = c(0, 1)) 
ggsave(filename = "specificity_p0_w0.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

ggplot(metrics_df_HM_wide, aes(x = country, y = PPV, fill = prior)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Positive predictive value of HealthMap according to prior choices", y = "PPV", x = "") +
  scale_x_discrete(limits = rev(levels(metrics_df_HM_wide$country))) +
  scale_y_continuous(limits = c(0, 1)) 
ggsave(filename = "PPV_p0_w0.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

ggplot(metrics_df_HM_wide, aes(x = country, y = frac_prevented, fill = prior)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Timeliness (prevented fraction) of HealthMap according to prior choices", y = "Prevented fraction", x = "") +
  scale_x_discrete(limits = rev(levels(metrics_df_HM_wide$country))) +
  scale_y_continuous(limits = c(0, 1)) 
ggsave(filename = "frac_prevented_p0_w0.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")
