library(ggplot2)
library(plyr)
library(readxl)
library(lubridate)
library(stringr)
library(dplyr)
library(bcp)
library(qcc)
library(surveillance)
library(cpm)
library(tidyr)

# load FluNet data (just for country_list vector)
setwd("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet")
FluNet_data <- list.files(path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet", pattern = ".csv") %>% 
  lapply(read.csv, stringsAsFactors=F) %>% bind_rows()

country_code <- c("Arg", "Aus", "Bra", "Bul", "Chn", "Cri", "Ecu", "Egy", "Fra", "Gbr", "Ger", "Grc",
                  "Ind", "Irn", "Mex", "Nig", "Rus", "Sau", "Swe", "Tha", "Ury", "Usa", "Vnm", "Zaf")
FluNet_data$Country <- revalue(FluNet_data$Country, c("Iran (Islamic Republic of)" = "Iran", "Russian Federation" = "Russia", 
                                                      "United Kingdom of Great Britain and Northern Ireland" = "United Kingdom",
                                                      "United States of America" = "United States", "Viet Nam" = "Vietnam"))

str(FluNet_data)

cols_F <- c("Country", "WHOREGION", "FLUREGION", "TITLE")
cols_P <- c("SDATE", "EDATE")  

FluNet_data[cols_F] <- lapply(FluNet_data[cols_F], as.factor)
FluNet_data[cols_P] <- lapply(FluNet_data[cols_P], as.POSIXct)

country_list <- levels(FluNet_data$Country)

# load EIOS data
EIOSreports <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/iris_flu_analysis/EIOSreports.csv", 
                        comment.char="#", stringsAsFactors=FALSE)
EIOSreports <- EIOSreports %>% rename(Country = mentionedCountries, CountryIso2 = mentionedCountriesIso2) %>% 
  select(-c(X, translatedFullText, translatedTitle))


# filter for 24 countries in target list
EIOSreports$Country[grep("Viet Nam", EIOSreports$Country)] <- "Vietnam"
EIOSreports$Country[grep("United States of America", EIOSreports$Country)] <- "United States"

EIOSreports <- EIOSreports %>% filter(Country %in% country_list)
EIOSreports$Country <- as.factor(EIOSreports$Country)

# format datetime columns correctly
EIOSreports$fetchDate <- gsub("Z", "", EIOSreports$fetchDate) 
EIOSreports$fetchDate <- gsub("T", " ", EIOSreports$fetchDate)
EIOSreports$fetchDate <- str_sub(EIOSreports$fetchDate, end = -9)
EIOSreports$fetchDate <- as.POSIXct(EIOSreports$fetchDate, format="%Y-%m-%d %H:%M:%S")

EIOSreports$importDate <- gsub("Z", "", EIOSreports$importDate) 
EIOSreports$importDate <- gsub("T", " ", EIOSreports$importDate)
EIOSreports$importDate <- as.POSIXct(EIOSreports$importDate, format="%Y-%m-%d %H:%M:%S")

EIOSreports <- EIOSreports[-which(is.na(EIOSreports$importDate)), ] # remove all observations without import date


# summarize into weekly data
EIOS_byweek <- EIOSreports %>% group_by(date = floor_date(importDate, "week", week_start = 1), country = Country, .drop = FALSE) %>%
  summarize(counts=n()) %>% as.data.frame()
EIOS_byweek <- EIOS_byweek[order(EIOS_byweek$country), ]

EIOS_byweek$count_smooth <- NA
for (i in seq_along(country_list)) { 
  EIOS_temp <- filter(EIOS_byweek, EIOS_byweek$country==country_list[i])
  count_smooth_temp <- loess(EIOS_temp$counts ~ as.numeric(EIOS_temp$date), span = 0.15, degree = 2)
  EIOS_byweek$count_smooth[EIOS_byweek$country==country_list[i]] <- count_smooth_temp$fitted
}


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
### start and end of epidemics in EIOS

p_prior <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)
w_prior <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)
p_prior_col <- paste("p_prior", p_prior, sep = "")
w_prior_col <- paste("w_prior", w_prior, sep = "")
EIOS_byweek <- cbind(EIOS_byweek, setNames(lapply(p_prior_col, function(x) x=NA), p_prior_col))
EIOS_byweek <- cbind(EIOS_byweek, setNames(lapply(w_prior_col, function(x) x=NA), w_prior_col))

bcp_list_p <- vector(mode = "list")
bcp_list_w <- vector(mode = "list")
bcp_postprob_list_p <- vector(mode = "list")
bcp_postprob_list_w <- vector(mode = "list")

for (i in seq_along(country_list)) { 
  EIOS_temp <- filter(EIOS_byweek, country == country_list[i])
  
  bcp_list_p <- lapply(p_prior, function(x) bcp(EIOS_temp$counts, p0 = x, burnin = 100, mcmc = 500, w0 = 0.2))
  bcp_postprob_list_p <- lapply(bcp_list_p, '[[', "posterior.prob")
  bcp_postprob_df_p <- data.frame(matrix(unlist(bcp_postprob_list_p), nrow = 109, byrow = FALSE))
  EIOS_byweek[EIOS_byweek$country == country_list[i], 5:10] <- bcp_postprob_df_p
  
  bcp_list_w <- lapply(w_prior, function(x) bcp(EIOS_temp$counts, w0 = x, burnin = 100, mcmc = 500, p0 = 0.2))
  bcp_postprob_list_w <- lapply(bcp_list_w, '[[', "posterior.prob")
  bcp_postprob_df_w <- data.frame(matrix(unlist(bcp_postprob_list_w), nrow = 109, byrow = FALSE))
  EIOS_byweek[EIOS_byweek$country == country_list[i], 11:16] <- bcp_postprob_df_w
}


start_col <- paste("bcp_start", rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)), sep = "_")
end_col <- paste("bcp_end", rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)), sep = "_")
EIOS_byweek <- cbind(EIOS_byweek, setNames(lapply(start_col, function(x) x=NA), start_col))
EIOS_byweek <- cbind(EIOS_byweek, setNames(lapply(end_col, function(x) x=NA), end_col))

bcp_start_list <- vector(mode = "list")
bcp_end_list <- vector(mode = "list")
bcp_start <- NA
bcp_end <- NA

for (i in seq_along(country_list)) {
  EIOS_temp <- filter(EIOS_byweek, country == country_list[i] & is.na(EIOS_byweek$p_prior0.05) == FALSE)
  
  # start of epidemics
  bcp_start_list <- lapply(5:16, function(x) epi_start(EIOS_temp, col = x))
  bcp_start_df <- data.frame(matrix(unlist(bcp_start_list), nrow = 109, byrow = FALSE))
  
  EIOS_byweek[EIOS_byweek$country == country_list[i], 17:28] <- bcp_start_df
  
  # end of epidemics
  bcp_end_list <- lapply(5:16, function(x) epi_end(EIOS_temp, col = x))
  bcp_end_df <- data.frame(matrix(unlist(bcp_end_list), nrow = 109, byrow = FALSE))
  
  EIOS_byweek[EIOS_byweek$country == country_list[i], 29:40] <- bcp_end_df
}

# replace all 0 with NA
EIOS_byweek[17:40] <- sapply(EIOS_byweek[17:40], na_if, 0)

# start end indicator column
patterns <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")

start_end_list <- vector(mode = "list")
for(pattern in patterns){
  EIOS_temp <- EIOS_byweek[, grep(names(EIOS_byweek), pattern = pattern)]
  start_end_list[[pattern]] <- ifelse(is.na(EIOS_temp[1]) == FALSE, "start", 
                                      ifelse(is.na(EIOS_temp[2]) == FALSE, "end", NA))
}

start_end_col <- paste("start_end", patterns, sep = "_")
EIOS_byweek <- cbind(EIOS_byweek, setNames(lapply(start_end_col, function(x) x=NA), start_end_col))
EIOS_byweek[41:52] <- data.frame(matrix(unlist(start_end_list), nrow = 2616, byrow = FALSE), stringsAsFactors = FALSE)

# for spikes: if 'start' is not followed by an 'end' within 30 weeks, add an 'end' directly after - this period should be optimized
# for(i in 1:nrow(EIOS_byweek)){
#   if(EIOS_byweek$date[i] < "2019-09-01"){
#     if(is.na(EIOS_byweek[i, 40:56]) == FALSE){
#       if(EIOS_byweek[i, 40:56] == "start" & sum(EIOS_byweek[((i+1):(i+30)), 40:56] == "end", na.rm = TRUE) == 0){
#         EIOS_byweek[(i+1), 40:56] <- "end"
#       }
#     }
#   }
# }


# epidemic indicator

EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.05))) %>%
  mutate(epidemic_p_0.05 = replace(start_end_p_0.05, first(start_end_p_0.05) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.1))) %>%
  mutate(epidemic_p_0.1 = replace(start_end_p_0.1, first(start_end_p_0.1) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.2))) %>%
  mutate(epidemic_p_0.2 = replace(start_end_p_0.2, first(start_end_p_0.2) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.3))) %>%
  mutate(epidemic_p_0.3 = replace(start_end_p_0.3, first(start_end_p_0.3) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.5))) %>%
  mutate(epidemic_p_0.5 = replace(start_end_p_0.5, first(start_end_p_0.5) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.8))) %>%
  mutate(epidemic_p_0.8 = replace(start_end_p_0.8, first(start_end_p_0.8) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)

EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.05))) %>%
  mutate(epidemic_w_0.05 = replace(start_end_w_0.05, first(start_end_w_0.05) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.1))) %>%
  mutate(epidemic_w_0.1 = replace(start_end_w_0.1, first(start_end_w_0.1) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.2))) %>%
  mutate(epidemic_w_0.2 = replace(start_end_w_0.2, first(start_end_w_0.2) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.3))) %>%
  mutate(epidemic_w_0.3 = replace(start_end_w_0.3, first(start_end_w_0.3) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.5))) %>%
  mutate(epidemic_w_0.5 = replace(start_end_w_0.5, first(start_end_w_0.5) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.8))) %>%
  mutate(epidemic_w_0.8 = replace(start_end_w_0.8, first(start_end_w_0.8) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)

EIOS_byweek<- EIOS_byweek %>% mutate_at(vars(53:64), ~replace(., is.na(.), FALSE))
EIOS_byweek<- EIOS_byweek %>% mutate_at(vars(53:64), ~replace(., . == "end", TRUE))


write.csv(EIOS_byweek, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_epidemic_bcp_p0_w0.csv")



#### calculate evaluation metrics #####
# load epidemic datasets with outbreak indicators
FluNet_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/FluNet_epidemic.csv")
FluNet_epidemic$SDATE <- as.POSIXct(FluNet_epidemic$SDATE)
EIOS_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_epidemic_bcp_p0_w0.csv") %>%
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
    
    detected_list[[j]] <- sapply(53:64, function(x){
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

EIOS_sens_outbreak <- data.frame(matrix(unlist(EIOS_sens_outbreak_list), ncol = 12, byrow = TRUE))
names(EIOS_sens_outbreak) <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")
EIOS_sens_outbreak$country <- country_list
EIOS_sens_outbreak <- select(EIOS_sens_outbreak, country, everything())


### sensitivity per week ###
EIOS_sens_week_list <- vector(mode = "list")

for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  EIOS_sens_list <- sapply(53:64, function(x){
    state_EIOS <- FluNet_EIOS_temp$epidemic
    alarm_EIOS <- EIOS_temp[, x]
    sens_EIOS <- sum(alarm_EIOS[state_EIOS == TRUE], na.rm = TRUE) / length(alarm_EIOS[state_EIOS == TRUE])
    return(sens_EIOS)
  })
  
  EIOS_sens_week_list[[i]] <- EIOS_sens_list
}

EIOS_sens_week <- data.frame(matrix(unlist(EIOS_sens_week_list), ncol = 12, byrow = TRUE))
names(EIOS_sens_week) <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")
EIOS_sens_week$country <- country_list
EIOS_sens_week <- select(EIOS_sens_week, country, everything())


### false alarm rate
EIOS_far_list <- vector(mode = "list")

for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  EIOS_far <- sapply(53:64, function(x){
    state_EIOS <- FluNet_EIOS_temp$epidemic
    alarm_EIOS <- EIOS_temp[, x]
    far_EIOS <- sum(alarm_EIOS[state_EIOS == FALSE], na.rm = TRUE) / length(alarm_EIOS[state_EIOS == FALSE])
    return(far_EIOS)
  })
  
  EIOS_far_list[[i]] <- EIOS_far
}

EIOS_false_alarm_rate <- data.frame(matrix(unlist(EIOS_far_list), ncol = 12, byrow = TRUE))
names(EIOS_false_alarm_rate) <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")
EIOS_false_alarm_rate$country <- country_list
EIOS_false_alarm_rate <- select(EIOS_false_alarm_rate, country, everything())


### PPV ####
EIOS_PPV_list <- vector(mode = "list")

for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  EIOS_ppv <- sapply(53:64, function(x){
    state_EIOS <- FluNet_EIOS_temp$epidemic
    alarm_EIOS <- EIOS_temp[, x]
    ppv_EIOS <- sum(alarm_EIOS[state_EIOS == FALSE], na.rm = TRUE) / sum(alarm_EIOS, na.rm = TRUE)
    return(ppv_EIOS)
  })
  
  EIOS_PPV_list[[i]] <- EIOS_ppv
}

EIOS_PPV <- data.frame(matrix(unlist(EIOS_PPV_list), ncol = 12, byrow = TRUE))
names(EIOS_PPV) <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")
EIOS_PPV$country <- country_list
EIOS_PPV <- select(EIOS_PPV, country, everything())


### timeliness #### 
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
    
    prevented_list[[j]] <- sapply(53:64, function(x){
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

EIOS_timeliness <- data.frame(matrix(unlist(EIOS_timeliness_list), ncol = 12, byrow = TRUE))
names(EIOS_timeliness) <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")
EIOS_timeliness$country <- country_list
EIOS_timeliness <- select(EIOS_timeliness, country, everything())


#### combine all metrics into one df ####
metrics_df_EIOS <- rbind(EIOS_false_alarm_rate, EIOS_PPV, EIOS_sens_outbreak, EIOS_sens_week, EIOS_timeliness)
metrics_df_EIOS$metric <- rep(c("FAR", "PPV", "sens_outbreak", "sens_week", "frac_prevented"), each = 24)
metrics_df_EIOS_long <- pivot_longer(metrics_df_EIOS, cols = 2:13, names_to = "prior", values_to = "value")
metrics_df_EIOS_wide <- pivot_wider(metrics_df_EIOS_long, names_from = "metric", values_from = "value")
metrics_df_EIOS_wide$country <- as.factor(metrics_df_EIOS_wide$country)
metrics_df_EIOS_wide$prior <- factor(metrics_df_EIOS_wide$prior, 
                                   levels = paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_"),
                                   labels = paste(rep(c("p0=", "w0="), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2),
                                                  rep(c(", w0=0.2", ", p0=0.2"), each = 6), sep = ""))
write.csv(metrics_df_EIOS_wide, "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_metrics_bcp_p0_w0.csv", row.names = FALSE)

### make plots ###
metrics_df_EIOS_wide <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_metrics_bcp_p0_w0.csv")

my_colors <- c("grey", "grey", "red", rep("grey", 5), "blue", "grey", "grey", "grey")

ggplot(metrics_df_EIOS_wide, aes(x = country, y = sens_outbreak, fill = prior)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Sensitivity per outbreak of EIOS according to prior choices", y = "Sensitivity per outbreak", x = "") +
  scale_x_discrete(limits = rev(levels(metrics_df_EIOS_wide$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = my_colors)
ggsave(filename = "EIOS_sens_per_outbreak_p0_w0.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

ggplot(metrics_df_EIOS_wide, aes(x = country, y = sens_week, fill = prior)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Sensitivity per week of EIOS according to prior choices", y = "Sensitivity per week", x = "") +
  scale_x_discrete(limits = rev(levels(metrics_df_EIOS_wide$country))) +
  scale_y_continuous(limits = c(0, 1))  +
  scale_fill_manual(values = my_colors)
ggsave(filename = "EIOS_sens_per_week_p0_w0.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

ggplot(metrics_df_EIOS_wide, aes(x = country, y = 1-FAR, fill = prior)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Specificity of EIOS according to prior choices", y = "Specificity", x = "") +
  scale_x_discrete(limits = rev(levels(metrics_df_EIOS_wide$country))) +
  scale_y_continuous(limits = c(0, 1)) 
ggsave(filename = "EIOS_specificity_p0_w0.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

ggplot(metrics_df_EIOS_wide, aes(x = country, y = PPV, fill = prior)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Positive predictive value of EIOS according to prior choices", y = "PPV", x = "") +
  scale_x_discrete(limits = rev(levels(metrics_df_EIOS_wide$country))) +
  scale_y_continuous(limits = c(0, 1)) 
ggsave(filename = "EIOS_PPV_p0_w0.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

ggplot(metrics_df_EIOS_wide, aes(x = country, y = frac_prevented, fill = prior)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Timeliness (prevented fraction) of EIOS according to prior choices", y = "Prevented fraction", x = "") +
  scale_x_discrete(limits = rev(levels(metrics_df_EIOS_wide$country))) +
  scale_y_continuous(limits = c(0, 1)) 
ggsave(filename = "EIOS_frac_prevented_p0_w0.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")



#### see if bcp converges #### 
EIOS_temp <- filter(EIOS_byweek, country == "Argentina")

bcp_list_conv_test <- vector(mode = "list")
for(i in 1:10){
  bcp_list_conv_test[[i]] <- bcp(EIOS_temp$counts, p0 = 0.2, burnin = 100, mcmc = 500, w0 = 0.2)
}

bcp_list_conv_test_postmean <- lapply(bcp_list_conv_test, '[[', "posterior.mean")
bcp_list_conv_test_postmean_df <- data.frame(matrix(unlist(bcp_list_conv_test_postmean), nrow = 109, byrow = FALSE))
bcp_list_conv_test_postmean_df$index <- 1:109

bcp_list_conv_test_postmean_df_long <- pivot_longer(bcp_list_conv_test_postmean_df, cols = 1:10, 
                                                    names_to = "replicate", values_to = "posterior_mean")
bcp_list_conv_test_postmean_df_long$replicate <- factor(bcp_list_conv_test_postmean_df_long$replicate, 
                                                        levels = paste("X", 1:10, sep = ""),
                                                        labels = 1:10)

bcp_list_conv_test_postprob <- lapply(bcp_list_conv_test, '[[', "posterior.prob")
bcp_list_conv_test_postprob_df <- data.frame(matrix(unlist(bcp_list_conv_test_postprob), nrow = 109, byrow = FALSE))
bcp_list_conv_test_postprob_df$index <- 1:109

bcp_list_conv_test_postprob_df_long <- pivot_longer(bcp_list_conv_test_postprob_df, cols = 1:10, 
                                                    names_to = "replicate", values_to = "posterior_prob")
bcp_list_conv_test_postprob_df_long$replicate <- factor(bcp_list_conv_test_postprob_df_long$replicate, 
                                                        levels = paste("X", 1:10, sep = ""),
                                                        labels = 1:10)


ggplot(bcp_list_conv_test_postmean_df_long, aes(x = index, y = posterior_mean, col = replicate)) + 
  geom_point() + 
  labs(title = "Posterior mean of 10 bcp replicates (EIOS), country = Argentina")

ggplot(bcp_list_conv_test_postprob_df_long, aes(x = index, y = posterior_prob, col = replicate)) + 
  geom_point() + 
  labs(title = "Posterior probability of 10 bcp replicates (EIOS), country = Argentina")


bcp_list_conv_test_postprob <- lapply(bcp_list_conv_test, '[[', "posterior.prob")
bcp_list_conv_test_postprob_df <- data.frame(matrix(unlist(bcp_list_conv_test_postprob), nrow = 109, byrow = FALSE))
bcp_list_conv_test_postprob_df$index <- 1:109

bcp_list_conv_test_postprob_df_long <- pivot_longer(bcp_list_conv_test_postprob_df, cols = 1:10, 
                                                    names_to = "replicate", values_to = "posterior_probability")

ggplot(bcp_list_conv_test_postprob_df_long, aes(x = index, y = posterior_probability, col = replicate)) + 
  geom_point() + 
  labs(title = "Posterior probablity of 10 bcp replicates (EIOS)")


# bcp.fun <- function(df, countries = country_list, p_prior, w_prior){
#   bcp.res <- vector(mode = "list")
#   
#   for(i in country_list){
#     df_temp <- filter(df, country == countries[i])
#     
#     .bcp.fun<- function(df, p_prior, w_prior){
#       bcp.res <- bcp(df$counts, p0 = p_prior, w0 = w_prior, burnin = 100, mcmc = 500)
#       return(bcp.res)
#     }
#     
#     bcp.res[[i]] <- mapply(.bcp.fun, p_prior, w_prior, MoreArgs = list(df = df))
#   }
#   return(bcp.res)
# }

p_prior <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8, rep(0.2, 6))
w_prior <- c(rep(0.2, 6), 0.05, 0.1, 0.2, 0.3, 0.5, 0.8)

bcp_list_p <- vector(mode = "list")
bcp_list_w <- vector(mode = "list")
bcp_postprob_list_p <- vector(mode = "list")
bcp_postprob_list_w <- vector(mode = "list")

for (i in seq_along(country_list)) { 
  EIOS_temp <- filter(EIOS_byweek, country == country_list[i])
  
  bcp_list_p <- lapply(p_prior, function(x) bcp(EIOS_temp$counts, p0 = x, burnin = 100, mcmc = 500, w0 = 0.2))
  bcp_postprob_list_p <- lapply(bcp_list_p, '[[', "posterior.prob")
  bcp_postprob_df_p <- data.frame(matrix(unlist(bcp_postprob_list_p), nrow = 109, byrow = FALSE))
  EIOS_byweek[EIOS_byweek$country == country_list[i], 5:10] <- bcp_postprob_df_p
  
  bcp_list_w <- lapply(w_prior, function(x) bcp(EIOS_temp$counts, w0 = x, burnin = 100, mcmc = 500, p0 = 0.2))
  bcp_postprob_list_w <- lapply(bcp_list_w, '[[', "posterior.prob")
  bcp_postprob_df_w <- data.frame(matrix(unlist(bcp_postprob_list_w), nrow = 109, byrow = FALSE))
  EIOS_byweek[EIOS_byweek$country == country_list[i], 11:16] <- bcp_postprob_df_w
}

bcp_list <- vector(mode = "list")
bcp_list <- bcp.fun(EIOS_byweek, countries = country_list, p_prior = p_prior, w_prior = w_prior)


bcp_list[[1]]["posterior.prob", 1]

plot(EIOS_epidemic[EIOS_epidemic$country == "United States", "p_prior0.2"], ylab = "Posterior probability",
     main = "Posterior probability of two identical bcp runs")
points(EIOS_epidemic[EIOS_epidemic$country == "United States", "w_prior0.2"], col = "red")
abline(0.5, 0, lty = 2)
legend(70, 1.04, legend = c("first run", "second run", "cutoff"), col = c("black", "red", "black"), 
       pch = c(21, 21, NA), lty = c(NA, NA, 2))

