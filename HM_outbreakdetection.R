library(plyr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(bcp)
library(readxl)
library(qcc)
library(surveillance)
library(cpm)


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


HM_byweek <- dataHM %>% group_by(date = floor_date(load_date, "week"), country = country, .drop = FALSE) %>%
  summarize(counts=n()) %>% as.data.frame()
HM_byweek <- HM_byweek[order(HM_byweek$country), ]

# load indicators df
indicators <- read.csv("country_indicators.csv")

# establish country groups
high_count_countries <- as.character(indicators$country[indicators$HM_total_cat == "high"])
medium_count_countries <- as.character(indicators$country[indicators$HM_total_cat == "medium"])
low_count_countries <- as.character(indicators$country[indicators$HM_total_cat == "low"])
country_list <- levels(HM_byweek$country)


### outbreak detection with bcp in high count countries ###
# smooth counts
HM_byweek$count_smooth <- NA
for (i in seq_along(country_list)) { 
  HM_temp <- filter(HM_byweek, HM_byweek$country==country_list[i])
  count_smooth_temp <- loess(HM_temp$counts ~ as.numeric(HM_temp$date), span = 0.07, degree = 2)
  HM_byweek$count_smooth[HM_byweek$country==country_list[i]] <- count_smooth_temp$fitted
}

# plot HM data for high count countries
for (i in seq_along(high_count_countries)) {
  HM_temp <- filter(HM_byweek, country == high_count_countries[i])
  
  plot <- ggplot(HM_temp, aes(x = date, y = counts, col = bcp.postprob)) + 
    geom_line(size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "HealthMap event counts", 
         title = paste("HealthMap data for", high_count_countries[i], "with smoothed time trend (span = 0.09)", sep = " ")) + 
    geom_line(aes(y = count_smooth), col = "red")
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet smoothing", country_list[i], ".jpeg", sep=' '))
}

# apply bcp algorithm
for (i in seq_along(high_count_countries)) { 
  HM_temp <- filter(HM_byweek, country == high_count_countries[i])
  
  bcp_temp <- bcp(HM_temp$counts, burnin = 100, mcmc = 5000)
  HM_byweek$bcp.postprob[HM_byweek$country==high_count_countries[i]] <- bcp_temp$posterior.prob
}

# apply criteria for start of epidemics
HM_byweek$bcp_start <- NA
for (i in seq_along(high_count_countries)) {
  HM_temp <- filter(HM_byweek, country == high_count_countries[i] & is.na(HM_byweek$bcp.postprob) == FALSE)
  
  for(j in 4:13){ # first few rows without looking back for previous outbreak
    if(HM_temp$bcp.postprob[j] >= 0.5 & HM_temp$bcp.postprob[j-1] < 0.5 &
       HM_temp$count_smooth[j] > HM_temp$count_smooth[j-1]){
      HM_temp$bcp_start[j] <- HM_temp$date[j]
    } else {
      HM_temp$bcp_start[j] <- NA
    }
  }
  for(j in 10:nrow(HM_temp)){
    if(HM_temp$bcp.postprob[j] >= 0.5 & HM_temp$bcp.postprob[(j-1)] < 0.5 & # transition from non-epidemic to epidemic
       HM_temp$count_smooth[j] > HM_temp$count_smooth[j-1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(HM_temp$bcp_start[(j-10):(j-1)])) == 0){ # No outbreak flagged during the previous 10 weeks
      HM_temp$bcp_start[j] <- HM_temp$date[j]
    } else {
      HM_temp$bcp_start[j] <- NA
    }
  }
  HM_byweek$bcp_start[HM_byweek$country==high_count_countries[i]] <- HM_temp$bcp_start
}

#  apply criteria for end of epidemics
HM_byweek$bcp_end <- NA
for (i in seq_along(high_count_countries)) { 
  HM_temp <- filter(HM_byweek, country == high_count_countries[i] & is.na(HM_byweek$bcp.postprob) == FALSE)
  
  for(j in nrow(HM_temp):1){# run in reverse because otherwise, third criterion cannot be recognized
    if(HM_temp$bcp.postprob[j] >= 0.5 & HM_temp$bcp.postprob[(j+1)] < 0.5 & # transition from non-epidemic to epidemic
       HM_temp$count_smooth[j] > HM_temp$count_smooth[j+1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(HM_temp$bcp_end[(j+1):(j+15)])) == 0){ # No outbreak flagged during the previous 15 weeks
       HM_temp$bcp_end[j] <- HM_temp$date[j]
    } else {
      HM_temp$bcp_end[j] <- NA
    }
  }
  HM_byweek$bcp_end[HM_byweek$country==high_count_countries[i]] <- HM_temp$bcp_end
}

# start and end indicators of epidemics
HM_byweek$startend <- ifelse(is.na(HM_byweek$bcp_start) == FALSE, "start", 
                               ifelse(is.na(HM_byweek$bcp_end) == FALSE, "end", NA))
# for spikes: if 'start' is not followed by an 'end' within 30 weeks, add an 'end' directly after - this period should be optimized
for(i in 1:nrow(HM_byweek)){
  if(is.na(HM_byweek$startend[i]) == FALSE){
    if(HM_byweek$startend[i] == "start" & sum(HM_byweek$startend[(i+1):(i+30)] == "end", na.rm = TRUE) == 0){
      HM_byweek$startend[i+1] <- "end"
    }
  }
}

# look at plots
for (i in seq_along(high_count_countries)) {
  HM_temp <- filter(HM_byweek, country == high_count_countries[i])
  
  plot <- ggplot(HM_temp, aes(x = date, y = counts, col = bcp.postprob)) + 
    geom_line(size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "HealthMap event counts", 
         title = paste("HealthMap data for", high_count_countries[i], "with smoothed time trend (span = 0.09)", sep = " ")) + 
    geom_vline(xintercept = na.omit(HM_temp$bcp_start[HM_temp$startend == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(HM_temp$bcp_end[HM_temp$startend == "end"]), lty = 2, col = "darkgreen")
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet smoothing", country_list[i], ".jpeg", sep=' '))
}

