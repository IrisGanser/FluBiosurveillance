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

# load indicators df
indicators <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/country_indicators.csv")


HM_counts <- HM_byweek %>% group_by(country) %>% summarise(mean_count = mean(counts), median_counts = median(counts), 
                                                           max_count = max(counts), total_count = sum(counts))

# establish country groups
high_count_countries <- as.character(indicators$country[indicators$HM_total_cat == "high"])
medium_count_countries <- as.character(indicators$country[indicators$HM_total_cat == "medium"])
low_count_countries <- as.character(indicators$country[indicators$HM_total_cat == "low"])
country_list <- levels(HM_byweek$country)


### outbreak detection with bcp in high count countries ###
# remove China from high count countries because there is just too much noise in data (no clear seasonality)
high_count_countries <- high_count_countries[high_count_countries != "China"]

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
  
  plot <- ggplot(HM_temp, aes(x = date, y = counts)) + 
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
  
  for(j in 2:13){ # first few rows without looking back for previous outbreak
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
  
  for(j in (nrow(HM_temp)-2):2){# run in reverse because otherwise, third criterion cannot be recognized
    if(HM_temp$bcp.postprob[j-1] >= 0.5 & HM_temp$bcp.postprob[j] < 0.5 & # transition from non-epidemic to epidemic
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
HM_byweek$startend_bcp <- ifelse(is.na(HM_byweek$bcp_start) == FALSE, "start", 
                               ifelse(is.na(HM_byweek$bcp_end) == FALSE, "end", NA))

# for spikes: if 'start' is not followed by an 'end' within 30 weeks, add an 'end' directly after - this period should be optimized
# for(i in 1:nrow(HM_byweek)){
#   if(HM_byweek$date[i] < "2019-09-01"){
#     if(is.na(HM_byweek$startend_bcp[i]) == FALSE){
#       if(HM_byweek$startend_bcp[i] == "start" & sum(HM_byweek$startend_bcp[(i+1):(i+40)] == "end", na.rm = TRUE) == 0){
#         HM_byweek$startend_bcp[i+1] <- "end"
#       }
#     }
#   }
# }

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

# repeat date copying into bcp_start and bcp_end columns to correct for inserted 'ends' (and manually curated values)
HM_byweek$bcp_end <- ifelse(HM_byweek$startend_bcp == "end", HM_byweek$date, NA)
HM_byweek$bcp_start <- ifelse(HM_byweek$startend_bcp == "start", HM_byweek$date, NA)

# look at plots
for (i in seq_along(high_count_countries)) {
  HM_temp <- filter(HM_byweek, country == high_count_countries[i])
  
  plot <- ggplot(HM_temp, aes(x = date, y = counts, col = bcp.postprob)) + 
    geom_line(size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "HealthMap event counts", 
         title = paste("HealthMap data for", high_count_countries[i], "with start and end of epidemics", sep = " ")) + 
    geom_vline(xintercept = na.omit(HM_temp$bcp_start[HM_temp$startend_bcp == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(HM_temp$bcp_end[HM_temp$startend_bcp == "end"]), lty = 2, col = "darkgreen")
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet smoothing", country_list[i], ".jpeg", sep=' '))
}




### apply cpm algorithms to medium count countries (countries with enough counts for cpm)###
medium_count_countries <- c(medium_count_countries, "China")

non_cpm_countries <- filter(HM_counts, max_count <= 5 | total_count < 60) %>% filter(total_count < 100) %>%
  select(country) %>% unlist(use.names = FALSE) %>% as.character()
# for these countries, cpm is not likely to work because not enough counts
cpm_countries <- c(medium_count_countries, low_count_countries[!low_count_countries %in% non_cpm_countries])
cpm_countries <- sort(cpm_countries) #sort alphabetically

# plot HM data for cpm countries
for (i in seq_along(cpm_countries)) {
  HM_temp <- filter(HM_byweek, country == cpm_countries[i])
  
  plot <- ggplot(HM_temp, aes(x = date, y = counts)) + 
    geom_line(size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "HealthMap event counts", 
         title = paste("HealthMap data for", cpm_countries[i], "with smoothed time trend (span = 0.09)", sep = " ")) + 
    geom_line(aes(y = count_smooth), col = "red", size = 0.75)
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet smoothing", country_list[i], ".jpeg", sep=' '))
}


# apply cpm algorithm with Mann-Whitney method
HM_byweek$cpm <- NA
for (i in seq_along(cpm_countries)) {
  HM_temp <- filter(HM_byweek, country == cpm_countries[i])
  
  cpm_temp <- processStream(HM_temp$counts, cpmType = "Mann-Whitney", startup = 10, ARL0 = 500)
  HM_temp$cpm[cpm_temp$changePoints] <- HM_temp$date[cpm_temp$changePoints]
  HM_byweek$cpm[HM_byweek$country == cpm_countries[i]] <- HM_temp$cpm
}

# apply criteria for start of epidemics
HM_byweek$cpm_start <- NA
for (i in seq_along(cpm_countries)) {
  HM_temp <- filter(HM_byweek, country == cpm_countries[i])
  
  for(j in 2:13){ # first few rows without looking back for previous outbreak
    if(is.na(HM_temp$cpm[j]) == FALSE &
       HM_temp$count_smooth[j] > HM_temp$count_smooth[j-1]){
      HM_temp$cpm_start[j] <- HM_temp$date[j]
    } else {
      HM_temp$cpm_start[j] <- NA
    }
  }
  for(j in 10:nrow(HM_temp)){
    if(is.na(HM_temp$cpm[j]) == FALSE & # change point
       HM_temp$count_smooth[j] > HM_temp$count_smooth[j-1] & # smoothed mean to ensure that curve is rising 
       sum(!is.na(HM_temp$cpm_start[(j-10):(j-1)])) == 0){ # No outbreak flagged during the previous 10 weeks
      HM_temp$cpm_start[j] <- HM_temp$date[j]
    } else {
      HM_temp$cpm_start[j] <- NA
    }
  }
  HM_byweek$cpm_start[HM_byweek$country==cpm_countries[i]] <- HM_temp$cpm_start
}



#  apply criteria for end of epidemics
HM_byweek$cpm_end <- NA
for (i in seq_along(cpm_countries)) { 
  HM_temp <- filter(HM_byweek, country == cpm_countries[i])
  
  for(j in nrow(HM_temp):1){# run in reverse because otherwise, third criterion cannot be recognized
    if(is.na(HM_temp$cpm[j]) == FALSE &
       HM_temp$count_smooth[j] > HM_temp$count_smooth[j+1] & 
       sum(!is.na(HM_temp$cpm_end[(j+1):(j+15)])) == 0){ 
      HM_temp$cpm_end[j] <- HM_temp$date[j]
    } else {
      HM_temp$cpm_end[j] <- NA
    }
  }
  HM_byweek$cpm_end[HM_byweek$country==cpm_countries[i]] <- HM_temp$cpm_end
}


# start and end indicators of epidemics
HM_byweek$startend_cpm <- ifelse(is.na(HM_byweek$cpm_start) == FALSE, "start", 
                             ifelse(is.na(HM_byweek$cpm_end) == FALSE, "end", NA))


# look at plots
for (i in seq_along(cpm_countries)) {
  HM_temp <- filter(HM_byweek, country == cpm_countries[i])
  
  plot <- ggplot(HM_temp, aes(x = date, y = counts)) + 
    geom_line(size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "HealthMap event counts", 
         title = paste("HealthMap data for", cpm_countries[i], "with start and end of epidemics", sep = " ")) + 
    geom_vline(xintercept = na.omit(HM_temp$cpm_start[HM_temp$startend_cpm == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(HM_temp$cpm_end[HM_temp$startend_cpm == "end"]), lty = 2, col = "darkgreen")
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet smoothing", country_list[i], ".jpeg", sep=' '))
}



### run EWMA on left countries ###
non_cpm_countries

HM_byweek$ewma <- NA
for (i in seq_along(non_cpm_countries)) {
  HM_temp <- filter(HM_byweek, country == non_cpm_countries[i])
  
  ewma_temp <- ewma(HM_temp$counts, lambda = 0.3, nsigmas = 3,
                    title = paste(non_cpm_countries[i], "EWMA", sep = " "))
  
  violations <- ewma_temp$violations
  HM_temp$ewma[violations] <- TRUE
  HM_byweek$ewma[HM_byweek$country == non_cpm_countries[i]] <- HM_temp$ewma
}

for(i in 1:nrow(HM_byweek)){
  if(HM_byweek$country[i] %in% non_cpm_countries){
    if(is.na(HM_byweek$ewma[i]) == TRUE){
      HM_byweek$ewma[i] <- FALSE
    } 
  }
}


### set binary 'epidemic' indicator

HM_byweek$epidemic <- ifelse(HM_byweek$country %in% non_cpm_countries, HM_byweek$ewma, NA)

# overall startend indicator
HM_byweek$startend <- NA
for(i in 1:nrow(HM_byweek)){
  if(HM_byweek$country[i] %in% high_count_countries){
    HM_byweek$startend[i] <- HM_byweek$startend_bcp[i]
  } else if(HM_byweek$country[i] %in% cpm_countries){
    HM_byweek$startend[i] <- HM_byweek$startend_cpm[i]
  } else{
    HM_byweek$startend[i] <- NA
  }
}


HM_epidemic <- HM_byweek %>% filter(country %in% high_count_countries | country %in% cpm_countries) %>% 
  group_by(country, grp = cumsum(!is.na(startend))) %>% 
  mutate(epidemic = replace(startend, first(startend) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp) 
HM_epidemic$epidemic[which(HM_epidemic$epidemic == "end")] <- TRUE
HM_epidemic$epidemic[which(is.na(HM_epidemic$epidemic) == TRUE)] <- FALSE

HM_epidemic <- rbind(filter(HM_byweek, country %in% non_cpm_countries), HM_epidemic) 
HM_epidemic <- HM_epidemic[order(HM_epidemic$country), ]


for (i in seq_along(country_list)) {
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  
  plot <- ggplot(HM_temp, aes(x = date, y = counts, col = epidemic)) + 
    geom_line(aes(group = 1), size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "HealthMap event counts", 
         title = paste("HealthMap data for", country_list[i], "with start and end of epidemics", sep = " ")) + 
    geom_vline(xintercept = na.omit(HM_temp$cpm_start[HM_temp$startend_cpm == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(HM_temp$cpm_end[HM_temp$startend_cpm == "end"]), lty = 2, col = "darkgreen") + 
    geom_vline(xintercept = na.omit(HM_temp$bcp_start[HM_temp$startend_bcp == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(HM_temp$bcp_end[HM_temp$startend_bcp == "end"]), lty = 2, col = "darkgreen") +
    scale_color_manual(values = c("#6e6868", "#e64040"))

  print(plot)
  #ggsave(plot = plot, file = paste("HealthMap outbreak", country_list[i], ".jpeg", sep=' '))
}

# write data in file
#write.csv(HM_epidemic, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HealthMap_epidemic.csv")



##### apply bcp to all countries ##### 

for (i in seq_along(country_list)) { 
  HM_temp <- filter(HM_byweek, country == country_list[i])
  
  bcp_temp <- bcp(HM_temp$counts, burnin = 100, mcmc = 500)
  HM_byweek$bcp.postprob[HM_byweek$country==country_list[i]] <- bcp_temp$posterior.prob
}

# apply criteria for start of epidemics
HM_byweek$bcp_start <- NA
for (i in seq_along(country_list)) {
  HM_temp <- filter(HM_byweek, country == country_list[i] & is.na(HM_byweek$bcp.postprob) == FALSE)
  
  for(j in 2:13){ # first few rows without looking back for previous outbreak
    if(HM_temp$bcp.postprob[j] >= 0.5 & HM_temp$bcp.postprob[j-1] < 0.5 &
       HM_temp$count_smooth[j] > HM_temp$count_smooth[j-1]){
      HM_temp$bcp_start[j] <- HM_temp$date[j]
    } else {
      HM_temp$bcp_start[j] <- NA
    }
  }
  for(j in 12:nrow(HM_temp)){
    if(HM_temp$bcp.postprob[j] >= 0.5 & HM_temp$bcp.postprob[(j-1)] < 0.5 & # transition from non-epidemic to epidemic
       HM_temp$count_smooth[j] > HM_temp$count_smooth[j-1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(HM_temp$bcp_start[(j-12):(j-1)])) == 0){ # No outbreak flagged during the previous 10 weeks
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

# look at plots
for (i in seq_along(country_list)) {
  HM_temp <- filter(HM_byweek, country == country_list[i])
  
  plot <- ggplot(HM_temp, aes(x = date, y = counts, col = bcp.postprob)) + 
    geom_line(size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "HealthMap event counts", 
         title = paste("HealthMap data for", country_list[i], "with start and end of epidemics", sep = " ")) + 
    geom_vline(xintercept = na.omit(HM_temp$bcp_start[HM_temp$startend_bcp == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(HM_temp$bcp_end[HM_temp$startend_bcp == "end"]), lty = 2, col = "darkgreen")
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet smoothing", country_list[i], ".jpeg", sep=' '))
}



# epidemic indicator
HM_epidemic_all_bcp <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(startend_bcp))) %>% 
  mutate(epidemic = replace(startend_bcp, first(startend_bcp) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp) 
HM_epidemic_all_bcp$epidemic[which(HM_epidemic_all_bcp$epidemic == "end")] <- TRUE
HM_epidemic_all_bcp$epidemic[which(is.na(HM_epidemic_all_bcp$epidemic) == TRUE)] <- FALSE

# write data in file
#write.csv(HM_epidemic_all_bcp, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HealthMap_epidemic_all_bcp.csv")

for (i in seq_along(country_list)) {
  HM_temp <- filter(HM_epidemic_all_bcp, country == country_list[i])
  
  plot <- ggplot(HM_temp, aes(x = date, y = counts, col = epidemic)) + 
    geom_line(aes(group = 1), size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "HealthMap event counts", 
         title = paste("HealthMap data for", country_list[i], "with start and end of epidemics", sep = " ")) + 
    geom_vline(xintercept = na.omit(HM_temp$bcp_start[HM_temp$startend_bcp == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(HM_temp$bcp_end[HM_temp$startend_bcp == "end"]), lty = 2, col = "darkgreen")  + 
    scale_color_manual(values = c("#6e6868", "#e64040"))
  
  print(plot)
  #ggsave(plot = plot, file = paste("HealthMap outbreak", country_list[i], ".jpeg", sep=' '))
}



##### 17 different cutoffs for bcp ##### 
### start and end of epidemics in HealthMap

start_col <- paste("bcp_start", seq(0.1, 0.9, 0.05), sep = "")
end_col <- paste("bcp_end", seq(0.1, 0.9, 0.05), sep = "")
HM_byweek <- cbind(HM_byweek, setNames(lapply(start_col, function(x) x=NA), start_col))
HM_byweek <- cbind(HM_byweek, setNames(lapply(end_col, function(x) x=NA), end_col))

bcp_start_list <- vector(mode = "list")
bcp_start <- rep(0, 341)
bcp_end_list <- vector(mode = "list")
bcp_end <- rep(0, 341)

for (i in seq_along(country_list)) {
  HM_temp <- filter(HM_byweek, country == country_list[i] & is.na(HM_byweek$bcp.postprob) == FALSE)
  
  # start of epidemics
  bcp_start_list <- lapply(seq(0.1, 0.9, 0.05), function(x) epi_start(HM_temp, cutoff = x))
  bcp_start_df <- data.frame(matrix(unlist(bcp_start_list), nrow = 341, byrow = FALSE))
  names(bcp_start_df) <- paste("bcp_start", seq(0.1, 0.9, 0.05), sep = "")
  
  HM_byweek[HM_byweek$country == country_list[i], 6:22] <- bcp_start_df
  
  # end of epidemics
  bcp_end_list <- lapply(seq(0.1, 0.9, 0.05), function(x) epi_end(HM_temp, cutoff = x))
  bcp_end_df <- data.frame(matrix(unlist(bcp_end_list), nrow = 341, byrow = FALSE))
  names(bcp_end_df) <- paste("bcp_end", seq(0.1, 0.9, 0.05), sep = "")
  
  HM_byweek[HM_byweek$country == country_list[i], 23:39] <- bcp_end_df
}

# replace all 0 with NA
HM_byweek[6:39] <- sapply(HM_byweek[6:39], na_if, 0)

# start end indicator column
patterns <- c("0.1$", "0.15", "0.2$", "0.25", "0.3$", "0.35", "0.4$", "0.45", "0.5$", "0.55", "0.6$", "0.65", "0.7$", "0.75",
              "0.8$", "0.85", "0.9$")
start_end_list <- vector(mode = "list")
for(pattern in patterns){
  HM_temp <- HM_byweek[, grep(names(HM_byweek), pattern = pattern)]
  start_end_list[[pattern]] <- ifelse(is.na(HM_temp[1]) == FALSE, "start", 
                                      ifelse(is.na(HM_temp[2]) == FALSE, "end", NA))
}

start_end_col <- paste("start_end", seq(0.1, 0.9, 0.05), sep = "")
HM_byweek <- cbind(HM_byweek, setNames(lapply(start_end_col, function(x) x=NA), start_end_col))
HM_byweek[40:56] <- data.frame(matrix(unlist(start_end_list), nrow = 8184, byrow = FALSE), stringsAsFactors = FALSE)

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

# somehow I can't get this code to work
# for(j in names(HM_byweek)[40:56]){
#   HM_temp <- HM_byweek[, c("country", j)]
#   HM_temp <- HM_temp %>% 
#     group_by(country, grp = cumsum(!is.na(HM_temp[2]))) %>%
#     mutate(epidemic = replace(HM_temp$j, first(HM_temp$j) == 'start', TRUE))
# 
#   HM$epidemic <- ifelse(!is.na(HM_temp[j]) & first(HM_temp[j]) == 'start', replace(HM_temp[j], first(HM_temp[j]) == 'start', TRUE), FALSE)
#   %>%
#     ungroup() %>%
#     select(-grp)
# 
# }


HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.1))) %>%
  mutate(epidemic0.1 = replace(start_end0.1, first(start_end0.1) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.15))) %>%
  mutate(epidemic0.15 = replace(start_end0.15, first(start_end0.15) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.2))) %>%
  mutate(epidemic0.2 = replace(start_end0.2, first(start_end0.2) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.25))) %>%
  mutate(epidemic0.25 = replace(start_end0.25, first(start_end0.25) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.3))) %>%
  mutate(epidemic0.3 = replace(start_end0.3, first(start_end0.3) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.35))) %>%
  mutate(epidemic0.35 = replace(start_end0.35, first(start_end0.35) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.4))) %>%
  mutate(epidemic0.4 = replace(start_end0.4, first(start_end0.4) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.45))) %>%
  mutate(epidemic0.45 = replace(start_end0.45, first(start_end0.45) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.5))) %>%
  mutate(epidemic0.5 = replace(start_end0.5, first(start_end0.5) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.55))) %>%
  mutate(epidemic0.55 = replace(start_end0.55, first(start_end0.55) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.65))) %>%
  mutate(epidemic0.65 = replace(start_end0.65, first(start_end0.65) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.6))) %>%
  mutate(epidemic0.6 = replace(start_end0.6, first(start_end0.6) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.7))) %>%
  mutate(epidemic0.7 = replace(start_end0.7, first(start_end0.7) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.75))) %>%
  mutate(epidemic0.75 = replace(start_end0.75, first(start_end0.75) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.8))) %>%
  mutate(epidemic0.8 = replace(start_end0.8, first(start_end0.8) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.85))) %>%
  mutate(epidemic0.85 = replace(start_end0.85, first(start_end0.85) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.9))) %>%
  mutate(epidemic0.9 = replace(start_end0.9, first(start_end0.9) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)


HM_byweek<- HM_byweek %>% mutate_at(vars(57:73), ~replace(., is.na(.), FALSE))
HM_byweek<- HM_byweek %>% mutate_at(vars(57:73), ~replace(., . == "end", TRUE))


write.csv(HM_byweek, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HM_epidemic_ROC.csv")
