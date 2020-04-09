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

# load indicators df
indicators <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/country_indicators.csv")

EIOS_counts <- EIOS_byweek %>% group_by(country) %>% summarise(mean_count = mean(counts), median_counts = median(counts), 
                                                           max_count = max(counts), total_count = sum(counts))

# establish country groups
high_count_countries <- as.character(indicators$country[indicators$EIOS_total_cat == "high"])
medium_count_countries <- as.character(indicators$country[indicators$EIOS_total_cat == "medium"])
low_count_countries <- as.character(indicators$country[indicators$EIOS_total_cat == "low"])

# smooth counts
EIOS_byweek$count_smooth <- NA
for (i in seq_along(country_list)) { 
  EIOS_temp <- filter(EIOS_byweek, EIOS_byweek$country==country_list[i])
  count_smooth_temp <- loess(EIOS_temp$counts ~ as.numeric(EIOS_temp$date), span = 0.15, degree = 2)
  EIOS_byweek$count_smooth[EIOS_byweek$country==country_list[i]] <- count_smooth_temp$fitted
}

bcp_countries <- filter(EIOS_counts, max_count > 50) %>% select(country) %>% unlist(use.names = FALSE) %>% as.character()
non_bcp_countries <- country_list[!country_list %in% bcp_countries]

#### outbreak detection with bcp in higher count countries ####

# plot EIOS data for high count countries
for (i in seq_along(bcp_countries)) {
  EIOS_temp <- filter(EIOS_byweek, country == bcp_countries[i])
  
  plot <- ggplot(EIOS_temp, aes(x = date, y = counts)) + 
    geom_line(size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "EIOS event counts", 
         title = paste("EIOS data for", bcp_countries[i], "with smoothed time trend (span = 0.15)", sep = " ")) + 
    geom_line(aes(y = count_smooth), col = "red")
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet smoothing", country_list[i], ".jpeg", sep=' '))
}

# apply bcp algorithm
for (i in seq_along(bcp_countries)) { 
  EIOS_temp <- filter(EIOS_byweek, country == bcp_countries[i])
  
  bcp_temp <- bcp(EIOS_temp$counts, burnin = 100, mcmc = 5000)
  EIOS_byweek$bcp.postprob[EIOS_byweek$country==bcp_countries[i]] <- bcp_temp$posterior.prob
}

# apply criteria for start of epidemics
EIOS_byweek$bcp_start <- NA
for (i in seq_along(bcp_countries)) {
  EIOS_temp <- filter(EIOS_byweek, country == bcp_countries[i] & is.na(EIOS_byweek$bcp.postprob) == FALSE)
  
  for(j in 2:13){ # first few rows without looking back for previous outbreak
    if(EIOS_temp$bcp.postprob[j] >= 0.5 & EIOS_temp$bcp.postprob[j-1] < 0.5 &
       EIOS_temp$count_smooth[j] > EIOS_temp$count_smooth[j-1]){
      EIOS_temp$bcp_start[j] <- EIOS_temp$date[j]
    } else {
      EIOS_temp$bcp_start[j] <- NA
    }
  }
  for(j in 10:nrow(EIOS_temp)){
    if(EIOS_temp$bcp.postprob[j] >= 0.5 & EIOS_temp$bcp.postprob[(j-1)] < 0.5 & # transition from non-epidemic to epidemic
       EIOS_temp$count_smooth[j] > EIOS_temp$count_smooth[j-1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(EIOS_temp$bcp_start[(j-10):(j-1)])) == 0){ # No outbreak flagged during the previous 10 weeks
      EIOS_temp$bcp_start[j] <- EIOS_temp$date[j]
    } else {
      EIOS_temp$bcp_start[j] <- NA
    }
  }
  EIOS_byweek$bcp_start[EIOS_byweek$country==bcp_countries[i]] <- EIOS_temp$bcp_start
}


#  apply criteria for end of epidemics
EIOS_byweek$bcp_end <- NA
for (i in seq_along(bcp_countries)) { 
  EIOS_temp <- filter(EIOS_byweek, country == bcp_countries[i] & is.na(EIOS_byweek$bcp.postprob) == FALSE)
  
  for(j in (nrow(EIOS_temp)-2):2){# run in reverse because otherwise, third criterion cannot be recognized
    if(EIOS_temp$bcp.postprob[j-1] >= 0.5 & EIOS_temp$bcp.postprob[j] < 0.5 & # transition from non-epidemic to epidemic
       EIOS_temp$count_smooth[j] > EIOS_temp$count_smooth[j+1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(EIOS_temp$bcp_end[(j+1):(j+15)])) == 0){ # No outbreak flagged during the previous 15 weeks
      EIOS_temp$bcp_end[j] <- EIOS_temp$date[j]
    } else {
      EIOS_temp$bcp_end[j] <- NA
    }
  }
  EIOS_byweek$bcp_end[EIOS_byweek$country==bcp_countries[i]] <- EIOS_temp$bcp_end
}

# start and end indicators of epidemics
EIOS_byweek$startend_bcp <- ifelse(is.na(EIOS_byweek$bcp_start) == FALSE, "start", 
                                 ifelse(is.na(EIOS_byweek$bcp_end) == FALSE, "end", NA))

#for spikes: if 'start' is not followed by an 'end' within 30 weeks, add an 'end' directly after - this period should be optimized
for(i in 1:nrow(EIOS_byweek)){
  if(EIOS_byweek$date[i] < "2019-09-01"){
    if(is.na(EIOS_byweek$startend_bcp[i]) == FALSE){
      if(EIOS_byweek$startend_bcp[i] == "start" & sum(EIOS_byweek$startend_bcp[(i+1):(i+30)] == "end", na.rm = TRUE) == 0){
        EIOS_byweek$startend_bcp[i+1] <- "end"
      }
    }
  }
}

# manually insert one 'start' in Iranw hich was omitted by the algorithm because of smoothing
which(EIOS_byweek$country == "Iran" & EIOS_byweek$startend_bcp == "end")
EIOS_byweek$startend_bcp[1373] <- "start"

# repeat date copying into bcp_start and bcp_end columns to correct for inserted 'ends' (and manually curated values)
EIOS_byweek$bcp_end <- ifelse(EIOS_byweek$startend_bcp == "end", EIOS_byweek$date, NA)
EIOS_byweek$bcp_start <- ifelse(EIOS_byweek$startend_bcp == "start", EIOS_byweek$date, NA)

# look at plots
for (i in seq_along(bcp_countries)) {
  EIOS_temp <- filter(EIOS_byweek, country == bcp_countries[i])
  
  plot <- ggplot(EIOS_temp, aes(x = date, y = counts, col = bcp.postprob)) + 
    geom_line(size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "EIOS event counts", 
         title = paste("EIOS data for", bcp_countries[i], "with start and end of epidemics", sep = " ")) + 
    geom_vline(xintercept = na.omit(EIOS_temp$bcp_start[EIOS_temp$startend_bcp == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(EIOS_temp$bcp_end[EIOS_temp$startend_bcp == "end"]), lty = 2, col = "darkgreen")  + 
    geom_line(aes(y = count_smooth), col = "red")
  
  print(plot)
  #ggsave(plot = plot, file = paste("EIOS smoothing", country_list[i], ".jpeg", sep=' '))
}


#### apply cpm with Mann-Whitney algorithm to non-bcp countries ####

# add Russia because bcp looked weird
cpm_countries <- c(non_bcp_countries, "Russia")

# plot EIOS data for low count countries
for (i in seq_along(cpm_countries)) {
  EIOS_temp <- filter(EIOS_byweek, country == cpm_countries[i])
  
  plot <- ggplot(EIOS_temp, aes(x = date, y = counts)) + 
    geom_line(size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "EIOS event counts", 
         title = paste("EIOS data for", cpm_countries[i], "with smoothed time trend (span = 0.15)", sep = " ")) + 
    geom_line(aes(y = count_smooth), col = "red")
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet smoothing", country_list[i], ".jpeg", sep=' '))
}

##cpm did not work for all of the countries


### apply EWMA to non bpc countries ### 
EIOS_byweek$ewma <- NA
for(i in seq_along(non_bcp_countries)){
  EIOS_temp <- filter(EIOS_byweek, country == non_bcp_countries[i])
  ewma_temp <- ewma(EIOS_temp$counts, lambda = 0.5, nsigmas = 3,
       title = paste(non_bcp_countries[i], "EWMA", sep = " "))
  
  violations <- ewma_temp$violations
  EIOS_temp$ewma[violations] <- TRUE
  EIOS_byweek$ewma[EIOS_byweek$country == non_bcp_countries[i]] <- EIOS_temp$ewma
}

for(i in 1:nrow(EIOS_byweek)){
  if(EIOS_byweek$country[i] %in% non_bcp_countries){
    if(is.na(EIOS_byweek$ewma[i]) == TRUE){
      EIOS_byweek$ewma[i] <- FALSE
    } 
  }
}

### apply EARS to non-bcp countries ###
for(i in seq_along(non_bcp_countries)){
  EIOS_temp <- filter(EIOS_byweek, country == non_bcp_countries[i])
  
  EIOS_temp_Counts <- EIOS_temp$counts
  EIOS_temp_Epoch <- as.Date(EIOS_temp$date)
  EIOS_temp_sts <- sts(observed = EIOS_temp_Counts, epoch = EIOS_temp_Epoch, epochAsDate = TRUE)
  
  
  EIOS_temp_outbreak_C2_7 <- earsC(EIOS_temp_sts, control = list(
    method = "C2", baseline = 7, minSigma = 1, alpha = 0.01
  ))
  plot(EIOS_temp_outbreak_C2_7, main = paste("EIOS data for", non_bcp_countries[i], "EARS C2, 7 weeks baseline", sep = " "))
}
## looks very ugly too

#### set binary 'epidemic' indicator ####
EIOS_byweek$epidemic <- NA

EIOS_byweek$epidemic <- ifelse(EIOS_byweek$country %in% non_bcp_countries, EIOS_byweek$ewma, NA)

EIOS_epidemic <- EIOS_byweek %>% filter(country %in% bcp_countries) %>%
  group_by(country, grp = cumsum(!is.na(startend_bcp))) %>% 
  mutate(epidemic = replace(startend_bcp, first(startend_bcp) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp) 
EIOS_epidemic$epidemic[which(EIOS_epidemic$epidemic == "end")] <- TRUE
EIOS_epidemic$epidemic[which(is.na(EIOS_epidemic$epidemic) == TRUE)] <- FALSE

EIOS_epidemic <- rbind(filter(EIOS_byweek, country %in% non_bcp_countries), EIOS_epidemic)
EIOS_epidemic <- EIOS_epidemic[order(EIOS_epidemic$country), ]


# plot
for (i in seq_along(country_list)) {
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  plot <- ggplot(EIOS_temp, aes(x = date, y = counts, col = epidemic)) + 
    geom_line(aes(group = 1), size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "EIOS event counts", 
         title = paste("EIOS data for", country_list[i], "with start and end of epidemics", sep = " ")) + 
    geom_vline(xintercept = na.omit(EIOS_temp$bcp_start[EIOS_temp$startend_bcp == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(EIOS_temp$bcp_end[EIOS_temp$startend_bcp == "end"]), lty = 2, col = "darkgreen") +
    scale_color_manual(values = c("#6e6868", "#e64040"))
  
  print(plot)
  ggsave(plot = plot, file = paste("EIOS outbreak", country_list[i], ".jpeg", sep=' '))
}

# write data in file
#write.csv(EIOS_epidemic, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_epidemic.csv")


##### apply bcp to every country and set 'epidemic' indicator ##### 
# apply bcp algorithm
for (i in seq_along(country_list)) { 
  EIOS_temp <- filter(EIOS_byweek, country == country_list[i])
  
  bcp_temp <- bcp(EIOS_temp$counts, burnin = 100, mcmc = 5000)
  EIOS_byweek$bcp.postprob[EIOS_byweek$country == country_list[i]] <- bcp_temp$posterior.prob
}

# apply criteria for start of epidemics
EIOS_byweek$bcp_start <- NA
for (i in seq_along(country_list)) {
  EIOS_temp <- filter(EIOS_byweek, country == country_list[i] & is.na(EIOS_byweek$bcp.postprob) == FALSE)
  
  for(j in 2:13){ # first few rows without looking back for previous outbreak
    if(EIOS_temp$bcp.postprob[j] >= 0.5 & EIOS_temp$bcp.postprob[j-1] < 0.5 &
       EIOS_temp$count_smooth[j] > EIOS_temp$count_smooth[j-1]){
      EIOS_temp$bcp_start[j] <- EIOS_temp$date[j]
    } else {
      EIOS_temp$bcp_start[j] <- NA
    }
  }
  for(j in 10:nrow(EIOS_temp)){
    if(EIOS_temp$bcp.postprob[j] >= 0.5 & EIOS_temp$bcp.postprob[(j-1)] < 0.5 & # transition from non-epidemic to epidemic
       EIOS_temp$count_smooth[j] > EIOS_temp$count_smooth[j-1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(EIOS_temp$bcp_start[(j-10):(j-1)])) == 0){ # No outbreak flagged during the previous 10 weeks
      EIOS_temp$bcp_start[j] <- EIOS_temp$date[j]
    } else {
      EIOS_temp$bcp_start[j] <- NA
    }
  }
  EIOS_byweek$bcp_start[EIOS_byweek$country==country_list[i]] <- EIOS_temp$bcp_start
}


#  apply criteria for end of epidemics
EIOS_byweek$bcp_end <- NA
for (i in seq_along(country_list)) { 
  EIOS_temp <- filter(EIOS_byweek, country == country_list[i] & is.na(EIOS_byweek$bcp.postprob) == FALSE)
  
  for(j in (nrow(EIOS_temp)-2):2){# run in reverse because otherwise, third criterion cannot be recognized
    if(EIOS_temp$bcp.postprob[j-1] >= 0.5 & EIOS_temp$bcp.postprob[j] < 0.5 & # transition from non-epidemic to epidemic
       EIOS_temp$count_smooth[j] > EIOS_temp$count_smooth[j+1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(EIOS_temp$bcp_end[(j+1):(j+15)])) == 0){ # No outbreak flagged during the previous 15 weeks
      EIOS_temp$bcp_end[j] <- EIOS_temp$date[j]
    } else {
      EIOS_temp$bcp_end[j] <- NA
    }
  }
  EIOS_byweek$bcp_end[EIOS_byweek$country==country_list[i]] <- EIOS_temp$bcp_end
}

# start and end indicators of epidemics
EIOS_byweek$startend_bcp <- ifelse(is.na(EIOS_byweek$bcp_start) == FALSE, "start", 
                                   ifelse(is.na(EIOS_byweek$bcp_end) == FALSE, "end", NA))

#for spikes: if 'start' is not followed by an 'end' within 30 weeks, add an 'end' directly after - this period should be optimized
for(i in 1:nrow(EIOS_byweek)){
  if(EIOS_byweek$date[i] < "2019-09-01"){
    if(is.na(EIOS_byweek$startend_bcp[i]) == FALSE){
      if(EIOS_byweek$startend_bcp[i] == "start" & sum(EIOS_byweek$startend_bcp[(i+1):(i+30)] == "end", na.rm = TRUE) == 0){
        EIOS_byweek$startend_bcp[i+1] <- "end"
      }
    }
  }
}

# manually insert one 'start' in Iranw hich was omitted by the algorithm because of smoothing
which(EIOS_byweek$country == "Iran" & EIOS_byweek$startend_bcp == "end")
EIOS_byweek$startend_bcp[1373] <- "start"

# repeat date copying into bcp_start and bcp_end columns to correct for inserted 'ends' (and manually curated values)
EIOS_byweek$bcp_end <- ifelse(EIOS_byweek$startend_bcp == "end", EIOS_byweek$date, NA)
EIOS_byweek$bcp_start <- ifelse(EIOS_byweek$startend_bcp == "start", EIOS_byweek$date, NA)

# look at plots
for (i in seq_along(country_list)) {
  EIOS_temp <- filter(EIOS_byweek, country == country_list[i])
  
  plot <- ggplot(EIOS_temp, aes(x = date, y = counts, col = bcp.postprob)) + 
    geom_line(size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "EIOS event counts", 
         title = paste("EIOS data for", country_list[i], "with start and end of epidemics", sep = " ")) + 
    geom_vline(xintercept = na.omit(EIOS_temp$bcp_start[EIOS_temp$startend_bcp == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(EIOS_temp$bcp_end[EIOS_temp$startend_bcp == "end"]), lty = 2, col = "darkgreen")  + 
    geom_line(aes(y = count_smooth), col = "red")
  
  print(plot)
  #ggsave(plot = plot, file = paste("EIOS smoothing", country_list[i], ".jpeg", sep=' '))
}

EIOS_byweek$epidemic <- NA

EIOS_epidemic_all_bcp <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(startend_bcp))) %>% 
  mutate(epidemic = replace(startend_bcp, first(startend_bcp) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp) 
EIOS_epidemic_all_bcp$epidemic[which(EIOS_epidemic_all_bcp$epidemic == "end")] <- TRUE
EIOS_epidemic_all_bcp$epidemic[which(is.na(EIOS_epidemic_all_bcp$epidemic) == TRUE)] <- FALSE


# write data in file
#write.csv(EIOS_epidemic_all_bcp, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_epidemic_all_bcp.csv")



##### 17 different cutoffs for bcp ##### 
### start and end of epidemics in EIOS
start_col <- paste("bcp_start", seq(0.1, 0.9, 0.05), sep = "")
end_col <- paste("bcp_end", seq(0.1, 0.9, 0.05), sep = "")
EIOS_byweek <- cbind(EIOS_byweek, setNames(lapply(start_col, function(x) x=NA), start_col))
EIOS_byweek <- cbind(EIOS_byweek, setNames(lapply(end_col, function(x) x=NA), end_col))

bcp_start_list <- vector(mode = "list")
bcp_start <- rep(0, 109)
bcp_end_list <- vector(mode = "list")
bcp_end <- rep(0, 109)

for (i in seq_along(country_list)) {
  EIOS_temp <- filter(EIOS_byweek, country == country_list[i] & is.na(EIOS_byweek$bcp.postprob) == FALSE)
  
  # start of epidemics
  bcp_start_list <- lapply(seq(0.1, 0.9, 0.05), function(x) epi_start(EIOS_temp, cutoff = x))
  bcp_start_df <- data.frame(matrix(unlist(bcp_start_list), nrow = 109, byrow = FALSE))
  names(bcp_start_df) <- paste("bcp_start", seq(0.1, 0.9, 0.05), sep = "")
  
  EIOS_byweek[EIOS_byweek$country == country_list[i], 6:22] <- bcp_start_df
  
  # end of epidemics
  bcp_end_list <- lapply(seq(0.1, 0.9, 0.05), function(x) epi_end(EIOS_temp, cutoff = x))
  bcp_end_df <- data.frame(matrix(unlist(bcp_end_list), nrow = 109, byrow = FALSE))
  names(bcp_end_df) <- paste("bcp_end", seq(0.1, 0.9, 0.05), sep = "")
  
  EIOS_byweek[EIOS_byweek$country == country_list[i], 23:39] <- bcp_end_df
}

# replace all 0 with NA
EIOS_byweek[6:39] <- sapply(EIOS_byweek[6:39], na_if, 0)

# start end indicator column
patterns <- c("0.1$", "0.15", "0.2$", "0.25", "0.3$", "0.35", "0.4$", "0.45", "0.5$", "0.55", "0.6$", "0.65", "0.7$", "0.75",
              "0.8$", "0.85", "0.9$")
start_end_list <- vector(mode = "list")
for(pattern in patterns){
  EIOS_temp <- EIOS_byweek[, grep(names(EIOS_byweek), pattern = pattern)]
  start_end_list[[pattern]] <- ifelse(is.na(EIOS_temp[1]) == FALSE, "start", 
                                      ifelse(is.na(EIOS_temp[2]) == FALSE, "end", NA))
}

start_end_col <- paste("start_end", seq(0.1, 0.9, 0.05), sep = "")
EIOS_byweek <- cbind(EIOS_byweek, setNames(lapply(start_end_col, function(x) x=NA), start_end_col))
EIOS_byweek[40:56] <- data.frame(matrix(unlist(start_end_list), nrow = 2616, byrow = FALSE), stringsAsFactors = FALSE)

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


EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.1))) %>%
  mutate(epidemic0.1 = replace(start_end0.1, first(start_end0.1) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.15))) %>%
  mutate(epidemic0.15 = replace(start_end0.15, first(start_end0.15) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.2))) %>%
  mutate(epidemic0.2 = replace(start_end0.2, first(start_end0.2) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.25))) %>%
  mutate(epidemic0.25 = replace(start_end0.25, first(start_end0.25) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.3))) %>%
  mutate(epidemic0.3 = replace(start_end0.3, first(start_end0.3) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.35))) %>%
  mutate(epidemic0.35 = replace(start_end0.35, first(start_end0.35) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.4))) %>%
  mutate(epidemic0.4 = replace(start_end0.4, first(start_end0.4) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.45))) %>%
  mutate(epidemic0.45 = replace(start_end0.45, first(start_end0.45) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.5))) %>%
  mutate(epidemic0.5 = replace(start_end0.5, first(start_end0.5) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.55))) %>%
  mutate(epidemic0.55 = replace(start_end0.55, first(start_end0.55) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.65))) %>%
  mutate(epidemic0.65 = replace(start_end0.65, first(start_end0.65) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.6))) %>%
  mutate(epidemic0.6 = replace(start_end0.6, first(start_end0.6) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.7))) %>%
  mutate(epidemic0.7 = replace(start_end0.7, first(start_end0.7) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.75))) %>%
  mutate(epidemic0.75 = replace(start_end0.75, first(start_end0.75) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.8))) %>%
  mutate(epidemic0.8 = replace(start_end0.8, first(start_end0.8) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.85))) %>%
  mutate(epidemic0.85 = replace(start_end0.85, first(start_end0.85) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
EIOS_byweek <- EIOS_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end0.9))) %>%
  mutate(epidemic0.9 = replace(start_end0.9, first(start_end0.9) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)


EIOS_byweek<- EIOS_byweek %>% mutate_at(vars(57:73), ~replace(., is.na(.), FALSE))
EIOS_byweek<- EIOS_byweek %>% mutate_at(vars(57:73), ~replace(., . == "end", TRUE))


write.csv(EIOS_byweek, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_epidemic_ROC.csv")
