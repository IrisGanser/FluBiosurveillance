library(plyr)
library(lubridate)
library(surveillance)
library(ggplot2)
library(dplyr)
library(bcp)
library(cpm)

# load FluNet data
setwd("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet")
FluNet_data <- list.files(path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet", pattern = ".csv") %>% 
  lapply(read.csv, stringsAsFactors=F) %>% bind_rows()

country_code <- c("Arg", "Aus", "Bra", "Bul", "Chn", "Cri", "Ecu", "Egy", "Fra", "Gbr", "Ger", "Grc",
                  "Ind", "Irn", "Mex", "Nig", "Rus", "Sau", "Swe", "Tha", "Ury", "Usa", "Vnm", "Zaf")
FluNet_data$Country <- revalue(FluNet_data$Country, c("Iran (Islamic Republic of)" = "Iran", "Russian Federation" = "Russia", 
                                                      "United Kingdom of Great Britain and Northern Ireland" = "United Kingdom",
                                                      "Viet Nam" = "Vietnam"))

str(FluNet_data)

cols_F <- c("Country", "WHOREGION", "FLUREGION", "TITLE")
cols_P <- c("SDATE", "EDATE")  

FluNet_data[cols_F] <- lapply(FluNet_data[cols_F], as.factor)
FluNet_data[cols_P] <- lapply(FluNet_data[cols_P], as.POSIXct)

FluNet_data$Country <- revalue(FluNet_data$Country, c("Iran (Islamic Republic of)" = "Iran", "Russian Federation" = "Russia", 
                                                      "United Kingdom of Great Britain and Northern Ireland" = "United Kingdom",
                                                      "United States of America" = "United States", "Viet Nam" = "Vietnam"))
country_list <- levels(FluNet_data$Country)



### smoothing according to data abundance 
# all countries with cubic  method looks good visually
FluNet_total <- FluNet_data %>% group_by(Country) %>% summarize(total = sum(ALL_INF), max = max(ALL_INF))
FluNet_total$total_cat <- cut(FluNet_total$total, 
                              breaks = c(0, quantile(FluNet_total$total, probs = 0.5), quantile(FluNet_total$total, probs = 0.75), Inf),
                              labels = c("low", "medium", "high"))

poor_quality_countries <- c("Nigeria", "Thailand", "Vietnam")
FluNet_low_countries <- as.character(FluNet_total$Country[FluNet_total$total_cat == "low"])
# FluNet_low_countries <- FluNet_low_countries[!FluNet_low_countries %in% poor_quality_countries]


FluNet_data$count_smooth <- NA
for (i in seq_along(country_list)) { 
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==country_list[i])
  count_smooth_temp <- loess(FluNet_temp$ALL_INF ~ as.numeric(FluNet_temp$SDATE), span = 0.09, degree = 2)
  FluNet_data$count_smooth[FluNet_data$Country==country_list[i]] <- count_smooth_temp$fitted
  
  # if(i %in% FluNet_low_countries){
  #   count_smooth_temp <- loess(FluNet_temp$ALL_INF ~ as.numeric(FluNet_temp$SDATE), span = 0.09, degree = 2)
  #   FluNet_data$count_smooth[FluNet_data$Country==country_list[i]] <- count_smooth_temp$fitted
  # }else{
  #   count_smooth_temp <- loess(FluNet_temp$ALL_INF ~ as.numeric(FluNet_temp$SDATE), span = 0.09, degree = 1)
  #   FluNet_data$count_smooth[FluNet_data$Country==country_list[i]] <- count_smooth_temp$fitted
  # }
}


for (i in seq_along(country_list)) {
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==country_list[i])
  
  plot <- ggplot(FluNet_temp, aes(x = SDATE, y = ALL_INF)) + 
    geom_line(size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "influenza case counts", 
         title = paste("WHO FluNet data for", country_list[i], "with smoothed time trend (span = 0.09)", sep = " ")) + 
    geom_line(aes(y = count_smooth), col = "red")
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet smoothing", country_list[i], ".jpeg", sep=' '))
}


## apply bcp and criteria for beginning/end of epidemics
# bcp
for (i in seq_along(country_list)) { 
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==country_list[i])
  
  bcp_temp <- bcp(FluNet_temp$ALL_INF, burnin = 100, mcmc = 5000)
  FluNet_data$bcp.postprob[FluNet_data$Country==country_list[i]] <- bcp_temp$posterior.prob
}


# apply criteria for start of epidemics
FluNet_data$bcp_start <- NA
for (i in seq_along(country_list)) { 
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==country_list[i] & is.na(FluNet_data$bcp.postprob) == FALSE)
  
  for(j in 4:13){ # first few rows without looking back for previous outbreak
    if(FluNet_temp$bcp.postprob[j] >= 0.5 & FluNet_temp$bcp.postprob[j-1] < 0.5 &
       FluNet_temp$count_smooth[j] > FluNet_temp$count_smooth[j-1]){
      FluNet_temp$bcp_start[j] <- FluNet_temp$SDATE[j]
    } else {
      FluNet_temp$bcp_start[j] <- NA
    }
  }
  for(j in 10:nrow(FluNet_temp)){
    if(FluNet_temp$bcp.postprob[j] >= 0.5 & FluNet_temp$bcp.postprob[(j-1)] < 0.5 & # transition from non-epidemic to epidemic
       FluNet_temp$count_smooth[j] > FluNet_temp$count_smooth[j-1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(FluNet_temp$bcp_start[(j-10):(j-1)])) == 0){ # No outbreak flagged during the previous 10 weeks
      FluNet_temp$bcp_start[j] <- FluNet_temp$SDATE[j]
    } else {
      FluNet_temp$bcp_start[j] <- NA
    }
  }
  FluNet_data$bcp_start[FluNet_data$Country==country_list[i]] <- FluNet_temp$bcp_start
}

# for (i in seq_along(country_list)) { 
#   FluNet_temp <- filter(FluNet_data, FluNet_data$Country==country_list[i] & is.na(FluNet_data$bcp.postprob) == FALSE)
#   
#   plot <- ggplot(FluNet_temp, aes(x = SDATE, y = ALL_INF)) + 
#     geom_line(aes(col = bcp.postprob), size = 0.75) +
#     scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
#     labs(x = "", y = "influenza case counts", 
#        title = paste("WHO FluNet data for", country_list[i], "with posterior probability of change point", sep = " ")) + 
#     geom_vline(xintercept = na.omit(FluNet_temp$bcp_start), lty = 2, col = "red")
# 
#   print(plot)
# }


#  apply criteria for end of epidemics
FluNet_data$bcp_end <- NA
for (i in seq_along(country_list)) { 
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==country_list[i] & is.na(FluNet_data$bcp.postprob) == FALSE)
  
  for(j in (nrow(FluNet_temp)- 2):2){# run in reverse because otherwise, third criterion cannot be recognized
    if(FluNet_temp$bcp.postprob[j-1] >= 0.5 & FluNet_temp$bcp.postprob[j] < 0.5 & # transition from non-epidemic to epidemic
       FluNet_temp$count_smooth[j] > FluNet_temp$count_smooth[j+1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(FluNet_temp$bcp_end[(j+1):(j+15)])) == 0){ # No outbreak flagged during the previous 15 weeks
      FluNet_temp$bcp_end[j] <- FluNet_temp$SDATE[j]
    } else {
      FluNet_temp$bcp_end[j] <- NA
    }
  }
  FluNet_data$bcp_end[FluNet_data$Country==country_list[i]] <- FluNet_temp$bcp_end
}

# for (i in seq_along(country_list)) { 
# FluNet_temp <- filter(FluNet_data, FluNet_data$Country==country_list[i] & is.na(FluNet_data$bcp.postprob) == FALSE)
# plot <- ggplot(FluNet_temp, aes(x = SDATE, y = ALL_INF)) + 
#         geom_line(aes(col = bcp.postprob), size = 0.75) +
#         scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") +
#           labs(x = "", y = "influenza case counts",
#           title = paste("WHO FluNet data for", country_list[i], "with posterior probability of change point", sep = " ")) +
#         geom_vline(xintercept = na.omit(FluNet_temp$bcp_end), lty = 2, col = "red")
# 
# print(plot)
# }


FluNet_data$startend <- ifelse(is.na(FluNet_data$bcp_start) == FALSE, "start", 
                              ifelse(is.na(FluNet_data$bcp_end) == FALSE, "end", NA))
# if the series starts with "end", remove these entries 
if(min(which(FluNet_data$startend == "end")) < min(which(FluNet_data$startend == "start"))){
  FluNet_data$startend[min(which(FluNet_data$startend == "end"))] <- NA
}

# if 'start' is not followed by an 'end' within 30 weeks, add an 'end' directly after - this period should be optimized
for(i in 1:(nrow(FluNet_data))){
  if(FluNet_data$SDATE[i] < "2019-09-01"){
    if(is.na(FluNet_data$startend[i]) == FALSE){
      if(FluNet_data$startend[i] == "start" & sum(FluNet_data$startend[(i+1):(i+40)] == "end", na.rm = TRUE) == 0){
        FluNet_data$startend[i+1] <- "end"
      }
    }
  }
}

# look at plots and manually delete some starts or ends that don't fit
which(FluNet_data$Country == "Australia" & FluNet_data$startend == "start")
FluNet_data$startend[671] <- NA

which(FluNet_data$Country == "Ecuador" & FluNet_data$startend == "end")
FluNet_data$startend[2212] <- NA

which(FluNet_data$Country == "Egypt" & FluNet_data$startend == "end")
FluNet_data$startend[2602] <- NA

which(FluNet_data$Country == "Germany" & FluNet_data$startend == "end")
FluNet_data$startend[3656] <- NA

which(FluNet_data$Country == "Saudi Arabia" & FluNet_data$startend == "end")
FluNet_data$startend[6058] <- NA
FluNet_data$startend[6096] <- NA

which(FluNet_data$Country == "South Africa" & FluNet_data$startend == "end")
FluNet_data$startend[7937] <- NA

which(FluNet_data$Country == "United Kingdom" & FluNet_data$startend == "end")
FluNet_data$startend[3358] <- NA

which(FluNet_data$Country == "Uruguay" & FluNet_data$startend == "end")
FluNet_data$startend[6959] <- NA

# repeat date copying into bcp_start and bcp_end columns to correct for inserted 'ends' (and manually curated values)
FluNet_data$bcp_end <- ifelse(FluNet_data$startend == "end", FluNet_data$SDATE, NA)
FluNet_data$bcp_start <- ifelse(FluNet_data$startend == "start", FluNet_data$SDATE, NA)


for (i in seq_along(country_list)) {
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==country_list[i] & is.na(FluNet_data$bcp.postprob) == FALSE)
  plot <- ggplot(FluNet_temp, aes(x = SDATE, y = ALL_INF)) +
            geom_line(aes(col = bcp.postprob), size = 0.75) +
            scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") +
              labs(x = "", y = "influenza case counts",
              title = paste("WHO FluNet data for", country_list[i], "with first and last epidemic points", sep = " ")) +
            geom_vline(xintercept = na.omit(FluNet_temp$bcp_start[FluNet_temp$startend == "start"]), lty = 2, col = "red") +
            geom_vline(xintercept = na.omit(FluNet_temp$bcp_end[FluNet_temp$startend == "end"]), lty = 2, col = "darkgreen")

  print(plot)
  #ggsave(plot = plot, file = paste("FluNet outbreak", country_list[i], ".jpeg", sep=' '))
}


## outbreak detection for Nigeria, Thailand, and Vietnam
FluNet_Nig <- filter(FluNet_data, Country == "Nigeria")
FluNet_Tha <- filter(FluNet_data, Country == "Thailand")
FluNet_Vnm <- filter(FluNet_data, Country == "Vietnam") 

#cpm with Mann-Whitney method
cpm_Nig <- processStream(FluNet_Nig$ALL_INF, cpmType = "Mann-Whitney", startup = 10, ARL0 = 500)
FluNet_Nig$cpm <- NA
FluNet_Nig$cpm[cpm_Nig$changePoints] <- FluNet_Nig$SDATE[cpm_Nig$changePoints]

FluNet_Nig$cpm_start <- NA
for(i in 10:nrow(FluNet_Nig)){
  if(is.na(FluNet_Nig$cpm[i]) == FALSE & 
     FluNet_Nig$count_smooth[i] > FluNet_Nig$count_smooth[i-1] &
     sum(!is.na(FluNet_Nig$cpm_start[(i-10):(i-1)])) == 0){
    FluNet_Nig$cpm_start[i] <- FluNet_Nig$SDATE[i]
  } else {
    FluNet_Nig$cpm_start[i] <- NA
  }
}

FluNet_Nig$cpm_end <- NA
for(i in (nrow(FluNet_Nig)-9):1){
  if(is.na(FluNet_Nig$cpm[i]) == FALSE & 
     FluNet_Nig$count_smooth[i] > FluNet_Nig$count_smooth[i+1] &
     sum(!is.na(FluNet_Nig$cpm_end[(i+1):(i+17)])) == 0){
    FluNet_Nig$cpm_end[i] <- FluNet_Nig$SDATE[i]
  } else {
    FluNet_Nig$cpm_end[i] <- NA
  }
}

ggplot(data = FluNet_Nig, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Nigeria with change points, cpm package, 'Mann-Whitney' algorithm") + 
  geom_vline(xintercept = na.omit(FluNet_Nig$cpm_start), lty = 2, col = "red") + 
  geom_vline(xintercept = na.omit(FluNet_Nig$cpm_end), lty = 2, col = "darkgreen")


# Thailand
cpm_Tha <- processStream(FluNet_Tha$ALL_INF, cpmType = "Mann-Whitney", startup = 10, ARL0 = 500)

FluNet_Tha$cpm <- NA
FluNet_Tha$cpm[cpm_Tha$changePoints] <- FluNet_Tha$SDATE[cpm_Tha$changePoints]

ggplot(data = FluNet_Tha, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Thailand with change points, cpm package, 'Mann-Whitney' algorithm") + 
  geom_vline(xintercept = na.omit(FluNet_Tha$cpm), lty = 2, col = "red") 

FluNet_Tha$cpm_start <- NA
for(i in 10:nrow(FluNet_Tha)){
  if(is.na(FluNet_Tha$cpm[i]) == FALSE & 
     FluNet_Tha$count_smooth[i] > FluNet_Tha$count_smooth[i-1] &
     sum(!is.na(FluNet_Tha$cpm_start[(i-10):(i-1)])) == 0){
    FluNet_Tha$cpm_start[i] <- FluNet_Tha$SDATE[i]
  } else {
    FluNet_Tha$cpm_start[i] <- NA
  }
}

FluNet_Tha$cpm_end <- NA
for(i in nrow(FluNet_Tha):1){
  if(is.na(FluNet_Tha$cpm[i]) == FALSE & 
     FluNet_Tha$count_smooth[i] > FluNet_Tha$count_smooth[i+1] &
     sum(!is.na(FluNet_Tha$cpm_end[(i+1):(i+17)])) == 0){
    FluNet_Tha$cpm_end[i] <- FluNet_Tha$SDATE[i]
  } else {
    FluNet_Tha$cpm_end[i] <- NA
  }
}

cpm_Tha$changePoints[15]
FluNet_Tha$cpm_start[232] <- FluNet_Tha$SDATE[232] # manually insert that one 'start' because somehow it is not recognized by my filtering

ggplot(data = FluNet_Tha, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Thailand with change points, cpm package, 'Mann-Whitney' algorithm") + 
  geom_vline(xintercept = na.omit(FluNet_Tha$cpm_start), lty = 2, col = "red") + 
  geom_vline(xintercept = na.omit(FluNet_Tha$cpm_end), lty = 2, col = "darkgreen")


# Vietnam
cpm_Vnm <- processStream(FluNet_Vnm$ALL_INF, cpmType = "Mann-Whitney", startup = 10, ARL0 = 500)

FluNet_Vnm$cpm <- NA
FluNet_Vnm$cpm[cpm_Vnm$changePoints] <- FluNet_Vnm$SDATE[cpm_Vnm$changePoints]


FluNet_Vnm$cpm_start <- NA
for(i in 10:nrow(FluNet_Vnm)){
  if(is.na(FluNet_Vnm$cpm[i]) == FALSE & 
     FluNet_Vnm$count_smooth[i] > FluNet_Vnm$count_smooth[i-1] &
     sum(!is.na(FluNet_Vnm$cpm_start[(i-10):(i-1)])) == 0){
    FluNet_Vnm$cpm_start[i] <- FluNet_Vnm$SDATE[i]
  } else {
    FluNet_Vnm$cpm_start[i] <- NA
  }
}

FluNet_Vnm$cpm_end <- NA
for(i in nrow(FluNet_Vnm):1){
  if(is.na(FluNet_Vnm$cpm[i]) == FALSE & 
     FluNet_Vnm$count_smooth[i] > FluNet_Vnm$count_smooth[i+1] &
     sum(!is.na(FluNet_Vnm$cpm_end[(i+1):(i+17)])) == 0){
    FluNet_Vnm$cpm_end[i] <- FluNet_Vnm$SDATE[i]
  } else {
    FluNet_Vnm$cpm_end[i] <- NA
  }
}

# remove two spike "starts"
which(!is.na(FluNet_Vnm$cpm_start))
FluNet_Vnm$cpm_start[c(282, 325)] <- NA

ggplot(data = FluNet_Vnm, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Vietnam with change points, cpm package, 'Mann-Whitney' algorithm") + 
  geom_vline(xintercept = na.omit(FluNet_Vnm$cpm_start), lty = 2, col = "red") + 
  geom_vline(xintercept = na.omit(FluNet_Vnm$cpm_end), lty = 2, col = "darkgreen")


# startend indicator in FluNet_data data frame
FluNet_data$startend[FluNet_data$Country == "Nigeria"] <- ifelse(is.na(FluNet_Nig$cpm_start) == FALSE, "start", 
                               ifelse(is.na(FluNet_Nig$cpm_end) == FALSE, "end", NA))
FluNet_data$startend[FluNet_data$Country == "Thailand"] <- ifelse(is.na(FluNet_Tha$cpm_start) == FALSE, "start", 
                                                                 ifelse(is.na(FluNet_Tha$cpm_end) == FALSE, "end", NA))
FluNet_data$startend[FluNet_data$Country == "Vietnam"] <- ifelse(is.na(FluNet_Vnm$cpm_start) == FALSE, "start", 
                                                                 ifelse(is.na(FluNet_Vnm$cpm_end) == FALSE, "end", NA))

FluNet_data$bcp_end <- ifelse(FluNet_data$startend == "end", FluNet_data$SDATE, NA)
FluNet_data$bcp_start <- ifelse(FluNet_data$startend == "start", FluNet_data$SDATE, NA)


## set binary 'epidemic' indicator
FluNet_data$epidemic <- NA
FluNet_data <- FluNet_data %>% 
  group_by(Country, grp = cumsum(!is.na(startend))) %>% 
  mutate(epidemic = replace(startend, first(startend) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp) 
FluNet_data$epidemic[which(FluNet_data$epidemic == "end")] <- TRUE
FluNet_data$epidemic[which(is.na(FluNet_data$epidemic) == TRUE)] <- FALSE



for (i in seq_along(country_list)) {
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==country_list[i] & is.na(FluNet_data$bcp.postprob) == FALSE)
  plot <- ggplot(FluNet_temp, aes(x = SDATE, y = ALL_INF, col = epidemic)) +
    geom_line(aes(group = 1), size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") +
    labs(x = "", y = "influenza case counts",
         title = paste("WHO FluNet data for", country_list[i], "with epidemics", sep = " ")) +
    geom_vline(xintercept = na.omit(FluNet_temp$bcp_start[FluNet_temp$startend == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(FluNet_temp$bcp_end[FluNet_temp$startend == "end"]), lty = 2, col = "darkgreen") + 
    scale_color_manual(values = c("#6e6868", "#e64040"))
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet outbreak", country_list[i], ".jpeg", sep=' '))
}

# complete dataframe with epidemic indicator
FluNet_epidemic <- FluNet_data %>% select(c("Country", "WHOREGION", "FLUREGION", "Year", "Week", "SDATE",
                                            "INF_A", "INF_B", "ALL_INF", "ALL_INF2",
                                            "count_smooth", "bcp.postprob", "bcp_start", "bcp_end", "startend", "epidemic"))

write.csv(FluNet_epidemic, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/FluNet_epidemic.csv")



FluNet_epidemic %>% group_by(Country) %>% summarize(n_outbreak_weeks = sum(epidemic == TRUE), 
                                                    n_non_outbreak_week = sum(epidemic == FALSE))

outbreak_length <- FluNet_epidemic %>% filter(epidemic == TRUE) %>% group_by(Country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n()) %>% mutate(outbreak_no = row_number()) %>% select(-grp)
