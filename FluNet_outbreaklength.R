library(plyr)
library(lubridate)
library(surveillance)
library(ggplot2)
library(dplyr)
library(bcp)
library(ecp)
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


# filter for USA and Ecuador as illustrative examples
FluNet_USA <- filter(FluNet_data, Country == "United States") %>% dplyr::select(c(SDATE, ALL_INF))

FluNet_Ecu <- filter(FluNet_data, Country == "Ecuador") %>% dplyr::select(c(SDATE, ALL_INF))


## bcp on time series to find beginning of outbreaks
USA_bcp <- bcp(FluNet_USA$ALL_INF, burnin = 100, mcmc = 5000)
FluNet_USA <- FluNet_USA %>% mutate(bcp.postprob = USA_bcp$posterior.prob)

FluNet_USA_allcp <- FluNet_USA %>% mutate(changepoint = ifelse(bcp.postprob > 0.5, SDATE, NA))
ggplot(data = FluNet_USA_allcp, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with posterior probability of change point") + 
  geom_vline(xintercept = na.omit(FluNet_USA_allcp$changepoint), lty = 2, col = "red")


FluNet_USA$bcp_start <- NA
for(i in 10:nrow(FluNet_USA)){
  if(FluNet_USA$bcp.postprob[i] >= 0.5 & FluNet_USA$bcp.postprob[(i-1)] < 0.5 & # transition from non-epidemic to epidemic
     mean(FluNet_USA$ALL_INF[(i-5):(i+5)]) > mean(FluNet_USA$ALL_INF [(i-6):(i+4)]) & # running mean to ensure that curve is rising
     sum(!is.na(FluNet_USA$bcp_start[(i-10):(i-1)])) == 0){ # No outbreak flagged during the previous 5 weeks
    FluNet_USA$bcp_start[i] <- FluNet_USA$SDATE[i]
  } else {
    FluNet_USA$bcp_start[i] <- NA
  }
}
ggplot(data = FluNet_USA, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with posterior probability of change point (first epidemic point)") + 
  geom_vline(xintercept = na.omit(FluNet_USA$bcp_start), lty = 2, col = "red")

#Ecuador
Ecu_bcp <- bcp(FluNet_Ecu$ALL_INF, burnin = 100, mcmc = 5000)
FluNet_Ecu <- FluNet_Ecu %>% mutate(bcp.postprob = Ecu_bcp$posterior.prob)

FluNet_Ecu_allcp <- FluNet_Ecu %>% mutate(changepoint = ifelse(bcp.postprob > 0.5, SDATE, NA))
ggplot(data = FluNet_Ecu_allcp, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with posterior probability of change point",
       caption = "Dashed vertical lines mark all change points") + 
  geom_vline(xintercept = na.omit(FluNet_Ecu_allcp$changepoint), lty = 2, col = "red")


FluNet_Ecu$bcp_start <- NA
for(i in 10:nrow(FluNet_Ecu)){
  if(FluNet_Ecu$bcp.postprob[i] >= 0.5 & FluNet_Ecu$bcp.postprob[i-1] < 0.5 & 
     mean(FluNet_Ecu$ALL_INF[(i-5):(i+5)], na.rm = TRUE) > mean(FluNet_Ecu$ALL_INF [(i-6):(i+4)], na.rm = TRUE) & 
     sum(!is.na(FluNet_Ecu$bcp_start[(i-10):(i-1)])) == 0){
    FluNet_Ecu$bcp_start[i] <- FluNet_Ecu$SDATE[i]
  } else {
    FluNet_Ecu$bcp_start[i] <- NA
  }
}
ggplot(data = FluNet_Ecu, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with posterior probability of change point") + 
  geom_vline(xintercept = na.omit(FluNet_Ecu$bcp_start), lty = 2, col = "red")


## reverse timeline and apply bcp to get end of the epidemic
# FluNet_USA$rev_count <- rev(FluNet_USA$ALL_INF)
# FluNet_USA$rev_SDATE <- rev(FluNet_USA$SDATE)
# 
# USA_rev_bcp <- bcp(FluNet_USA$rev_count, burnin = 100, mcmc = 5000)
# FluNet_USA <- FluNet_USA %>% mutate(rev_bcp.postprob = USA_rev_bcp$posterior.prob)
# 
# FluNet_USA_allcp <- FluNet_USA_allcp %>% mutate(rev_changepoint = ifelse(rev_bcp.postprob > 0.5, rev_SDATE, NA))
# ggplot(data = FluNet_USA_rev_allcp, aes(x = SDATE, y = ALL_INF)) + 
#   geom_line(aes(col = rev_bcp.postprob), size = 0.75) +
#   scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
#   labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with posterior probability of change point (reverse)") + 
#   geom_vline(xintercept = na.omit(FluNet_USA_allcp$rev_changepoint), lty = 2, col = "red")
# 
# # bcp finds exactly the same change points, regardless of the order of the timeline (forward or reverse)
# rev(FluNet_USA_allcp$rev_changepoint) == FluNet_USA_allcp$changepoint


# find last changepoint in time series 
FluNet_USA$bcp_end <- NA
for(i in (nrow(FluNet_USA)-9):4){ # run in reverse because otherwise, third criterion cannot be recognized
  if(FluNet_USA$bcp.postprob[i] >= 0.5 & FluNet_USA$bcp.postprob[(i+1)] < 0.5 & # transition from epidemic to non-epidemic
     mean(FluNet_USA$ALL_INF[(i-3):(i+3)], na.rm = TRUE) > mean(FluNet_USA$ALL_INF [(i-2):(i+4)], na.rm = TRUE) & # running mean to ensure that curve is falling
     sum(!is.na(FluNet_USA$bcp_end[(i+1):(i+9)])) == 0){ # No outbreak flagged during the following 9 weeks
    FluNet_USA$bcp_end[i] <- FluNet_USA$SDATE[i]
  } else {
    FluNet_USA$bcp_end[i] <- NA
  }
}
ggplot(data = FluNet_USA, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with posterior probability of change point (last epidemic point)") + 
  geom_vline(xintercept = na.omit(FluNet_USA$bcp_end), lty = 2, col = "red")


# plot time series with start and end point of epidemic 
ggplot(data = FluNet_USA, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with first and last epidemic points") + 
  geom_vline(xintercept = na.omit(FluNet_USA$bcp_start), lty = 2, col = "red") +
  geom_vline(xintercept = na.omit(FluNet_USA$bcp_end), lty = 2, col = "darkgreen")


## do the same for Ecuador
FluNet_Ecu$bcp_end <- NA
for(i in (nrow(FluNet_Ecu)-9):4){ # run in reverse because otherwise, third criterion cannot be recognized
  if(FluNet_Ecu$bcp.postprob[i] >= 0.5 & FluNet_Ecu$bcp.postprob[(i+1)] < 0.5 & # transition from epidemic to non-epidemic
     mean(FluNet_Ecu$ALL_INF[(i-3):(i+3)], na.rm = TRUE) > mean(FluNet_Ecu$ALL_INF [(i-2):(i+4)], na.rm = TRUE) & # running mean to ensure that curve is falling
     sum(!is.na(FluNet_Ecu$bcp_end[(i+1):(i+9)])) == 0){ # No outbreak flagged during the following 9 weeks
    FluNet_Ecu$bcp_end[i] <- FluNet_Ecu$SDATE[i]
  } else {
    FluNet_Ecu$bcp_end[i] <- NA
  }
}
ggplot(data = FluNet_Ecu, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with posterior probability of change point (last epidemic point)") + 
  geom_vline(xintercept = na.omit(FluNet_Ecu$bcp_end), lty = 2, col = "red")


# plot time series with start and end point of epidemic 
ggplot(data = FluNet_Ecu, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with first and last epidemic points") + 
  geom_vline(xintercept = na.omit(FluNet_Ecu$bcp_start), lty = 2, col = "red") +
  geom_vline(xintercept = na.omit(FluNet_Ecu$bcp_end), lty = 2, col = "darkgreen")




### assign binary epidemic variable to all time points between start and end dates
FluNet_USA$startend <- ifelse(is.na(FluNet_USA$bcp_start) == FALSE, "start", 
                              ifelse(is.na(FluNet_USA$bcp_end) == FALSE, "end", NA))
# if the series starts with "end" and/or ends with "start", remove these entries 
if(min(which(FluNet_USA$startend == "end")) < min(which(FluNet_USA$startend == "start"))){
  FluNet_USA$startend[min(which(FluNet_USA$startend == "end"))] <- NA
}
if(max(which(FluNet_USA$startend == "start")) > max(which(FluNet_USA$startend == "end"))){
  FluNet_USA$startend[max(which(FluNet_USA$startend == "start"))] <- NA
}

FluNet_USA$epidemic <- NA
FluNet_USA <- FluNet_USA %>% 
  group_by(grp = cumsum(!is.na(startend))) %>% 
  mutate(epidemic = replace(startend, first(startend) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp) 
FluNet_USA$epidemic[which(FluNet_USA$epidemic == "end")] <- TRUE
FluNet_USA$epidemic[which(is.na(FluNet_USA$epidemic) == TRUE)] <- FALSE


## see if it also works with a little more messy Ecuador data
# no, need to clean the data better in advance
FluNet_Ecu$startend <- ifelse(is.na(FluNet_Ecu$bcp_start) == FALSE, "start", 
                              ifelse(is.na(FluNet_Ecu$bcp_end) == FALSE, "end", NA))
# if the series starts with "end" and/or ends with "start", remove these entries 
if(min(which(FluNet_Ecu$startend == "end")) < min(which(FluNet_Ecu$startend == "start"))){
  FluNet_Ecu$startend[min(which(FluNet_Ecu$startend == "end"))] <- NA
}
if(max(which(FluNet_Ecu$startend == "start")) > max(which(FluNet_Ecu$startend == "end"))){
  FluNet_Ecu$startend[max(which(FluNet_Ecu$startend == "start"))] <- NA
}

FluNet_Ecu$epidemic <- NA
FluNet_Ecu <- FluNet_Ecu %>% 
  group_by(grp = cumsum(!is.na(startend))) %>% 
  mutate(epidemic = replace(startend, first(startend) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp) 
FluNet_Ecu$epidemic[which(FluNet_Ecu$epidemic == "end")] <- TRUE
FluNet_Ecu$epidemic[which(is.na(FluNet_Ecu$epidemic) == TRUE)] <- FALSE
