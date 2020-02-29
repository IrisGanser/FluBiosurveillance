library(plyr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(bcp)
library(readxl)
library(qcc)
library(surveillance)
library(AnomalyDetection)
library(BreakoutDetection)
library(cpm)


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

years <- 2012:2020
HM_byweek$season <- cut(HM_byweek$date, 
                          breaks=as.POSIXct(paste(years,"-07-01",sep="")),
                          labels=paste(years[-length(years)],years[-length(years)]+1,sep="/"))



# filter for USA and Ecuador
HM_USA <- filter(HM_byweek, country == "United States")
HM_Ecu <- filter(HM_byweek, country == "Ecuador")


##### outbreak detection with surveillance package algorithms #####

# construct sts objects
HM_USA_Counts <- HM_USA$counts
HM_USA_Epoch <- as.Date(HM_USA$date)
USA_sts <- sts(observed = HM_USA_Counts, epoch = HM_USA_Epoch, epochAsDate = TRUE)

HM_Ecu_Counts <- HM_Ecu$counts
HM_Ecu_Epoch <- as.Date(HM_Ecu$date)
Ecu_sts <- sts(observed = HM_Ecu_Counts, epoch = HM_Ecu_Epoch, epochAsDate = TRUE)

HM_USA_2013 <- filter(HM_USA, date > "2013-04-01" & date < "2014-01-01") 
USA_2013_sts <- sts(observed = HM_USA_2013$counts, epochAsDate = FALSE, start = c(2013, 4), frequency = 52)

HM_Ecu_2013 <- filter(HM_Ecu, date > "2013-04-01" & date < "2014-01-01") 
Ecu_2013_sts <- sts(observed = HM_Ecu_2013$counts, epochAsDate = FALSE, start = c(2013, 4), frequency = 52)

USA_disProg <- sts2disProg(USA_sts)
Ecu_disProg <- sts2disProg(Ecu_sts)

## EARS C2 algorithm
USA_outbreak_C2_7 <- earsC(USA_sts, control = list(
  method = "C2", baseline = 7, minSigma = 1, alpha = 0.01
))
plot(USA_outbreak_C2_7, main = "USA HM data EARS C2, 7 weeks baseline")

USA_outbreak_C2_14 <- earsC(USA_sts, control = list(
  method = "C2", baseline = 14, minSigma = 1, alpha = 0.01
))
plot(USA_outbreak_C2_14, main = "USA HM data EARS C2, 14 weeks baseline")

USA_outbreak_C2_28 <- earsC(USA_sts, control = list(
  method = "C2", baseline = 28, minSigma = 1, alpha = 0.01
))
plot(USA_outbreak_C2_28, main = "USA HM data EARS C2, 28 weeks baseline")
# algorithm is very sensitive and picks up on mid-season variations. Very dependent on choice of baseline length

Ecu_outbreak_C2_7 <- earsC(Ecu_sts, control = list(
  method = "C2", baseline = 7, minSigma = 1, alpha = 0.01
))
plot(Ecu_outbreak_C2_7, main = "Ecuador HM EARS C2, 7 weeks baseline")

Ecu_outbreak_C2_14 <- earsC(Ecu_sts, control = list(
  method = "C2", baseline = 14, minSigma = 1, alpha = 0.01
))
plot(Ecu_outbreak_C2_14, main = "Ecuador HM EARS C2, 7 weeks baseline")
# Ecuador data are very sparse. At least it doesn't alert at one report per week, and one alert per spike in first or second week
# actually doesn't look too bad, might be an option for countries with sparse data


## CUSUM method
kh <- find.kh(ARLa = 500, ARLr = 7)

#USA
USA_outbreak_cusum <- surveillance::cusum(USA_sts, control = list(
  range = 1:length(observed(USA_sts)), k = kh$k, h= kh$h, trans = "rossi", m = NULL
))
plot(USA_outbreak_cusum, main = "USA HM CUSUM with Rossi method")

USA_2013_cusum <- surveillance::cusum(USA_2013_sts, control= list(
  range = 5:length(observed(USA_2013_sts)), k = kh$k, h= kh$h, trans = "rossi"))
plot(USA_2013_cusum, main = "USA HM CUSUM with Rossi method")

# Ecuador
Ecu_outbreak_cusum <- surveillance::cusum(Ecu_sts, control = list(
  range = 1:length(observed(Ecu_sts)), k = kh$k, h= kh$h, trans = "rossi", m = NULL
))
plot(Ecu_outbreak_cusum, main = "Ecuador HM CUSUM with Rossi method")

# for one season only
Ecu_2013_cusum <- surveillance::cusum(Ecu_2013_sts, control= list(
  range = 1:length(observed(Ecu_2013_sts)), k = kh$k, h= kh$h, trans = "rossi"))
plot(Ecu_2013_cusum, main = "Ecuador CUSUM with Rossi method")

# CUSUM does not work at all! For USA, the threshold is too high with all data and too low with only 2013 data
# Ecuador counts are too low for CUSUM to have any threshold apparently
# CUSUM can't handle 0 counts


## Bayes algorithms
USA_outbreak_b1 <- algo.bayes(USA_disProg, control = list(range = 10:length(observed(USA_sts)), alpha = 0.05, 
                                                          b = 0, w = 6))
plot(USA_outbreak_b1, main = "USA HM with bayes(6, 6, 0)")

Ecu_outbreak_b1 <- algo.bayes(Ecu_disProg, control = list(range = 10:length(observed(Ecu_sts)), alpha = 0.05, 
                                                          b = 0, w = 6))
plot(Ecu_outbreak_b1, main = "Ecuador HM with bayes(6, 6, 0)")
# ooes not look too bad for countries with sparse data, baseline length must be adjusted


## outbreakP algorithm (apparently especially suited for influenza outbreak detection)
# USA
USA_outbreakP <- outbreakP(USA_sts, control = list(
  range = 1:length(observed(USA_sts)), k = 100, ret = "cases"
))
plot(USA_outbreakP)
# not suitable for long time series with repeated outbreaks

USA_2013_outbreakP <- outbreakP(USA_2013_sts, control = list(
  range = 1:length(observed(USA_2013_sts)), k = 100, ret = "cases"
))
plot(USA_2013_outbreakP, main = "USA HM OutbreakP")
# too sensitive

# Ecuador

Ecu_outbreakP <- outbreakP(Ecu_sts, control = list(
  range = 1:length(observed(Ecu_sts)), k = 100, ret = "cases"
))
plot(Ecu_outbreakP)

Ecu_2013_outbreakP <- outbreakP(Ecu_2013_sts, control = list(
  range = 1:length(observed(Ecu_2013_sts)), k = 100, ret = "cases"
))
plot(Ecu_2013_outbreakP, main = "Ecuador HM OutbreakP")
# could be ok


## EWMA chart
USA_ewma <- ewma(HM_USA$counts, lambda = 0.3, nsigmas = 3,
                 title = "USA HM EWMA")
summary(USA_ewma)
# threshold for outbreak detection is too high even with HM data

# calculate mean and standard deviation during non-epidemic period (May-October for temperate regions in Northern hemisphere)
# baseline could also be all observations below the 30th percentile (or other value)
HM_USA_nep <- filter(HM_USA, month(HM_USA$date) %in% 5:10)
USA_meancount_nep <- mean(HM_USA_nep$counts)
USA_sdcount_nep <- sd(HM_USA_nep$counts)

# new EWMA chart with non-epidemic season mean and sd
USA_ewma_nep <- ewma(HM_USA$counts, lambda = 0.3, nsigmas = 3, center = USA_meancount_nep, std.dev = USA_sdcount_nep, 
                     title = "USA HM EWMA with baseline mean and SD")
summary(USA_ewma_nep)
# results look good when taking baseline mean and sd

# Ecuador
Ecu_ewma <- ewma(HM_Ecu$counts, lambda = 0.3, nsigmas = 3,
                 title = "Ecuador HM EWMA")
summary(Ecu_ewma)
# sparse data countries don't look too bad, because baseline is inherently low, but not too low


## AnomalyDetection library
# prepare the data
HM_USA_vec <- HM_USA$counts

anomalies_USA = AnomalyDetectionVec(HM_USA_vec, max_anoms = 0.2, direction = "pos", alpha = 0.05, plot = TRUE, 
                                period = 52, longterm_period = NULL)
anomalies_USA$plot
# this does not work really well, detects too few anomalies

#Ecuador
HM_Ecu_vec <- HM_Ecu$counts

anomalies_Ecu = AnomalyDetectionVec(HM_Ecu_vec, max_anoms = 0.1, direction = "pos", alpha = 0.05, plot = TRUE, 
                                    period = 52, longterm_period = NULL)
anomalies_Ecu$plot
# does not work at all, since every count is an anomaly now


##### change point analysis #####
## bcp
USA_bcp <- bcp(HM_USA$counts, burnin = 100, mcmc = 5000)
plot(USA_bcp)

HM_USA <- HM_USA %>% mutate(bcp.postprob = USA_bcp$posterior.prob)

first_cp <- HM_USA %>% group_by(season) %>% mutate(changepoint = ifelse(bcp.postprob > 0.5, date, NA)) %>% 
  filter(!is.na(changepoint)) %>% filter(date == min(date)) %>% ungroup() %>% select(date, changepoint)
HM_USA_cp <- left_join(HM_USA, first_cp, by = "date")

ggplot(data = HM_USA_cp, aes(x = date, y = counts)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "HealthMap data for USA with posterior probability of change point",
       caption = "Dashed vertical lines mark first change point per influenza season (July 1st - June 30th)") + 
  geom_vline(xintercept = na.omit(HM_USA_cp$changepoint), lty = 2, col = "red")


HM_USA_allcp <- HM_USA %>% mutate(changepoint = ifelse(bcp.postprob > 0.5, date, NA))
ggplot(data = HM_USA_allcp, aes(x = date, y = counts)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "HealthMap data for USA with posterior probability of change point",
       caption = "Dashed vertical lines mark all change points") + 
  geom_vline(xintercept = na.omit(HM_USA_allcp$changepoint), lty = 2, col = "red")


# with established outbreak criteria
HM_USA$bcp_criteria <- NA
for(i in 10:nrow(HM_USA)){
  if(HM_USA$bcp.postprob[i] >= 0.5 & HM_USA$bcp.postprob[(i-1)] < 0.5 & # transition from non-epidemic to epidemic
     mean(HM_USA$counts[(i-3):(i+3)]) > mean(HM_USA$counts[(i-4):i]) & # running mean to ensure that curve is rising
     sum(!is.na(HM_USA$bcp_criteria[(i-10):(i-1)])) == 0){ # No outbreak flagged during the previous 5 weeks
    HM_USA$bcp_criteria[i] <- HM_USA$date[i]
  } else {
    HM_USA$bcp_criteria[i] <- NA
  }
}
ggplot(data = HM_USA, aes(x = date, y = counts)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "Healthmap data for USA with posterior probability of change point") + 
  geom_vline(xintercept = na.omit(HM_USA$bcp_criteria), lty = 2, col = "red")


# Ecuador
Ecu_bcp <- bcp(HM_Ecu$counts, burnin = 100, mcmc = 5000)
plot(Ecu_bcp)

HM_Ecu <- HM_Ecu %>% mutate(bcp.postprob = Ecu_bcp$posterior.prob)

first_cp_Ecu <- HM_Ecu %>% group_by(season) %>% mutate(changepoint = ifelse(bcp.postprob > 0.5, date, NA)) %>% 
  filter(!is.na(changepoint)) %>% filter(date == min(date)) %>% ungroup() %>% select(date, changepoint)

HM_Ecu_cp <- left_join(HM_Ecu, first_cp_Ecu, by = "date")
ggplot(data = HM_Ecu_cp, aes(x = date, y = counts)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "HealthMap data for Ecuador with posterior probability of change point",
       caption = "Dashed vertical lines mark first change point per influenza season (July 1st - June 30th)") + 
  geom_vline(xintercept = na.omit(HM_Ecu_cp$changepoint), lty = 2, col = "red")

HM_Ecu_allcp <- HM_Ecu %>% mutate(changepoint = ifelse(bcp.postprob > 0.5, date, NA))
ggplot(data = HM_Ecu_allcp, aes(x = date, y = counts)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "HealthMap data for Ecuador with posterior probability of change point",
       caption = "Dashed vertical lines mark all change points") + 
  geom_vline(xintercept = na.omit(HM_Ecu_allcp$changepoint), lty = 2, col = "red")

## cpm
cpm_USA <- processStream(HM_USA$counts, cpmType = "Exponential", startup = 10, ARL0 = 500)
# exponential method works, Lepage, Mood, Kolmogorov-Smirnov, Cramer-von-Mises too give all very similar results except Mood)
HM_USA$cpm <- NA
HM_USA$cpm[cpm_USA$changePoints] <- HM_USA$date[cpm_USA$changePoints]

ggplot(data = HM_USA, aes(x = date, y = counts)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "Event counts", title = "HealthMap data for USA with change points, cpm package") + 
  geom_vline(xintercept = na.omit(HM_USA$cpm), lty = 2, col = "red")

# only CPs in rising curves
HM_USA$cpm_criteria <- NA
for(i in 4:nrow(HM_USA)){
  # use running mean of +/- 3 observations to ensure that we are in the inreasing part of the curve
  if(is.na(HM_USA$cpm[i]) == FALSE & mean(HM_USA$counts[(i-3):(i+3)]) > mean(HM_USA$counts[(i-4):(i+2)])){
    HM_USA$cpm_criteria[i] <- HM_USA$date[i]
  } else {
    HM_USA$cpm_criteria[i] <- NA
  }
}
ggplot(data = HM_USA, aes(x = date, y = counts)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "Event counts", title = "HealthMap data for USA with change points, cpm package") + 
  geom_vline(xintercept = na.omit(HM_USA$cpm_criteria), lty = 2, col = "red") + 
  geom_hline(yintercept = quantile(HM_USA$counts, 0.4))
 
#Ecuador
cpm_Ecu <- processStream(HM_Ecu$counts, cpmType = "Mann-Whitney", startup = 10, ARL0 = 500)
# Mann-Whitney works well, Lepage, Mood, Kolmogorov-Smirnov, Cramer-von-Mises all give weird regular signals which makes no sense at all
HM_Ecu$cpm <- NA
HM_Ecu$cpm[cpm_Ecu$changePoints] <- HM_Ecu$date[cpm_Ecu$changePoints]

ggplot(data = HM_Ecu, aes(x = date, y = counts)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "Event counts", title = "HealthMap data for Ecuador with change points, cpm package") + 
  geom_vline(xintercept = na.omit(HM_Ecu$cpm), lty = 2, col = "red")

# only CPs in rising curves
HM_Ecu$cpm_criteria <- NA
for(i in 6:nrow(HM_Ecu)){
  # use running mean of +/- 5 observations to ensure that we are in the inreasing part of the curve
  if(is.na(HM_Ecu$cpm[i]) == FALSE & mean(HM_Ecu$counts[(i-5):(i+5)]) > mean(HM_Ecu$counts[(i-6):(i+4)])){
    HM_Ecu$cpm_criteria[i] <- HM_Ecu$date[i]
  } else {
    HM_Ecu$cpm_criteria[i] <- NA
  }
}
ggplot(data = HM_Ecu, aes(x = date, y = counts)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "Event counts", title = "HealthMap data for Ecuador with change points, cpm package") + 
  geom_vline(xintercept = na.omit(HM_Ecu$cpm_criteria), lty = 2, col = "red")




##### Twitter BreakoutDetection#####
breakout_data_USA <- select(HM_USA, c(counts, date)) %>% rename(count = counts, timestamp = date)
breakout_USA <- breakout(breakout_data_USA, min.size = 5, method = "multi", plot = TRUE, beta = 0.001)
breakout_USA$plot


breakout_data_Ecu <- select(HM_Ecu, c(counts, date)) %>% rename(count = counts, timestamp = date)
breakout_Ecu <- breakout(breakout_data_Ecu, min.size = 5, method = "multi", plot = TRUE, beta = 0.001)
breakout_Ecu$plot
# does not work at all