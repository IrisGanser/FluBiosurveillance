library(plyr)
library(lubridate)
library(surveillance)
library(ggplot2)
library(MASS)
library(dplyr)
library(msm)
library(qcc)
library(AnomalyDetection)
library(bcp)
library(R2WinBUGS)
library(ecp)
library(cpm)
library(BreakoutDetection)


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

years <- 2012:2020
FluNet_data$season <- cut(FluNet_data$SDATE, 
                          breaks=as.POSIXct(paste(years,"-07-01",sep="")),
                          labels=paste(years[-length(years)],years[-length(years)]+1,sep="/"))



##### Periodic outbreak detection method #####

# look at France and at different pruning values
FluNet_France <- filter(FluNet_data, Country == "France") %>% select(c(SDATE, ALL_INF))
FluNet_France_Counts <- FluNet_France$ALL_INF
FluNet_France_Epoch <- as.Date(FluNet_France$SDATE)
prune_France <- quantile(FluNet_France$ALL_INF, probs = c(0.35, 0.4, 0.45, 0.5))


ggplot(data = FluNet_France, aes(x = SDATE, y = ALL_INF)) + 
  geom_line() +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for France") + 
  geom_hline(yintercept = prune_France, col = "red", lty = "dashed")


ggplot(data = FluNet_France, aes(x = SDATE, y = ALL_INF)) + 
  geom_line() +
  geom_point() + 
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for France (cut off so potential pruning limits can be shown)") + 
  geom_hline(yintercept = prune_France, col = "red", lty = "dashed") + 
  scale_y_continuous(limits = c(0, 100)) + 
  annotate("text", x = as.POSIXct("2013-01-01"), y = c(19.6, 27, 43.1, 77.6), 
           label = c("35th perc.", "40th perc.", "45th perc.", "50th perc."), col = "red", size = 3.4)
# 35th percentile cutoff might work best for France data, but is probably too low for other countries

ggplot(data = FluNet_France, aes(x = SDATE, y = ALL_INF)) + 
  geom_line() +
  geom_point() + 
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "laboratory-confirmed influenza case counts", title = "WHO FluNet data for France (zoom)") + 
  scale_y_continuous(limits = c(0, 100))

for (i in seq_along(country_list)) { 
  prune <- quantile(subset(FluNet_data$ALL_INF, FluNet_data$Country == country_list[i]), probs = c(0.35, 0.4, 0.45, 0.5))
  plot <- ggplot(data = subset(FluNet_data, FluNet_data$Country == country_list[i]), aes(x = SDATE, y = ALL_INF)) + 
    geom_line() +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "influenza case counts", title = paste("WHO FluNet data for", country_list[i], "(pruned)", sep=' ')) + 
    geom_hline(yintercept = prune, col = "red", lty = "dashed") + 
    scale_y_continuous(limits = c(0, prune[4] + 0.2*prune[4])) + 
    annotate("text", x = as.POSIXct("2013-01-01"), y = prune - 2.5, 
             label = c("35th perc.", "40th perc.", "45th perc.", "50th perc."), col = "red", size = 3.4)
  print(plot)
}

# problem 1: Some countries report counts only in 'epidemic season'
# problem 2: Epidemic threshold (pruning percentile) is different for every country, more of an eyeballing process


# filter for US, Nigeria and Ecuador as illustrative examples
FluNet_USA <- filter(FluNet_data, Country == "United States") %>% dplyr::select(c(SDATE, ALL_INF, season))
FluNet_USA_Counts <- FluNet_USA$ALL_INF
FluNet_USA_Epoch <- as.Date(FluNet_USA$SDATE)

FluNet_Nig <- filter(FluNet_data, Country == "Nigeria") %>% dplyr::select(c(SDATE, ALL_INF, season))
FluNet_Nig_Counts <- FluNet_Nig$ALL_INF
FluNet_Nig_Epoch <- as.Date(FluNet_Nig$SDATE)

FluNet_Ecu <- filter(FluNet_data, Country == "Ecuador") %>% dplyr::select(c(SDATE, ALL_INF, season))
FluNet_Ecu_Counts <- FluNet_Ecu$ALL_INF
FluNet_Ecu_Epoch <- as.Date(FluNet_Ecu$SDATE)


# construct sts objects
USA_sts <- sts(observed = FluNet_USA_Counts, epoch = FluNet_USA_Epoch, epochAsDate = TRUE)
plot(USA_sts)
head(epoch(USA_sts))
head(observed(USA_sts))

USA_disProg <- sts2disProg(USA_sts)


Nig_sts <- sts(observed = FluNet_Nig_Counts, epoch = FluNet_Nig_Epoch, epochAsDate = TRUE)
Nig_disProg <- sts2disProg(Nig_sts)
plot(Nig_sts)

Ecu_sts <- sts(observed = FluNet_Ecu_Counts, epoch = FluNet_Ecu_Epoch, epochAsDate = TRUE)
Ecu_disProg <- sts2disProg(Ecu_sts)
plot(Ecu_sts)

##### Percentage methods #####
### Neuzil method (weeks with > 1% of annual positive tests) and Izurieta method (weeks with > 5% of annual positive tests))

# years <- 2012:2020
# FluNet_USA$season <- cut(FluNet_USA$SDATE, 
#                          breaks=as.POSIXct(paste(years,"-07-01",sep="")),
#                          labels=paste(years[-length(years)],years[-length(years)]+1,sep="/"))
FluNet_USA <- FluNet_USA %>%  group_by(season) %>% mutate(annual_flu_pos = sum(ALL_INF)) %>% ungroup()
FluNet_USA <- FluNet_USA %>% mutate(perc_flu_pos = ALL_INF/annual_flu_pos, Neuzil = ifelse(perc_flu_pos > 0.01, TRUE, FALSE)) %>% 
  mutate(Neuzil_case_limit = annual_flu_pos*0.01)

ggplot(data = FluNet_USA, aes(x = SDATE, y = ALL_INF)) + 
  geom_line() +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with seasonal Neuzil outbreak limit", 
       caption = "Neuzil outbreak limit is weeks with > 1% of annual positive tests \n
       Influenza season: July 1st - June 30th") +
  geom_line(aes(x = SDATE, y = Neuzil_case_limit), col = "red", lty = 2, size = 1) 
# alerts very late, problem with missing data in non-epidemic season and multi-country scale

# Nigeria
FluNet_Nig <- FluNet_Nig %>%  group_by(season) %>% mutate(annual_flu_pos = sum(ALL_INF)) %>% ungroup()
FluNet_Nig <- FluNet_Nig %>% mutate(perc_flu_pos = ALL_INF/annual_flu_pos, Neuzil = ifelse(perc_flu_pos > 0.01, TRUE, FALSE)) %>% 
  mutate(Neuzil_case_limit = annual_flu_pos*0.01)

ggplot(data = FluNet_Nig, aes(x = SDATE, y = ALL_INF)) + 
  geom_line() +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Nigeria with seasonal Neuzil outbreak limit", 
       caption = "Neuzil outbreak limit is weeks with > 1% of annual positive tests \n
       Influenza season: July 1st - June 30th") +
  geom_line(aes(x = SDATE, y = Neuzil_case_limit), col = "red", lty = 2, size = 1) 

# Ecuador
FluNet_Ecu <- FluNet_Ecu %>%  group_by(season) %>% mutate(annual_flu_pos = sum(ALL_INF)) %>% ungroup()
FluNet_Ecu <- FluNet_Ecu %>% mutate(perc_flu_pos = ALL_INF/annual_flu_pos, Neuzil = ifelse(perc_flu_pos > 0.01, TRUE, FALSE)) %>% 
  mutate(Neuzil_case_limit = annual_flu_pos*0.01)

ggplot(data = FluNet_Ecu, aes(x = SDATE, y = ALL_INF)) + 
  geom_line() +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with seasonal Neuzil outbreak limit", 
       caption = "Neuzil outbreak limit is weeks with > 1% of annual positive tests \n
       Influenza season: July 1st - June 30th") +
  geom_line(aes(x = SDATE, y = Neuzil_case_limit), col = "red", lty = 2, size = 1) 


### proportion of positive isolates (Cowling method)
FluNet_USA <- filter(FluNet_data, Country == "United States") %>% mutate(Pos.Prop = ALL_INF/SPEC_PROCESSED_NB) %>%
  select(c(SDATE, ALL_INF, Pos.Prop, season))
FluNet_USA <- FluNet_USA %>% group_by(season) %>% mutate(max_Pos.Prop = max(Pos.Prop), threshold30 = 0.3*max_Pos.Prop, 
                                                         threshold20 = 0.2*max_Pos.Prop, threshold10 = 0.1*max_Pos.Prop) %>% 
  ungroup()

ggplot(data = FluNet_USA, aes(x = SDATE)) +
  geom_line(aes(y = Pos.Prop, col = "pos. sample prop")) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", 
       title = "WHO FluNet data for USA with % threshold of maximum positive sample proportion") +
  geom_line(aes(y = threshold30, col = "30% threshold"), lty = 2, size = 1) +
  geom_line(aes(y = threshold20, col = "20% threshold"), lty = 2, size = 1) +
  geom_line(aes(y = threshold10, col = "10% threshold"), lty = 2, size = 1) + 
  scale_color_brewer(palette = "Set1")

FluNet_Ecu <- filter(FluNet_data, Country == "Ecuador") %>% mutate(Pos.Prop = ALL_INF/SPEC_PROCESSED_NB) %>%
  select(c(SDATE, ALL_INF, Pos.Prop, season))
FluNet_Ecu <- FluNet_Ecu %>% group_by(season) %>% mutate(max_Pos.Prop = max(Pos.Prop), threshold30 = 0.3*max_Pos.Prop, 
                                                         threshold20 = 0.2*max_Pos.Prop, threshold10 = 0.1*max_Pos.Prop) %>% 
  ungroup()

ggplot(data = FluNet_Ecu, aes(x = SDATE)) +
  geom_line(aes(y = Pos.Prop, col = "pos. sample prop")) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", 
       title = "WHO FluNet data for Ecuador with % threshold of maximum positive sample proportion") +
  geom_line(aes(y = threshold30, col = "30% threshold"), lty = 2, size = 1) +
  geom_line(aes(y = threshold20, col = "20% threshold"), lty = 2, size = 1) +
  geom_line(aes(y = threshold10, col = "10% threshold"), lty = 2, size = 1) + 
  scale_color_brewer(palette = "Set1")


#### outbreak detection with surveillance pacakge #####
### Farrington algorithm
USA_outbreak_Farr <- farringtonFlexible(sts = USA_sts, control = list(
  b = 2, w = 3, weightsThreshold = 1, noPeriods = 1,
  pastWeeksNotIncluded = 26, pThresholdTrend = 0.05,
  thresholdMethod = "delta", trend = TRUE, thresholdMethod = "nbPlugin"
))
plot(USA_outbreak_Farr, main = "USA improved Farrington algorithm")
# not useful because takes seasonality into account

# Ecu_outbreak_Farr <- farringtonFlexible(sts = Ecu_sts, control = list(
#   b = 2, w = 3, weightsThreshold = 1, noPeriods = 1,
#   pastWeeksNotIncluded = 26, pThresholdTrend = 0.05,
#   thresholdMethod = "delta", trend = TRUE, thresholdMethod = "nbPlugin"
# ))
# plot(Ecu_outbreak_Farr, main = "USA improved Farrington algorithm")

## EARS C1-3 algorithms
# USA
USA_outbreak_C1 <- earsC(USA_sts, control = list(
  method = "C1", baseline = 7, minSigma = 1, alpha = 0.05
))
plot(USA_outbreak_C1, main = "USA EARS C1")
sum(alarms(USA_outbreak_C1))


USA_outbreak_C2_7 <- earsC(USA_sts, control = list(
  method = "C2", baseline = 7, minSigma = 1, alpha = 0.01
))
plot(USA_outbreak_C2_7, main = "USA EARS C2, 7 weeks baseline")

USA_outbreak_C2_14 <- earsC(USA_sts, control = list(
  method = "C2", baseline = 14, minSigma = 1, alpha = 0.01
))
plot(USA_outbreak_C2_14, main = "USA EARS C2, 14 weeks baseline")

USA_outbreak_C2_28 <- earsC(USA_sts, control = list(
  method = "C2", baseline = 28, minSigma = 1, alpha = 0.01
))
plot(USA_outbreak_C2_28, main = "USA EARS C2, 28 weeks baseline")

USA_outbreak_C3 <- earsC(USA_sts, control = list(
  method = "C3", baseline = 7, minSigma = 1, alpha = 0.01
))
plot(USA_outbreak_C3, main = "USA EARS 3")
sum(alarms(USA_outbreak_C3))
# tuning of parameters required, raises alarm very early in season, many false positives

# Ecuador
Ecu_outbreak_C2_7 <- earsC(Ecu_sts, control = list(
  method = "C2", baseline = 7, minSigma = 1, alpha = 0.01
))
plot(Ecu_outbreak_C2_7, main = "Ecuador EARS C2, 7 weeks baseline")

Ecu_outbreak_C2_14 <- earsC(Ecu_sts, control = list(
  method = "C2", baseline = 14, minSigma = 1, alpha = 0.01
))
plot(Ecu_outbreak_C2_14, main = "Ecuador EARS C2, 14 weeks baseline")

Ecu_outbreak_C2_28 <- earsC(Ecu_sts, control = list(
  method = "C2", baseline = 28, minSigma = 1, alpha = 0.01
))
plot(Ecu_outbreak_C2_28, main = "Ecuador EARS C2, 28 weeks baseline")



## CUSUM method
kh <- find.kh(ARLa = 500, ARLr = 7)

#USA
USA_outbreak_cusum <- surveillance::cusum(USA_sts, control = list(
  range = 1:length(observed(USA_sts)), k = kh$k, h= kh$h, trans = "rossi", m = NULL
))
plot(USA_outbreak_cusum, main = "USA CUSUM with Rossi method")

# for one season only
FluNet_USA_2013 <- filter(FluNet_USA, SDATE > "2013-04-01" & SDATE < "2014-01-01") 
USA_2013_sts <- sts(observed = FluNet_USA_2013$ALL_INF, epochAsDate = FALSE, start = c(2013, 4), frequency = 52)

USA_2013_cusum <- surveillance::cusum(USA_2013_sts, control= list(
  range = 5:length(observed(USA_2013_sts)), k = kh$k, h= kh$h, trans = "rossi"))
plot(USA_2013_cusum, main = "USA CUSUM with Rossi method")

#Ecuador
Ecu_outbreak_cusum <- surveillance::cusum(Ecu_sts, control = list(
  range = 1:length(observed(Ecu_sts)), k = kh$k, h= kh$h, trans = "rossi", m = NULL
))
plot(Ecu_outbreak_cusum, main = "Ecuador CUSUM with Rossi method")

# for one season only
FluNet_Ecu_2013 <- filter(FluNet_Ecu, SDATE > "2013-04-01" & SDATE < "2014-01-01") 
Ecu_2013_sts <- sts(observed = FluNet_Ecu_2013$ALL_INF, epochAsDate = FALSE, start = c(2013, 4), frequency = 52)

Ecu_2013_cusum <- surveillance::cusum(Ecu_2013_sts, control= list(
  range = 1:length(observed(Ecu_2013_sts)), k = kh$k, h= kh$h, trans = "rossi"))
plot(Ecu_2013_cusum, main = "Ecuador CUSUM with Rossi method")



# glrnb method
# estimate overdispersion parameter
USA_qp <- glm(formula = ALL_INF ~ 1, family = "quasipoisson", data = FluNet_USA)
summary(USA_qp)
# Overdispersion parameter is estimated as alpha = 8254.8, which is huge! mean is 7.95

USA_outbreak_glrnb <- glrnb(USA_sts, control = list(
  range = 30:length(observed(USA_sts)), c.ARL = 5, mu0 = NULL, alpha = FALSE, Mtilde=1, M=-1, change="intercept", 
  dir = "inc", theta=NULL, ret="value"
)) 
plot(USA_outbreak_glrnb)

USA_2013_glrnb <- glrnb(USA_2013_sts, control = list(
  range = 10:length(observed(USA_2013_sts)), c.ARL = 5, mu0 = NULL, alpha = FALSE, Mtilde=1, M=-1, change="intercept", 
  dir = "inc", theta=NULL, ret="value"
)) 
plot(USA_2013_glrnb)
# don't understand what's going on here, but it does not look very useful. 

# Hidden Markov model with the surveillance package
USA_outbreak_HMM <- algo.hmm(USA_disProg, control = list(
  range = 15:364, Mtilde = 2*52, noStates = 2, trend = FALSE, noHarmonics = 1
))
# hrows errors, I think I need a state indicator in disProj object

# Bayes model
USA_outbreak_b1 <- algo.bayes1(USA_disProg, control = list(range = 10:364, alpha = 0.05))
plot(USA_outbreak_b1)

USA_outbreak_b2 <- algo.bayes2(USA_disProg, control = list(range = 65:364, alpha = 0.05))
plot(USA_outbreak_b2)

USA_outbreak_b3 <- algo.bayes3(USA_disProg, control = list(range = 120:364, alpha = 0.05))
plot(USA_outbreak_b3)
# only b1 would be good here, b2 and b3 take into account data from previous years
# b1 is very sensitive and gives a lot of alarms


USA_outbreak_rki1 <- algo.rki1(USA_disProg, control = list(range = 8:364))
plot(USA_outbreak_rki1)
# only rki1 would be good here, rki2 and rki3 take into account data from previous years


# boda algorithm 
# control_boda <- list(range = 1:length(observed(USA_sts)), X = NULL, season = FALSE, prior = "iid", alpha = 0.05, samplingMethod = "marginals")
# USA_outbreak_boda <- boda(USA_sts, control = control_boda)
# plot(USA_outbreak_boda)
# took forever to calculate although the number of generated samples was very low



## outbreakP algorithm (apparently especially suited for influenza outbreak detection)
# USA
USA_outbreakP <- outbreakP(USA_sts, control = list(
  range = 1:length(observed(USA_sts)), k = 100, ret = "cases"
))
plot(USA_outbreakP)

USA_2013_outbreakP <- outbreakP(USA_2013_sts, control = list(
  range = 1:length(observed(USA_2013_sts)), k = 100, ret = "cases"
))
plot(USA_2013_outbreakP, main = "USA OutbreakP")

# Ecuador
Ecu_2013_outbreakP <- outbreakP(Ecu_2013_sts, control = list(
  range = 1:length(observed(Ecu_2013_sts)), k = 100, ret = "cases"
))
plot(Ecu_2013_outbreakP, main = "Ecuador OutbreakP")



## EWMA chart
USA_ewma <- ewma(FluNet_USA$ALL_INF, lambda = 0.3, nsigmas = 3,
                 title = "USA EWMA")
summary(USA_ewma)
# threshold for outbreak detection is way too high. 
# Somehow, the threshold must be brought down! Take only non-epidemic period for baseline calculation?

# calculate mean and standard deviation during non-epidemic period (May-October for temperate regions in Northern hemisphere)
FluNet_USA_nep <- filter(FluNet_USA, month(FluNet_USA$SDATE) %in% 5:10)
USA_meancount_nep <- mean(FluNet_USA_nep$ALL_INF)
USA_sdcount_nep <- sd(FluNet_USA_nep$ALL_INF)

# new EWMA chart with non-epidemic season mean and sd
USA_ewma_nep <- ewma(FluNet_USA$ALL_INF, lambda = 0.3, nsigmas = 3, center = USA_meancount_nep, std.dev = USA_sdcount_nep, 
                     title = "USA EWMA with baseline mean and SD")
summary(USA_ewma_nep)
# gives more reasonable and accurate results, tuning of lambda required, but 0.3 seems to be good

# Ecuador
Ecu_ewma <- ewma(FluNet_Ecu$ALL_INF, lambda = 0.3, nsigmas = 3,
                 title = "Ecuador EWMA")
summary(Ecu_ewma)






## AnomalyDetection library
# prepare the data
FluNet_USA_vec <- FluNet_USA$ALL_INF

anomalies = AnomalyDetectionVec(FluNet_USA_vec, max_anoms = 0.1, direction = "pos", alpha = 0.05, plot = TRUE, 
                               period = 52, longterm_period = NULL)
anomalies$plot
# this does not work at all

# try anomaly detection with only one season of data
FluNet_USA_2013 <- filter(FluNet_USA, SDATE > "2013-04-01" & SDATE < "2014-01-01") %>% select(c(SDATE, ALL_INF))
anomalies_2013 = AnomalyDetectionTs(FluNet_USA_2013, max_anoms = 0.1, direction = "pos", alpha = 0.05, plot = TRUE,
                                    longterm = FALSE)
anomalies_2013$plot
# nope! FluNet data are not suited for this algorithm, might be worth a try with HM and EIOS data




## Bayesian change point analysis
USA_bcp <- bcp(FluNet_USA$ALL_INF, burnin = 100, mcmc = 5000)
plot(USA_bcp)

USA_bcp_2013 <- bcp(FluNet_USA_2013$ALL_INF, burnin = 100, mcmc = 5000)
plot(USA_bcp_2013)


FluNet_USA <- FluNet_USA %>% mutate(bcp.postprob = USA_bcp$posterior.prob)

first_cp <- FluNet_USA %>% group_by(season) %>% mutate(changepoint = ifelse(bcp.postprob > 0.5, SDATE, NA)) %>% 
  filter(!is.na(changepoint)) %>% filter(SDATE == min(SDATE)) %>% ungroup() %>% dplyr::select(SDATE, changepoint)

FluNet_USA_cp <- left_join(FluNet_USA, first_cp, by = "SDATE")

# Ctrl-Shift-C for commenting whole sections of code!
# for(i in levels(FluNet_USA2$season)){
#   for(j in 1:length(FluNet_USA2$season[FluNet_USA2$season == i])){
#     index = j
#     print(index)
#     if(FluNet_USA2$bcp.postprob[j] >= 0.5){
#       FluNet_USA2$changepoint[j] <- as.POSIXct(FluNet_USA2$SDATE[j])
#     #   break
#     # } else {
#     #   FluNet_USA2$changepoint[j] <- 1
#     }
#   }
# }
# 
# first_cp <- c(2, 48, 98, 160, 206, 257, 310, 359)
# FluNet_USA$changepoint <- NA
# FluNet_USA$changepoint[first_cp] <- FluNet_USA$SDATE[first_cp]


ggplot(data = FluNet_USA_cp, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with posterior probability of change point",
       caption = "Dashed vertical lines mark first change point per influenza season (July 1st - June 30th)") + 
  geom_vline(xintercept = na.omit(FluNet_USA_cp$changepoint), lty = 2, col = "red")


FluNet_USA_allcp <- FluNet_USA %>% mutate(changepoint = ifelse(bcp.postprob > 0.5, SDATE, NA))
ggplot(data = FluNet_USA_allcp, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with posterior probability of change point") + 
  geom_vline(xintercept = na.omit(FluNet_USA_allcp$changepoint), lty = 2, col = "red")


# cp determined by special criteria

FluNet_USA$bcp_criteria <- NA

for(i in 10:nrow(FluNet_USA)){
  if(FluNet_USA$bcp.postprob[i] >= 0.5 & FluNet_USA$bcp.postprob[(i-1)] < 0.5 & # transition from non-epidemic to epidemic
     mean(FluNet_USA$ALL_INF[(i-5):(i+5)]) > mean(FluNet_USA$ALL_INF [(i-6):(i+4)]) & # running mean to ensure that curve is rising
     sum(!is.na(FluNet_USA$bcp_criteria[(i-10):(i-1)])) == 0){ # No outbreak flagged during the previous 5 weeks
    FluNet_USA$bcp_criteria[i] <- FluNet_USA$SDATE[i]
  } else {
    FluNet_USA$bcp_criteria[i] <- NA
  }
}
ggplot(data = FluNet_USA, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with posterior probability of change point") + 
  geom_vline(xintercept = na.omit(FluNet_USA$bcp_criteria), lty = 2, col = "red")

#Ecuador
Ecu_bcp <- bcp(FluNet_Ecu$ALL_INF, burnin = 100, mcmc = 5000)
plot(Ecu_bcp)

FluNet_Ecu <- FluNet_Ecu %>% mutate(bcp.postprob = Ecu_bcp$posterior.prob)

first_cp_Ecu <- FluNet_Ecu %>% group_by(season) %>% mutate(changepoint = ifelse(bcp.postprob > 0.5, SDATE, NA)) %>% 
  filter(!is.na(changepoint)) %>% filter(SDATE == min(SDATE)) %>% ungroup() %>% select(SDATE, changepoint)

FluNet_Ecu_cp <- left_join(FluNet_Ecu, first_cp_Ecu, by = "SDATE")
ggplot(data = FluNet_Ecu_cp, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with posterior probability of change point",
       caption = "Dashed vertical lines mark first change point per influenza season (July 1st - June 30th)") + 
  geom_vline(xintercept = na.omit(FluNet_Ecu_cp$changepoint), lty = 2, col = "red")

FluNet_Ecu_allcp <- FluNet_Ecu %>% mutate(changepoint = ifelse(bcp.postprob > 0.5, SDATE, NA))
ggplot(data = FluNet_Ecu_allcp, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with posterior probability of change point",
       caption = "Dashed vertical lines mark all change points") + 
  geom_vline(xintercept = na.omit(FluNet_Ecu_allcp$changepoint), lty = 2, col = "red")


FluNet_Ecu$bcp_criteria <- NA

for(i in 10:nrow(FluNet_Ecu)){
  if(FluNet_Ecu$bcp.postprob[i] >= 0.5 & FluNet_Ecu$bcp.postprob[i-1] < 0.5 & 
     mean(FluNet_Ecu$ALL_INF[(i-5):(i+5)]) > mean(FluNet_Ecu$ALL_INF [(i-6):(i+4)]) & 
     sum(!is.na(FluNet_Ecu$bcp_criteria[(i-10):(i-1)])) == 0){
    FluNet_Ecu$bcp_criteria[i] <- FluNet_Ecu$SDATE[i]
  } else {
    FluNet_Ecu$bcp_criteria[i] <- NA
  }
}
ggplot(data = FluNet_Ecu, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with posterior probability of change point") + 
  geom_vline(xintercept = na.omit(FluNet_Ecu$bcp_criteria), lty = 2, col = "red")



##### apply bcp to all countries #####
country_list <- levels(FluNet_data$Country)
country_code <- c("Arg", "Aus", "Bra", "Bul", "Chn", "Cri", "Ecu", "Egy", "Fra", "Gbr", "Ger", "Grc",
                  "Ind", "Irn", "Mex", "Nig", "Rus", "Sau", "Swe", "Tha", "Ury", "Usa", "Vnm", "Zaf")



for (i in seq_along(country_list)) { 
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==country_list[i])
  
  bcp_temp <- bcp(FluNet_temp$ALL_INF, burnin = 100, mcmc = 5000)
  plot(bcp_temp)
  
  FluNet_data$bcp.postprob[FluNet_data$Country==country_list[i]] <- bcp_temp$posterior.prob
}
FluNet_data <- FluNet_data %>% mutate(changepoint = ifelse(bcp.postprob > 0.5, SDATE, NA))

for (i in seq_along(country_list)) {
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==country_list[i])
  
  plot <- ggplot(FluNet_temp, aes(x = SDATE, y = ALL_INF)) + 
    geom_line(aes(col = bcp.postprob), size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "influenza case counts", 
         title = paste("WHO FluNet data for", country_list[i], "with posterior probability of change point", sep = " "),
         caption = "Dashed vertical lines mark all change points") + 
    geom_vline(xintercept = na.omit(FluNet_temp$changepoint), lty = 2, col = "red")
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet bcp", country_list[i], ".jpeg", sep=' '))
}



##### ecp package #####
obsUSA <- as.matrix(FluNet_USA$ALL_INF)
ecp_USA <- e.divisive(obsUSA, R = 499, k = NULL, min.size = 5, alpha = 1)
ecp_USA$estimates
cp_USA <- ecp_USA$estimates[c(-1, -length(ecp_USA$estimates))]

FluNet_USA$ecp <- NA
FluNet_USA$ecp[cp_USA] <- FluNet_USA$SDATE[cp_USA]

ggplot(data = FluNet_USA, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with change points") + 
  geom_vline(xintercept = na.omit(FluNet_USA$ecp), lty = 2, col = "red")

# with criteria
FluNet_USA$ecp_criteria <- NA
for(i in 10:nrow(FluNet_USA)){
  if(is.na(FluNet_USA$ecp[i]) == FALSE & 
     mean(FluNet_USA$ALL_INF[(i-5):(i+5)]) > mean(FluNet_USA$ALL_INF [(i-6):(i+4)]) &
     sum(!is.na(FluNet_USA$ecp_criteria[(i-10):(i-1)])) == 0){
    FluNet_USA$ecp_criteria[i] <- FluNet_USA$SDATE[i]
  } else {
    FluNet_USA$ecp_criteria[i] <- NA
  }
}
ggplot(data = FluNet_USA, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with change points, ecp package") + 
  geom_vline(xintercept = na.omit(FluNet_USA$ecp_criteria), lty = 2, col = "red")
# ecp alerts very late

# Ecuador

obsEcu <- as.matrix(FluNet_Ecu$ALL_INF)
ecp_Ecu <- e.divisive(obsEcu, R = 499, k = NULL, min.size = 5, alpha = 1)
ecp_Ecu$estimates
cp_Ecu <- ecp_Ecu$estimates[c(-1, -length(ecp_Ecu$estimates))]

FluNet_Ecu$ecp <- NA
FluNet_Ecu$ecp[cp_Ecu] <- FluNet_Ecu$SDATE[cp_Ecu]

ggplot(data = FluNet_Ecu, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with change points, ecp package") + 
  geom_vline(xintercept = na.omit(FluNet_Ecu$ecp), lty = 2, col = "red")

# with criteria
FluNet_Ecu$ecp_criteria <- NA
for(i in 10:nrow(FluNet_Ecu)){
  if(is.na(FluNet_Ecu$ecp[i]) == FALSE & 
     mean(FluNet_Ecu$ALL_INF[(i-5):(i+5)]) > mean(FluNet_Ecu$ALL_INF [(i-6):(i+4)]) &
     sum(!is.na(FluNet_Ecu$ecp_criteria[(i-10):(i-1)])) == 0){
    FluNet_Ecu$ecp_criteria[i] <- FluNet_Ecu$SDATE[i]
  } else {
    FluNet_Ecu$ecp_criteria[i] <- NA
  }
}
ggplot(data = FluNet_Ecu, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with change points, ecp package") + 
  geom_vline(xintercept = na.omit(FluNet_Ecu$ecp_criteria), lty = 2, col = "red")


## cpm package
cpm_USA <- processStream(FluNet_USA$ALL_INF, cpmType = "Exponential", startup = 10, ARL0 = 500)
# exponential method clearly worked best
FluNet_USA$cpm <- NA
FluNet_USA$cpm[cpm_USA$changePoints] <- FluNet_USA$SDATE[cpm_USA$changePoints]

ggplot(data = FluNet_USA, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with change points, cpm package, 'Exponential' algorithm") + 
  geom_vline(xintercept = na.omit(FluNet_USA$cpm), lty = 2, col = "red")

FluNet_USA$cpm_criteria <- NA
for(i in 10:nrow(FluNet_USA)){
  if(is.na(FluNet_USA$cpm[i]) == FALSE & 
     mean(FluNet_USA$ALL_INF[(i-5):(i+5)]) > mean(FluNet_USA$ALL_INF[(i-6):(i+4)]) &
     sum(!is.na(FluNet_USA$cpm_criteria[(i-10):(i-1)])) == 0){
    FluNet_USA$cpm_criteria[i] <- FluNet_USA$SDATE[i]
  } else {
    FluNet_USA$cpm_criteria[i] <- NA
  }
}
ggplot(data = FluNet_USA, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with change points, cpm package, 'Exponential' algorithm") + 
  geom_vline(xintercept = na.omit(FluNet_USA$cpm_criteria), lty = 2, col = "red")


# Ecuador
cpm_Ecu <- processStream(FluNet_Ecu$ALL_INF, cpmType = "Exponential", startup = 10, ARL0 = 500)
# exponential method clearly worked best
FluNet_Ecu$cpm <- NA
FluNet_Ecu$cpm[cpm_Ecu$changePoints] <- FluNet_Ecu$SDATE[cpm_Ecu$changePoints]

ggplot(data = FluNet_Ecu, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with change points, cpm package, 'Exponential' algorithm") + 
  geom_vline(xintercept = na.omit(FluNet_Ecu$cpm), lty = 2, col = "red")

FluNet_Ecu$cpm_criteria <- NA
for(i in 10:nrow(FluNet_Ecu)){
  if(is.na(FluNet_Ecu$cpm[i]) == FALSE & 
     mean(FluNet_Ecu$ALL_INF[(i-5):(i+5)]) > mean(FluNet_Ecu$ALL_INF[(i-6):(i+4)]) &
     sum(!is.na(FluNet_Ecu$cpm_criteria[(i-10):(i-1)])) == 0){
    FluNet_Ecu$cpm_criteria[i] <- FluNet_Ecu$SDATE[i]
  } else {
    FluNet_Ecu$cpm_criteria[i] <- NA
  }
}
ggplot(data = FluNet_Ecu, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with change points, cpm package, 'Exponential' algorithm") + 
  geom_vline(xintercept = na.omit(FluNet_Ecu$cpm_criteria), lty = 2, col = "red")


### summary of the bcp, ecp, and cpm packages
# USA
ggplot(data = FluNet_USA, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with change points", color = "R package") + 
  geom_vline(aes(xintercept = FluNet_USA$bcp_criteria, col = "BCP"), lty = 2, size = 0.75) + 
  geom_vline(aes(xintercept = FluNet_USA$ecp_criteria, col = "ECP"), lty = 2, size = 0.75) + 
  geom_vline(aes(xintercept = FluNet_USA$cpm_criteria, col = "CPM"), lty = 2, size = 0.75) + 
  scale_color_manual(values = c("#cf0007", "#ff7903", "#006100")) + 
  theme(legend.key.size = unit(1, "cm"))

# Ecuador
ggplot(data = FluNet_Ecu, aes(x = SDATE, y = ALL_INF)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for Ecuador with change points", color = "R package") + 
  geom_vline(aes(xintercept = FluNet_Ecu$bcp_criteria, col = "BCP"), lty = 2, size = 0.75) + 
  geom_vline(aes(xintercept = FluNet_Ecu$ecp_criteria, col = "ECP"), lty = 2, size = 0.75) + 
  geom_vline(aes(xintercept = FluNet_Ecu$cpm_criteria, col = "CPM"), lty = 2, size = 0.75) + 
  scale_color_manual(values = c("#cf0007", "#ff7903", "#006100")) + 
  theme(legend.key.size = unit(1, "cm"))



##### Twitter BreakoutDetection#####
breakout_data_USA <- select(FluNet_USA, c(ALL_INF, SDATE)) %>% rename(count = ALL_INF, timestamp = SDATE)
breakout_USA <- breakout(breakout_data_USA, min.size = 10, method = "multi", plot = TRUE)
breakout_USA$plot
# does not work at all, maybe better with HM or EIOS


# potentially problematic countries because of data quality
problematic_FluNet <- c("Brazil", "Costa Rica", "Egypt", "France", "India", "Nigeria", "Saudi Arabia", "Thailand", "Uruguay", "Vietnam")
#'%!in%' <- function(x,y)!('%in%'(x,y))
nonproblematic_FluNet <- subset(levels(FluNet_data$Country), !(levels(FluNet_data$Country) %in% problematic_FluNet))


# apply bcp (with outbreak criteria) to all non-problematic countries 
for (i in seq_along(nonproblematic_FluNet)) { 
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==nonproblematic_FluNet[i])
  
  bcp_temp <- bcp(FluNet_temp$ALL_INF, burnin = 100, mcmc = 5000)
  plot(bcp_temp)
  
  
  FluNet_data$bcp.postprob[FluNet_data$Country==nonproblematic_FluNet[i]] <- bcp_temp$posterior.prob
}

for (i in seq_along(nonproblematic_FluNet)) { 
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==nonproblematic_FluNet[i] & is.na(FluNet_data$bcp.postprob) == FALSE)
  
  FluNet_temp$bcp_criteria <- NA
  for(j in 10:nrow(FluNet_temp)){
    if(FluNet_temp$bcp.postprob[j] >= 0.5 & FluNet_temp$bcp.postprob[(j-1)] < 0.5 & # transition from non-epidemic to epidemic
       mean(FluNet_temp$ALL_INF[(j-3):(j+3)], na.rm = TRUE) > mean(FluNet_temp$ALL_INF[(j-4):j], na.rm = TRUE) & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(FluNet_temp$bcp_criteria[(j-10):(j-1)])) == 0){ # No outbreak flagged during the previous 10 weeks
      FluNet_temp$bcp_criteria[j] <- FluNet_temp$SDATE[j]
    } else {
      FluNet_temp$bcp_criteria[j] <- NA
    }
  }

  plot <- ggplot(FluNet_temp, aes(x = SDATE, y = ALL_INF)) + 
    geom_line(aes(col = bcp.postprob), size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "influenza case counts", 
         title = paste("WHO FluNet data for", nonproblematic_FluNet[i], "with posterior probability of change point", sep = " ")) + 
    geom_vline(xintercept = na.omit(FluNet_temp$bcp_criteria), lty = 2, col = "red")
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet bcp", country_list[i], ".jpeg", sep=' '))
}



# do the same with problematic countries and look at the resuls
for (i in seq_along(problematic_FluNet)) { 
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==problematic_FluNet[i])
  
  bcp_temp <- bcp(FluNet_temp$ALL_INF, burnin = 100, mcmc = 5000)
  plot(bcp_temp)
  
  
  FluNet_data$bcp.postprob[FluNet_data$Country==problematic_FluNet[i]] <- bcp_temp$posterior.prob
}

for (i in seq_along(problematic_FluNet)) { 
  FluNet_temp <- filter(FluNet_data, FluNet_data$Country==problematic_FluNet[i] & is.na(FluNet_data$bcp.postprob) == FALSE)
  
  FluNet_temp$bcp_criteria <- NA
  for(j in 10:nrow(FluNet_temp)){
    if(FluNet_temp$bcp.postprob[j] >= 0.5 & FluNet_temp$bcp.postprob[(j-1)] < 0.5 & # transition from non-epidemic to epidemic
       mean(FluNet_temp$ALL_INF[(j-3):(j+3)], na.rm = TRUE) > mean(FluNet_temp$ALL_INF[(j-4):j], na.rm = TRUE) & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(FluNet_temp$bcp_criteria[(j-10):(j-1)])) == 0){ # No outbreak flagged during the previous 10 weeks
      FluNet_temp$bcp_criteria[j] <- FluNet_temp$SDATE[j]
    } else {
      FluNet_temp$bcp_criteria[j] <- NA
    }
  }
  
  plot <- ggplot(FluNet_temp, aes(x = SDATE, y = ALL_INF)) + 
    geom_line(aes(col = bcp.postprob), size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "influenza case counts", 
         title = paste("WHO FluNet data for", problematic_FluNet[i], "with posterior probability of change point", sep = " ")) + 
    geom_vline(xintercept = na.omit(FluNet_temp$bcp_criteria), lty = 2, col = "red")
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet bcp", country_list[i], ".jpeg", sep=' '))
}

                                