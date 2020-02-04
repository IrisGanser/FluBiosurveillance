library(plyr)
library(lubridate)
library(surveillance)
library(ggplot2)
library(MASS)
library(dplyr)
library(msm)

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
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for France (cut off so potential pruning limits can be shown)") + 
  geom_hline(yintercept = prune_France, col = "red", lty = "dashed") + 
  scale_y_continuous(limits = c(0, 100)) + 
  annotate("text", x = as.POSIXct("2013-01-01"), y = c(19.6, 27, 43.1, 77.6), 
           label = c("35th perc.", "40th perc.", "45th perc.", "50th perc."), col = "red", size = 3.4)
# 35th percentile cutoff might work best for France data, but is probably too low for other countries

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



#### outbreak detection with surveillance pacakge #####
# filter for US as illustrative example
FluNet_USA <- filter(FluNet_data, Country == "United States") %>% select(c(SDATE, ALL_INF))
FluNet_USA_Counts <- FluNet_USA$ALL_INF
FluNet_USA_Epoch <- as.Date(FluNet_USA$SDATE)

# construct sts object
USA_sts <- sts(observed = FluNet_USA_Counts, epoch = FluNet_USA_Epoch, epochAsDate = TRUE)
plot(USA_sts)
head(epoch(USA_sts))
head(observed(USA_sts))

USA_disProg <- sts2disProg(USA_sts)

# try Farrington algorithm
USA_outbreak_Farr <- farringtonFlexible(sts = USA_sts, control = list(
  b = 2, w = 3, weightsThreshold = 1, noPeriods = 1,
  pastWeeksNotIncluded = 26, pThresholdTrend = 0.05,
  thresholdMethod = "delta", trend = TRUE, thresholdMethod = "nbPlugin"
))

plot(USA_outbreak_Farr)
# not useful because takes seasonality into account


# try EARS C1-3 algorithms
USA_outbreak_C1 <- earsC(USA_sts, control = list(
  method = "C1", baseline = 7, minSigma = 1, alpha = 0.05
))
plot(USA_outbreak_C1)
sum(alarms(USA_outbreak_C1))

USA_outbreak_C2 <- earsC(USA_sts, control = list(
  method = "C2", baseline = 7, minSigma = 1, alpha = 0.01
))
plot(USA_outbreak_C2)
sum(alarms(USA_outbreak_C2))

USA_outbreak_C3 <- earsC(USA_sts, control = list(
  method = "C3", baseline = 7, minSigma = 1, alpha = 0.01
))
plot(USA_outbreak_C3)
sum(alarms(USA_outbreak_C3))
# works fine, tuning of parameters required, raises alarm very early in season


# cusum method
USA_outbreak_cusum <- cusum(USA_sts, control = list(
  range = 1:length(observed(USA_sts)), k = 1.04, h = 2.26, trans = "standard", m = NULL
))
plot(USA_outbreak_cusum)
# I don't really understand what's going on here, which parameters to tune and why the thresholds look like this


# glrnb method
# estimate overdispersion parameter
USA_qp <- glm(formula = ALL_INF ~ 1, family = "quasipoisson", data = FluNet_USA)
summary(USA_qp)
# Overdispersion parameter is estimated as alpha = 8254.8, which is huge! mean is 7.95

USA_outbreak_glrnb <- glrnb(USA_sts, control = list(
  range = 1:length(observed(USA_sts)), c.ARL = 5, mu0 = NULL, alpha = 8254.8, Mtilde=1, M=-1, change="intercept", 
  dir = "inc", theta=NULL, ret="value"
)) 
plot(USA_outbreak_glrnb)
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

control_boda <- list(range = 1:length(observed(USA_sts)), X = NULL, season = FALSE, prior = "iid", alpha = 0.05, samplingMethod = "marginals")
# USA_outbreak_boda <- boda(USA_sts, control = control_boda)
plot(USA_outbreak_boda)
#took forever to calculate although the number of generated samples was very low


# Neuzil method (weeks with > 1% of annual positive tests) and Izurieta method (weeks with > 5% of annual positive tests))
years <- 2012:2020
FluNet_USA$season <- cut(FluNet_USA$SDATE, 
                         breaks=as.POSIXct(paste(years,"-07-01",sep="")),
                         labels=paste(years[-length(years)],years[-length(years)]+1,sep="/"))
FluNet_USA <- FluNet_USA %>%  group_by(season) %>% mutate(annual_flu_pos = sum(ALL_INF)) %>% ungroup()
FluNet_USA <- FluNet_USA %>% mutate(perc_flu_pos = ALL_INF/annual_flu_pos, Neuzil = ifelse(perc_flu_pos > 0.01, TRUE, FALSE), 
                                    Izurieta = ifelse(perc_flu_pos > 0.05, TRUE, FALSE)) 
FluNet_USA <- FluNet_USA %>% mutate(Neuzil_case_limit = annual_flu_pos*0.01, Izurieta_case_limit = annual_flu_pos * 0.05)

ggplot(data = FluNet_USA, aes(x = SDATE, y = ALL_INF)) + 
  geom_line() +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet data for USA with Neuzil outbreak limit", 
       caption = "Neuzil outbreak limit is weeks with > 1% of annual positive tests") +
  geom_line(aes(x = SDATE, y = Neuzil_case_limit), col = "red", lty = 2) 

FluNet_data$season <- cut(FluNet_data$SDATE, 
                          breaks=as.POSIXct(paste(years,"-07-01",sep="")),
                          labels=paste(years[-length(years)],years[-length(years)]+1,sep="/"))
# alerts very late, problem with missing data in non-epidemic season and multi-country scale
