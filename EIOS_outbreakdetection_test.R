library(ggplot2)
library(plyr)
library(readxl)
library(lubridate)
library(stringr)
library(dplyr)
library(bcp)
library(qcc)
library(surveillance)
library(AnomalyDetection)
library(BreakoutDetection)
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
EIOS_byweek <- EIOSreports %>% group_by(date = floor_date(importDate, "week"), country = Country, .drop = FALSE) %>%
  summarize(counts=n()) %>% as.data.frame()


# filter for USA and Ecuador
EIOS_USA <- filter(EIOS_byweek, country == "United States")
EIOS_Ecu <- filter(EIOS_byweek, country == "Ecuador")

# construct sts objects
EIOS_USA_Counts <- EIOS_USA$counts
EIOS_USA_Epoch <- as.Date(EIOS_USA$date)
USA_sts <- sts(observed = EIOS_USA_Counts, epoch = EIOS_USA_Epoch, epochAsDate = TRUE)

EIOS_Ecu_Counts <- EIOS_Ecu$counts
EIOS_Ecu_Epoch <- as.Date(EIOS_Ecu$date)
Ecu_sts <- sts(observed = EIOS_Ecu_Counts, epoch = EIOS_Ecu_Epoch, epochAsDate = TRUE)

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
# 7 weeks baseline looks good, 14 too long

# Ecuador
Ecu_outbreak_C2_7 <- earsC(Ecu_sts, control = list(
  method = "C2", baseline = 7, minSigma = 1, alpha = 0.01
))
plot(Ecu_outbreak_C2_7, main = "Ecuador HM EARS C2, 7 weeks baseline")

Ecu_outbreak_C2_14 <- earsC(Ecu_sts, control = list(
  method = "C2", baseline = 14, minSigma = 1, alpha = 0.01
))
plot(Ecu_outbreak_C2_14, main = "Ecuador HM EARS C2, 14 weeks baseline")


## CUSUM method
kh <- find.kh(ARLa = 500, ARLr = 7)

#USA
USA_outbreak_cusum <- surveillance::cusum(USA_sts, control = list(
  range = 1:length(observed(USA_sts)), k = kh$k, h= kh$h, trans = "rossi", m = NULL
))
plot(USA_outbreak_cusum, main = "USA HM CUSUM with Rossi method")

# Ecuador
Ecu_outbreak_cusum <- surveillance::cusum(Ecu_sts, control = list(
  range = 1:length(observed(Ecu_sts)), k = kh$k, h= kh$h, trans = "rossi", m = NULL
))
plot(Ecu_outbreak_cusum, main = "Ecuador HM CUSUM with Rossi method")


## Bayes algorithms
USA_outbreak_b1 <- algo.bayes(USA_disProg, control = list(range = 7:length(observed(USA_sts)), alpha = 0.05, 
                                                          b = 0, w = 6))
plot(USA_outbreak_b1, main = "USA HM with bayes(6, 6, 0)")

Ecu_outbreak_b1 <- algo.bayes(Ecu_disProg, control = list(range = 7:length(observed(Ecu_sts)), alpha = 0.05, 
                                                          b = 0, w = 6))
plot(Ecu_outbreak_b1, main = "Ecuador HM with bayes(6, 6, 0)")
# too sensitive for both countries, might be still acceptable for Ecuador


## outbreakP algorithm (apparently especially suited for influenza outbreak detection)
# USA
USA_outbreakP <- outbreakP(USA_sts, control = list(
  range = 1:length(observed(USA_sts)), k = 100, ret = "cases"
))
plot(USA_outbreakP)
# not suitable for long time series with repeated outbreaks

Ecu_outbreakP <- outbreakP(Ecu_sts, control = list(
  range = 1:length(observed(Ecu_sts)), k = 100, ret = "cases"
))
plot(Ecu_outbreakP)


## EWMA chart
USA_ewma <- ewma(EIOS_USA$counts, lambda = 0.3, nsigmas = 3,
                 title = "USA HM EWMA")
summary(USA_ewma)
# threshold for outbreak detection is too high even with EIOS data

# calculate mean and standard deviation during non-epidemic period (May-October for temperate regions in Northern hemisphere)
# baseline could also be all observations below the 30th percentile (or other value)
EIOS_USA_nep <- filter(EIOS_USA, month(EIOS_USA$date) %in% 5:10)
USA_meancount_nep <- mean(EIOS_USA_nep$counts)
USA_sdcount_nep <- sd(EIOS_USA_nep$counts)

# new EWMA chart with non-epidemic season mean and sd
USA_ewma_nep <- ewma(EIOS_USA$counts, lambda = 0.3, nsigmas = 3, center = USA_meancount_nep, std.dev = USA_sdcount_nep, 
                     title = "USA HM EWMA with baseline mean and SD")
summary(USA_ewma_nep)
# results are now not sensitive enough because of some high counts in baseline calculation

# Ecuador
Ecu_ewma <- ewma(EIOS_Ecu$counts, lambda = 0.3, nsigmas = 3,
                 title = "Ecuador HM EWMA")
summary(Ecu_ewma)
# not very sensitive, but looks not too bad



##### change point analysis #####
## bcp
USA_bcp <- bcp(EIOS_USA$counts, burnin = 100, mcmc = 5000)
plot(USA_bcp)

EIOS_USA <- EIOS_USA %>% mutate(bcp.postprob = USA_bcp$posterior.prob)

EIOS_USA_allcp <- EIOS_USA %>% mutate(all_cp = ifelse(bcp.postprob > 0.5, date, NA))
ggplot(data = EIOS_USA_allcp, aes(x = date, y = counts)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "EIOS data for USA with posterior probability of change point",
       caption = "Dashed vertical lines mark all change points") + 
  geom_vline(xintercept = na.omit(EIOS_USA_allcp$all_cp), lty = 2, col = "red")

# cp determined by special criteria
EIOS_USA$bcp_criteria <- NA

for(i in 4:nrow(EIOS_USA)){
  if(EIOS_USA$bcp.postprob[i] >= 0.5 & EIOS_USA$bcp.postprob[(i-1)] < 0.5 & # transition from non-epidemic to epidemic
     mean(EIOS_USA$counts[(i-2):(i+2)]) > mean(EIOS_USA$counts[(i-3):i]) & # running mean to ensure that curve is rising
     sum(!is.na(EIOS_USA$bcp_criteria[(i-4):(i-1)])) == 0){ # No outbreak flagged during the previous 6 weeks
    EIOS_USA$bcp_criteria[i] <- EIOS_USA$date[i]
  } else {
    EIOS_USA$bcp_criteria[i] <- NA
  }
}
ggplot(data = EIOS_USA, aes(x = date, y = counts)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "4 months") + 
  labs(x = "", y = "influenza case counts", title = "EIOS data for USA with posterior probability of change point") + 
  geom_vline(xintercept = na.omit(EIOS_USA$bcp_criteria), lty = 2, col = "red")

# Ecuador
Ecu_bcp <- bcp(EIOS_Ecu$counts, burnin = 100, mcmc = 5000)
plot(Ecu_bcp)

EIOS_Ecu <- EIOS_Ecu %>% mutate(bcp.postprob = Ecu_bcp$posterior.prob)

EIOS_Ecu_allcp <- EIOS_Ecu %>% mutate(all_cp = ifelse(bcp.postprob > 0.5, date, NA))
ggplot(data = EIOS_Ecu_allcp, aes(x = date, y = counts)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "EIOS data for Ecuador with posterior probability of change point",
       caption = "Dashed vertical lines mark all change points") + 
  geom_vline(xintercept = na.omit(EIOS_Ecu_allcp$all_cp), lty = 2, col = "red")

# cp determined by special criteria
EIOS_Ecu$bcp_criteria <- NA

for(i in 4:nrow(EIOS_Ecu)){
  if(EIOS_Ecu$bcp.postprob[i] >= 0.5 & EIOS_Ecu$bcp.postprob[(i-1)] < 0.5 & # transition from non-epidemic to epidemic
     mean(EIOS_Ecu$counts[(i-2):(i+2)]) > mean(EIOS_Ecu$counts[(i-3):i]) & # running mean to ensure that curve is rising
     sum(!is.na(EIOS_Ecu$bcp_criteria[(i-4):(i-1)])) == 0){ # No outbreak flagged during the previous 6 weeks
    EIOS_Ecu$bcp_criteria[i] <- EIOS_Ecu$date[i]
  } else {
    EIOS_Ecu$bcp_criteria[i] <- NA
  }
}
ggplot(data = EIOS_Ecu, aes(x = date, y = counts)) + 
  geom_line(aes(col = bcp.postprob), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "4 months") + 
  labs(x = "", y = "influenza case counts", title = "EIOS data for Ecuador with posterior probability of change point") + 
  geom_vline(xintercept = na.omit(EIOS_Ecu$bcp_criteria), lty = 2, col = "red")


## cpm
cpm_USA <- processStream(EIOS_USA$counts, cpmType = "Kolmogorov-Smirnov", startup = 5, ARL0 = 500)
# Mann-Whitney Kolmogorov-Smirnov, Cramer-von-Mises give all very similar results
EIOS_USA$cpm <- NA
EIOS_USA$cpm[cpm_USA$changePoints] <- EIOS_USA$date[cpm_USA$changePoints]

ggplot(data = EIOS_USA, aes(x = date, y = counts)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "Event counts", title = "EIOS data for USA with change points, cpm package") + 
  geom_vline(xintercept = na.omit(EIOS_USA$cpm), lty = 2, col = "red")

# only CPs in rising curves
EIOS_USA$cpm_criteria <- NA
for(i in 4:nrow(EIOS_USA)){
  # use running mean of +/- 3 observations to ensure that we are in the inreasing part of the curve
  if(is.na(EIOS_USA$cpm[i]) == FALSE & 
     mean(EIOS_USA$counts[(i-3):i]) > mean(EIOS_USA$counts[(i-4):(i-1)])){
    EIOS_USA$cpm_criteria[i] <- EIOS_USA$date[i]
  } else {
    EIOS_USA$cpm_criteria[i] <- NA
  }
}
ggplot(data = EIOS_USA, aes(x = date, y = counts)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "Event counts", title = "EIOS data for USA with change points, cpm package") + 
  geom_vline(xintercept = na.omit(EIOS_USA$cpm_criteria), lty = 2, col = "red")


#Ecuador
cpm_Ecu <- processStream(EIOS_Ecu$counts, cpmType = "Mann-Whitney", startup = 5, ARL0 = 500)
# Mann-Whitney works well, Lepage, Mood, Kolmogorov-Smirnov, Cramer-von-Mises all give weird regular signals which makes no sense at all
EIOS_Ecu$cpm <- NA
EIOS_Ecu$cpm[cpm_Ecu$changePoints] <- EIOS_Ecu$date[cpm_Ecu$changePoints]

ggplot(data = EIOS_Ecu, aes(x = date, y = counts)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "Event counts", title = "EIOS data for Ecuador with change points, cpm package") + 
  geom_vline(xintercept = na.omit(EIOS_Ecu$cpm), lty = 2, col = "red")

# only CPs in rising curves
EIOS_Ecu$cpm_criteria <- NA
for(i in 6:nrow(EIOS_Ecu)){
  # use running mean of +/- 5 observations to ensure that we are in the inreasing part of the curve
  if(is.na(EIOS_Ecu$cpm[i]) == FALSE & mean(EIOS_Ecu$counts[(i-3):i]) > mean(EIOS_Ecu$counts[(i-4):(i-1)])){
    EIOS_Ecu$cpm_criteria[i] <- EIOS_Ecu$date[i]
  } else {
    EIOS_Ecu$cpm_criteria[i] <- NA
  }
}
ggplot(data = EIOS_Ecu, aes(x = date, y = counts)) + 
  geom_line(size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "Event counts", title = "EIOS data for Ecuador with change points, cpm package") + 
  geom_vline(xintercept = na.omit(EIOS_Ecu$cpm_criteria), lty = 2, col = "red")
