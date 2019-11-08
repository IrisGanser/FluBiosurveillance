library(ggplot2)
library(plyr)
library(readxl)
library(dplyr)
library(lubridate)


## USA
FluNetUSA <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet/FluNet_Usa.csv")
str(FluNetUSA)
FluNetUSA$SDATE <- as.Date(FluNetUSA$SDATE)
FluNetUSA$EDATE <- as.Date(FluNetUSA$EDATE)

seasons <- seq(2014, 2019)
season.n <- NULL
for (season in seasons) {
  n.tests <- sum(FluNetUSA$ALL_INF[FluNetUSA$Year == season])
  season.n <- c(season.n, FluNetUSA$ALL_INF[FluNetUSA$Year == season] > n.tests*0.01)
}
FluNetUSA$neuzil <- season.n

n.2018 <- sum(FluNetUSA$ALL_INF[FluNetUSA$Year==2018])
season.2018 <- FluNetUSA$ALL_INF[FluNetUSA$Year==2018] > n.2018*0.01

ggplot(FluNetUSA, aes(x = SDATE)) + 
  geom_line(aes(y = ALL_INF, col = "Total influenza")) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "FluNet data USA 2014-2019") + 
  geom_line(aes(y = INF_A, col = "Influenza A")) + 
  geom_line(aes(y = INF_B, col = "Influenza B")) +
  scale_color_discrete(name = "")


## Australia
FluNetAus <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet/FluNet_Aus.csv")
FluNetAus$SDATE <- as.Date(FluNetAus$SDATE)
FluNetAus$EDATE <- as.Date(FluNetAus$EDATE)

ggplot(FluNetAus, aes(x = SDATE)) + 
  geom_line(aes(y = ALL_INF, col = "Total influenza")) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "FluNet data Australia 2014-2019") + 
  geom_line(aes(y = INF_A, col = "Influenza A")) + 
  geom_line(aes(y = INF_B, col = "Influenza B")) +
  scale_color_discrete(name = "")


## Bulgaria
FluNetBul <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet/FluNet_Bul.csv")
FluNetBul$SDATE <- as.Date(FluNetBul$SDATE)
FluNetBul$EDATE <- as.Date(FluNetBul$EDATE)

ggplot(FluNetBul, aes(x = SDATE)) + 
  geom_line(aes(y = ALL_INF, col = "Total influenza")) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "FluNet data Bulgaria") + 
  geom_line(aes(y = INF_A, col = "Influenza A")) + 
  geom_line(aes(y = INF_B, col = "Influenza B")) +
  scale_color_discrete(name = "")


## Nigeria
FluNetNgr <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet/FluNet_Nig.csv")
FluNetNgr$SDATE <- as.Date(FluNetNgr$SDATE)
FluNetNgr$EDATE <- as.Date(FluNetNgr$EDATE)

ggplot(FluNetNgr, aes(x = SDATE)) + 
  geom_line(aes(y = ALL_INF, col = "Total influenza")) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 year") + 
  labs(x = "", y = "influenza case counts", title = "FluNet data Nigeria") + 
  geom_line(aes(y = INF_A, col = "Influenza A")) + 
  geom_line(aes(y = INF_B, col = "Influenza B")) +
  scale_color_discrete(name = "")


# all four pilot countries
FluNetdata <- rbind(FluNetAus, FluNetBul, FluNetNgr, FluNetUSA)
ggplot(FluNetdata, aes(x = SDATE)) + 
  geom_line(aes(y = ALL_INF, col = "Total influenza")) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 years") + 
  labs(x = "", y = "influenza case counts", title = "WHO FluNet Data") + 
  geom_line(aes(y = INF_A, col = "Influenza A")) + 
  geom_line(aes(y = INF_B, col = "Influenza B")) +
  scale_color_discrete(name = "") + 
  facet_wrap(facets = ~Country, scales = "free_y")