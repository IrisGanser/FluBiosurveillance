library(psych)
library(ggplot2)
library(plyr)
library(readxl)
library(lubridate)
library(stringr)
library(dplyr)
library(tidyr)
library(irr)

## load all data
# load FluNet data
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


# summarize into weekly counts and filter out overseas territories
overseas_territories <- c("Bermuda [UK]", "CollectivitÃ© d'outre-mer de Saint BarthÃ©lemy, France", "Cayman Islands [UK]", 
                          "Pays d'outre-mer de French Polynesia, France", "RÃ©gion d'outre-mer de Mayotte, France", 
                          "RÃ©gion d'outre-mer de French Guiana, France", "RÃ©gion d'outre-mer de RÃ©union, France", 
                          "RÃ©gion d'outre-mer de Guadeloupe, France", "RÃ©gion d'outre-mer de Martinique, France", 
                          "American Samoa [USA]", "Northern Mariana Islands [United States]", "Guam [USA]")
dataHM <- filter(dataHM, !(country %in% overseas_territories | place_name %in% overseas_territories))
dataHM$country <- factor(dataHM$country)

HM_byweek <- dataHM %>% group_by(date = floor_date(load_date, "week", week_start = 1), country = country, .drop = FALSE) %>%
  summarize(counts=n()) %>% as.data.frame()



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

# summarize into weekly data
EIOS_byweek <- EIOSreports %>% group_by(date = floor_date(importDate, "week", week_start = 1), country = Country, .drop = FALSE) %>%
  summarize(counts=n()) %>% as.data.frame()


#### ICC ####
HM_USA <- filter(HM_byweek, country == "United States")
FluNet_USA <- filter(FluNet_data, Country == "United States")
EIOS_USA <- filter(EIOS_byweek, country == "United States")

icc_USA_HM <- left_join(HM_USA, FluNet_USA, by = c("date" = "SDATE")) %>% select(counts, ALL_INF)
icc_USA_EIOS <- left_join(EIOS_USA, FluNet_USA, by = c("date" = "SDATE")) %>% select(counts, ALL_INF)

icc(ratings = icc_USA_HM, model = "twoway", type = "consistency", unit = "single")
icc(ratings = icc_USA_EIOS, model = "twoway", type = "consistency", unit = "single")

b <- ICC(icc_USA_EIOS)
b$results$ICC[3]


ICC_HM <- numeric(length = 24)
ICC_EIOS <- numeric(length = 24)
for(i in seq_along(country_list)){
  HM_temp <- filter(HM_byweek, country == country_list[i])
  FluNet_temp <- filter(FluNet_data, Country == country_list[i])
  EIOS_temp <- filter(EIOS_byweek, country == country_list[i])
  
  icc_df_HM <- left_join(HM_temp, FluNet_temp, by = c("date" = "SDATE")) %>% select(counts, ALL_INF)
  icc_df_EIOS <- left_join(EIOS_temp, FluNet_temp, by = c("date" = "SDATE")) %>% select(counts, ALL_INF)
  
  icc_HM <- ICC(icc_df_HM)
  icc_EIOS <- ICC(icc_df_EIOS)
  
  ICC_HM[i] <- icc_HM$results$ICC[3]
  ICC_EIOS[i] <- icc_EIOS$results$ICC[3]
}

ICC_plot <- data.frame(country = rep(country_list, 2), ICC = c(ICC_HM, ICC_EIOS), 
                       source = c(rep("HealthMap", 24), rep("EIOS", 24)))
ICC_plot$source <- factor(ICC_plot$source, levels = c("HealthMap", "EIOS"))
ggplot(ICC_plot, aes(x = country, y = ICC, fill = source)) + 
  geom_col(position = "dodge") + 
  coord_flip() + 
  scale_fill_brewer(palette = "Set1") + 
  labs(x = "", title = "Intraclass Correlation Coefficient between FluNet and EBS systems")
ggsave("ICC.jpeg")
