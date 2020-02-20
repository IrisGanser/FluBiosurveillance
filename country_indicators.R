library(ggplot2)
library(plyr)
library(readxl)
library(lubridate)
library(stringr)
library(dplyr)
library(tidyr)

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

HM_byweek <- dataHM %>% group_by(date = floor_date(load_date, "week"), country = country, .drop = FALSE) %>%
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
EIOS_byweek <- EIOSreports %>% group_by(date = floor_date(importDate, "week"), country = Country, .drop = FALSE) %>%
  summarize(counts=n()) %>% as.data.frame()


## visualize total counts of every source
FluNet_total <- FluNet_data %>% group_by(Country) %>% summarize(total = sum(ALL_INF), max = max(ALL_INF))
ggplot(FluNet_total, aes(x = Country, y = total)) + 
  geom_col(fill = "darkorange3") + 
  geom_text(aes(label=total), position = position_stack(0.5)) + 
  labs(title = "Total number of FluNet Influenza counts from January 2013 - December 2019", x = "") + 
  coord_flip() +
  scale_x_discrete(limits = rev(levels(FluNet_data$Country)))

ggplot(FluNet_total, aes(x = Country, y = max)) + 
  geom_col(fill = "darkorange3") + 
  geom_text(aes(label=max), position = position_stack(0.5)) + 
  labs(title = "Maximum weekly FluNet Influenza count from January 2013 - December 2019", x = "") + 
  coord_flip() +
  scale_x_discrete(limits = rev(levels(FluNet_data$Country)))

#FluNet_total_long <- FluNet_total %>% pivot_longer(cols = c(total, max), names_to = "count_category", values_to = "counts")

HM_total <- HM_byweek %>% group_by(country) %>% summarize(total = sum(counts), max = max(counts))
ggplot(HM_total, aes(x = country, y = total)) + 
  geom_col(fill = "darkorange3") + 
  geom_text(aes(label=total), position = position_stack(0.5)) + 
  labs(title = "Total number of HealthMap events from January 2013 - July 2019", x = "") + 
  coord_flip() +
  scale_x_discrete(limits = rev(levels(dataHM$country))) + 
  scale_y_continuous(limits = c(0, 9300))

ggplot(HM_total, aes(x = country, y = max)) + 
  geom_col(fill = "darkorange3") + 
  geom_text(aes(label=max), position = position_stack(0.5)) + 
  labs(title = "Maximum weekly number of HealthMap events from January 2013 - July 2019", x = "") + 
  coord_flip() +
  scale_x_discrete(limits = rev(levels(dataHM$country)))

EIOS_total <- EIOS_byweek %>% group_by(country) %>% summarize(total = sum(counts), max = max(counts))
ggplot(EIOS_total, aes(x = country, y = total)) + 
  geom_col(fill = "darkorange3") + 
  geom_text(aes(label=total), position = position_stack(0.5)) + 
  labs(title = "Total number of EIOS events from November 2017 - December 2019", x = "") + 
  coord_flip() +
  scale_x_discrete(limits = rev(levels(dataHM$country))) 

ggplot(EIOS_total, aes(x = country, y = max)) + 
  geom_col(fill = "darkorange3") + 
  geom_text(aes(label=max), position = position_stack(0.5)) + 
  labs(title = "Maximum weekly number of EIOS events from November 2017 - December 2019", x = "") + 
  coord_flip() +
  scale_x_discrete(limits = rev(levels(dataHM$country))) 


# divide into categories of data abundance ('high', 'medium', 'low')
FluNet_total$total_cat <- cut(FluNet_total$total, 
                              breaks = c(0, quantile(FluNet_total$total, probs = 0.5), quantile(FluNet_total$total, probs = 0.75), Inf),
                              labels = c("low", "medium", "high"))

HM_total$total_cat <- cut(HM_total$total, 
                          breaks = c(0, quantile(HM_total$total, probs = 0.5), quantile(HM_total$total, probs = 0.75), Inf),
                          labels = c("low", "medium", "high"))

EIOS_total$total_cat <- cut(EIOS_total$total, 
                            breaks = c(0, quantile(EIOS_total$total, probs = 0.5), quantile(EIOS_total$total, probs = 0.75), Inf),
                            labels = c("low", "medium", "high"))


## influenza transmission zone
influenza_zones <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/InfluenzaTransmissionZones.csv", header = TRUE, sep = ";", 
                            stringsAsFactors = FALSE) %>% subset(select = 1:3)
influenza_zones$influenza_transmission_zone <- as.factor(influenza_zones$influenza_transmission_zone)

influenza_zones$country_name <- revalue(influenza_zones$country_name, 
                                        c("Iran (Islamic Republic of)" = "Iran", "Russian Federation (the)" = "Russia", 
                                          "United Kingdom of Great Britain and Northern Ireland (the)" = "United Kingdom",
                                          "United States of America (the)" = "United States", "Viet Nam" = "Vietnam"))
influenza_zones$country_name <- as.factor(influenza_zones$country_name)


indicators <- data.frame(country = FluNet_total$Country, FluNet_total_cat = FluNet_total$total_cat, 
                         HM_total_cat = HM_total$total_cat, EIOS_total_cat = EIOS_total$total_cat)
indicators <- left_join(indicators, influenza_zones, by = c("country" = "country_name"))
indicators$country <- as.factor(indicators$country)
indicators$influenza_transmission_zone <- as.factor(indicators$influenza_transmission_zone)

indicators$global_region <- NA
for(i in 1:nrow(indicators)){
  if(indicators$influenza_transmission_zone[i] %in% c("Central America Caribbean", "Central Asia", "Eastern Asia", "Eastern Europe",
                                                   "North America", "Northern Africa", "Northern Europe", "South West Europe", 
                                                   "Western Asia")){
    indicators$global_region[i] <- "temperate Northern hemisphere"
  } else if(indicators$influenza_transmission_zone[i] %in% c("Eastern Africa", "Middle Africa", "South East Asia", "Southern Asia",
                                                             "Tropical South America", "Western Africa")){
    indicators$global_region[i] <- "tropical"
  } else{
    indicators$global_region[i] <- "temperate Southern hemisphere"
  }
}

indicators$english <- NA
for(i in 1:nrow(indicators)){
  if(indicators$country[i] %in% c("Australia", "United Kingdom", "United States", "Nigeria")){
    indicators$english[i] <- TRUE
  }else{
    indicators$english[i] <- FALSE
  }
}

languages <- c("Spanish", "English", "Portuguese", "Bulgarian", "Mandarin", "Spanish", "Spanish", "Arabic", "French",
               "German", "Greek", "Hindi", "Farsi", "Spanish", "English", "Russian", "Arabic", "Afrikaans", "Swedish",
               "Thai", "English", "English", "Spanish", "Vietnamese")
indicators$language <- languages


indicators$problematic_FluNet <- ifelse(indicators$FluNet_total_cat == "low" | indicators$global_region == "tropical", 
                                        TRUE, FALSE)
indicators$problematic_FluNet[indicators$country == "France"] <- TRUE
not_problematic <- c("Germany", "Ecuador", "Bulgaria", "Greece", "Iran", "South Africa")
indicators$problematic_FluNet[indicators$country %in% not_problematic] <- FALSE
indicators$country[indicators$problematic_FluNet == TRUE]

indicators$problematic_EBS <- ifelse(indicators$HM_total_cat == "low" | indicators$EIOS_total_cat == "low", TRUE, FALSE)
indicators$country[indicators$problematic_EBS == TRUE]
