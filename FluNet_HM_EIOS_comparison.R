library(ggplot2)
library(plyr)
library(readxl)
library(dplyr)
library(lubridate)
library(stringr)


## load FluNet data
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


HM_byweek <- dataHM %>% group_by(date = floor_date(load_date, "week"), country = country, .drop = FALSE) %>%
  summarize(counts=n()) %>% as.data.frame()

HM_byweek <- dplyr::filter(HM_byweek, !(country %in% overseas_territories))
HM_byweek$country <- factor(HM_byweek$country)


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



## plot all three in one plot
for (i in seq_along(country_list)) { 
  FluNet_data_temp <- subset(FluNet_data, FluNet_data$Country==country_list[i])
  EIOS_byweek_temp <- subset(EIOS_byweek, EIOS_byweek$country==country_list[i])
  HM_byweek_temp <- subset(HM_byweek, HM_byweek$country==country_list[i])
  
  
  counts_temp <- c(HM_byweek_temp$counts, EIOS_byweek_temp$counts, FluNet_data_temp$ALL_INF)
  source_temp <- c(rep("HealthMap", nrow(HM_byweek_temp)), rep("EIOS", nrow(EIOS_byweek_temp)), 
                   rep("FluNet", nrow(FluNet_data_temp)))
  country_temp <- rep(country_list[i], nrow(HM_byweek_temp) + nrow(EIOS_byweek_temp) + nrow(FluNet_data_temp))
  date_temp <- c(HM_byweek_temp$date, EIOS_byweek_temp$date, FluNet_data_temp$SDATE)
  
  temp_df <- data.frame(date_temp, country_temp, counts_temp, source_temp)
  
  plot <- ggplot(temp_df, aes(x = date_temp, y = counts_temp)) + 
    geom_area(fill="#69b3a2", alpha=0.5) +
    geom_line(color="#69b3a2", size = 1) + 
    facet_wrap(facets = ~source_temp, nrow = 3, scales = "free_y", strip.position = "left", 
               labeller = as_labeller(c(HealthMap = "HealthMap events", EIOS = "EIOS events", FluNet = "FluNet counts"))) +
    ylab(NULL) +
    theme(strip.background = element_blank(), strip.placement = "outside") + 
    labs(title = paste("Comparison of HealthMap, EIOS and FluNet counts for", country_list[i], sep = " ")) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "12 months", date_minor_breaks = "2 months") + 
    xlab("")
  
  print(plot)
  ggsave(plot = plot, file = paste("HM_EIOS_FluNet_comparison", country_list[i], ".jpeg", sep=' '))
}
