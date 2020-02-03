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


# barplot of healthmap event counts in total over 5 years
ggplot(EIOSreports, aes(x = Country)) + 
  geom_bar(fill = "darkorange3") + 
  geom_text(aes(label=..count..), stat="count", position = position_stack(0)) + 
  labs(title = "Total number of EIOS events from November 2017 - November 2019", x = "") + 
  coord_flip() +
  scale_x_discrete(limits = rev(levels(EIOSreports$Country))) +
  scale_y_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000), limits = c(0, 28000))

#ggsave("EIOS total counts.jpeg")



# summarize into weekly data
EIOS_byweek <- EIOSreports %>% group_by(date = floor_date(importDate, "week"), country = Country, .drop = FALSE) %>%
  summarize(counts=n()) %>% as.data.frame()


## plot EIOS event data to compare with FluNet data
for (i in seq_along(country_list)) { 
  FluNet_data_temp <- subset(FluNet_data, FluNet_data$Country==country_list[i] & FluNet_data$SDATE > "2017-11-05")
  EIOS_byweek_temp <- subset(EIOS_byweek, EIOS_byweek$country==country_list[i])
  
  counts_temp <- c(EIOS_byweek_temp$counts, FluNet_data_temp$ALL_INF)
  source_temp <- c(rep("EIOS", nrow(EIOS_byweek_temp)), rep("FluNet", nrow(FluNet_data_temp)))
  country_temp <- rep(country_list[i], nrow(EIOS_byweek_temp) + nrow(FluNet_data_temp))
  date_temp <- c(EIOS_byweek_temp$date, FluNet_data_temp$SDATE)
  
  temp_df <- data.frame(date_temp, country_temp, counts_temp, source_temp)
  
  plot <- ggplot(temp_df, aes(x = date_temp, y = counts_temp)) + 
    geom_line() + 
    facet_wrap(facets = ~source_temp, nrow = 2, scales = "free_y", strip.position = "left", 
               labeller = as_labeller(c(EIOS = "EIOS events", FluNet = "FluNet counts"))) +
    ylab(NULL) +
    theme(strip.background = element_blank(), strip.placement = "outside") + 
    labs(title = paste("Comparison of EIOS and FluNet counts for", country_list[i], sep = " ")) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "6 months") + 
    xlab("")
  
  print(plot)
  #ggsave(plot = plot, file = paste("EIOS FluNet comparison", country_list[i], ".jpeg", sep=' '))
}


