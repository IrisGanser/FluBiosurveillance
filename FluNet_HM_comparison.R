library(ggplot2)
library(plyr)
library(readxl)
library(dplyr)
library(lubridate)
library(plotly)
library(htmlwidgets)
# library(hrbrthemes)

# 24 countries
## load data
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

## make individual plots 
country_list <- levels(FluNet_data$Country)

for (i in seq_along(country_list)) { 
  plot <- ggplot(subset(FluNet_data, FluNet_data$Country == country_list[i]), aes(x = SDATE)) + 
    geom_line(aes(y = ALL_INF, col = "Total influenza")) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "influenza case counts", title = paste("WHO FluNet data for", country_list[i], sep=' ')) + 
    geom_line(aes(y = INF_A, col = "Influenza A")) + 
    geom_line(aes(y = INF_B, col = "Influenza B")) +
    scale_color_discrete(name = "")
  print(plot)
  #ggsave(plot = plot, file = paste("WHO FluNet plot", country_list[i], ".jpeg", sep=' '))
}

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

# barplot of healthmap event counts in total over 5 years
dataHM <- filter(dataHM, !(country %in% overseas_territories | place_name %in% overseas_territories))
dataHM$country <- factor(dataHM$country)

ggplot(dataHM, aes(x = country)) + 
  geom_bar(fill = "darkorange3") + 
  geom_text(aes(label=..count..), stat="count", position = position_stack(0.5)) + 
  labs(title = "Total number of HealthMap events from January 2013 - July 2019", x = "") + 
  coord_flip() +
  scale_x_discrete(limits = rev(levels(dataHM$country))) + 
  scale_y_continuous(limits = c(0, 9300))

#ggsave("HM total counts.jpeg")


# summarize into weekly data
HM_byweek <- dataHM %>% group_by(date = floor_date(load_date, "week"), country = country, .drop = FALSE) %>%
  summarize(counts=n()) %>% as.data.frame()

HM_byweek <- dplyr::filter(HM_byweek, !(country %in% overseas_territories))
HM_byweek$country <- factor(HM_byweek$country)

## combine FluNet with HealthMap plots
# for Argentina
FluNet_data_Arg <- filter(FluNet_data, Country == "Argentina")
HM_byweek_Arg <- filter(HM_byweek, country == "Argentina")
countsArg <- c(HM_byweek_Arg$counts, FluNet_data_Arg$ALL_INF)
sourceArg <- c(rep("HealthMap", nrow(HM_byweek_Arg)), rep("WHO", nrow(FluNet_data_Arg)))
countryArg <- rep("Argentina", nrow(HM_byweek_Arg) + nrow(FluNet_data_Arg))
dateArg <- c(HM_byweek_Arg$date, FluNet_data_Arg$SDATE)

temp_data_Arg <- data.frame(dateArg, countryArg, countsArg, sourceArg) %>% filter(dateArg < "2019-06-30")

p <- ggplot(temp_data_Arg, aes(x = dateArg, y = countsArg)) + 
  geom_area(fill="#69b3a2", alpha=0.5) +
  geom_line(color="#69b3a2", size = 1) + 
  facet_wrap(facets = ~sourceArg, nrow = 2, scales = "free_y", strip.position = "left", 
             labeller = as_labeller(c(HealthMap = "HealthMap events", WHO = "WHO counts"))) +
  ylab(NULL) +
  theme(strip.background = element_blank(), strip.placement = "outside") + 
  labs(title = "Argentina comparison of HealthMap and WHO counts", caption = "Cor = 0.573") + 
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "12 months") + 
  xlab("") 
p 
p_ggplotly <- ggplotly(p)
p_ggplotly


cor(temp_data_Arg$countsArg[1:339], temp_data_Arg$countsArg[340:678], method = "spearman")
plot(temp_data_Arg$countsArg[1:339], temp_data_Arg$countsArg[340:678])


# loop for correlation calculation between HM and FluNet data
country_list <- levels(HM_byweek$country)
cor.value <- rep(NA, length(country_list))

for (i in seq_along(country_list)) { 
  FluNet_data_temp <- subset(FluNet_data, FluNet_data$Country==country_list[i]) %>% filter(SDATE < "2019-06-30")
  HM_byweek_temp <- subset(HM_byweek, HM_byweek$country==country_list[i]) %>% filter(date < "2019-06-30")
  
  counts_HM <- HM_byweek_temp$counts
  counts_WHO <- FluNet_data_temp$ALL_INF
  
  if(length(counts_HM) >= length(counts_WHO)){
    counts_HM <- counts_HM[1: length(counts_WHO)]
  } else {
    counts_WHO <- counts_WHO[1:length(counts_HM)]
  } # trim the longer of the two vectors to the size of the shorter vector
  
  cor.value[i] <- round(cor(counts_HM, counts_WHO, method = "spearman"), 3)
  
}

names(cor.value) <- country_list
cor.value


# build a loop for correlation plots for all 24 countries


for (i in seq_along(country_list)) { 
  FluNet_data_temp <- subset(FluNet_data, FluNet_data$Country==country_list[i])
  HM_byweek_temp <- subset(HM_byweek, HM_byweek$country==country_list[i])
  
  counts_temp <- c(HM_byweek_temp$counts, FluNet_data_temp$ALL_INF)
  source_temp <- c(rep("HealthMap", nrow(HM_byweek_temp)), rep("WHO", nrow(FluNet_data_temp)))
  country_temp <- rep(country_list[i], nrow(HM_byweek_temp) + nrow(FluNet_data_temp))
  date_temp <- c(HM_byweek_temp$date, FluNet_data_temp$SDATE)
  
  temp_df <- data.frame(date_temp, country_temp, counts_temp, source_temp)
  
  plot <- ggplot(temp_df, aes(x = date_temp, y = counts_temp)) + 
    geom_line() + 
    facet_wrap(facets = ~source_temp, nrow = 2, scales = "free_y", strip.position = "left", 
               labeller = as_labeller(c(HealthMap = "HealthMap events", WHO = "WHO counts"))) +
    ylab(NULL) +
    theme(strip.background = element_blank(), strip.placement = "outside") + 
    labs(title = paste("Comparison of HealthMap and WHO counts for", country_list[i], sep = " ")) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "12 months") + 
    xlab("")
    #labs(title = paste("Comparison of HealthMap and WHO counts for", country_list[i], sep = " "), 
         #caption = paste("Cor =", cor.value[i], sep = " ")) + 
    
  
  print(plot)
  
  #ggsave(plot = plot, file = paste("WHO HM comparison", country_list[i], ".jpeg", sep=' '))
}




