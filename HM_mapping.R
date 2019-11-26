library(sf)
library(sp)
library(mapview)
library(ggmap)
library(ggplot2)
library(plyr)
library(readxl)
library(dplyr)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(RColorBrewer)


setwd("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet")

# load HealthMap data
dataHM <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/mcgill flu iris 20190711144158.csv", header = F)
colnames(dataHM) <- c("place_name", "country", "disease_name", "species_name", "alert_id", "summary", "href", "issue_date", "load_date",   
                      "smooshed", "descr", "long_name", "lon", "lat")

dataHM$country <- as.factor(dataHM$country)
dataHM$long_name <- as.factor(dataHM$long_name)
dataHM$load_date <- as.POSIXct(dataHM$load_date)
dataHM$issue_date <- as.POSIXct(dataHM$issue_date)


# de-duplication based on alert IDs
dataHM <- dataHM[!duplicated(dataHM$alert_id), ]

# filter out overseas territories
overseas_territories <- c("Bermuda [UK]", "CollectivitÃ© d'outre-mer de Saint BarthÃ©lemy, France", "Cayman Islands [UK]", 
                          "Pays d'outre-mer de French Polynesia, France", "RÃ©gion d'outre-mer de Mayotte, France", 
                          "RÃ©gion d'outre-mer de French Guiana, France", "RÃ©gion d'outre-mer de RÃ©union, France", 
                          "RÃ©gion d'outre-mer de Guadeloupe, France", "RÃ©gion d'outre-mer de Martinique, France")
dataHM <- filter(dataHM, !(country %in% overseas_territories | place_name %in% overseas_territories))
dataHM$country <- factor(dataHM$country)

# make locations ready fror plotting
HM_events_sf <- st_as_sf(dataHM, coords = c("lon", "lat"), crs = 4326)
mapview <- mapview(HM_events_sf, cex = 3)
mapview
# mapshot(mapview, url = "HMevents_map.html", selfcontained = FALSE)
# mapshot(mapview, file = "HMevents_map.png", remove_controls = c("zoomControl", "layersControl", "homeButton"), selfcontained = FALSE)



## make a map with sp package - obsolete 
#countries_sp <- ne_countries(scale = "medium")
#coords <- select(dataHM, lon, lat)

# Create SpatialPoints object with coords and CRS
#points_sp <- SpatialPoints(coords = coords,
 #                          proj4string = CRS("+proj=longlat +datum=WGS84"))
#points_spdf <- SpatialPointsDataFrame(coords = coords,
  #                                    data = dataHM,  
   #                                   proj4string = CRS("+proj=longlat +datum=WGS84"))


#plot(countries_sp, col = gray(0.8), border = gray(0.7))
#plot(points_spdf, pch = 20, add = TRUE)
#box()


countries_sf <- ne_countries(scale = "medium", returnclass = "sf")
plot(countries_sf["iso_a3"], axes = TRUE, main = "")


#ggplot() +
 # geom_sf(data = countries_sf, aes(fill = region_un), alpha = 0.8) +
  #geom_sf(data = HM_events_sf, alpha = 0.5, col = "black") +
  #labs(title = "HealthMap events from 2013 - 2019", fill = "UN region") + 
  #theme(legend.position = "bottom")


influenza_zones <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/InfluenzaTransmissionZones.csv", header = TRUE, sep = ";", 
                            stringsAsFactors = FALSE) %>% subset(select = 1:3)
influenza_zones$influenza_transmission_zone <- as.factor(influenza_zones$influenza_transmission_zone)

countries_sf <- left_join(countries_sf, influenza_zones, by = c("iso_a3" = "ISO"))

countries_sf$geounit[which(is.na(countries_sf$influenza_transmission_zone))]
countries_sf$influenza_transmission_zone[which(is.na(countries_sf$influenza_transmission_zone))] <- c("Oceania Melanesia and Polynesia", 
                                                                                                      "Western Asia", "Southern Asia", "Southern Asia", 
                                                                                                      "Eastern Europe", "Eastern Africa")
countries_sf$influenza_transmission_zone[which(countries_sf$influenza_transmission_zone == "Antarctica (none)")] <- NA


colourCount = length(levels(countries_sf$influenza_transmission_zone))
getPalette = colorRampPalette(brewer.pal(19, "Set1"))

ggplot() +
  geom_sf(data = countries_sf, aes(fill = influenza_transmission_zone), alpha = 0.8) +
  scale_fill_manual(values = getPalette(colourCount)) +
  geom_sf(data = HM_events_sf, alpha = 0.5, col = "black") +
  labs(title = "Spatial distribution of HealthMap events from 2013 - 2019", fill = "Influenza transmission zone") + 
  theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

#ggsave(filename = "mapHM_CBSet1.png", scale = 1.5)
  

ggplot() +
  geom_sf(data = countries_sf, aes(fill = influenza_transmission_zone), alpha = 0.8) +
  #scale_fill_manual(values = getPalette(colourCount)) +
  geom_sf(data = HM_events_sf, alpha = 0.5, col = "black") +
  labs(title = "Spatial distribution of HealthMap events from 2013 - 2019", fill = "Influenza transmission zone") + 
  theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

#ggsave(filename = "mapHM_defaultcol.png", scale = 1.5)
