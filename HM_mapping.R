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
locations_sf <- st_as_sf(dataHM, coords = c("lon", "lat"), crs = 4326)
mapview(locations_sf, cex = 3)



countries_sp <- ne_countries(scale = "medium")


coords <- select(dataHM, lon, lat)

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
plot(countries_sf)
plot(locations_sf)

ggplot() +
  geom_sf(data = countries_sf) +
  geom_sf(data = locations_sf, alpha = 0.5, col = "darkgreen") +
  labs(title = "HealthMap events from 2013 - 2019")
