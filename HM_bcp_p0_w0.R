library(plyr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(bcp)
library(readxl)
library(qcc)
library(surveillance)
library(cpm)
library(tidyr)

setwd("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet")
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

overseas_territories <- c("Bermuda [UK]", "CollectivitÃ© d'outre-mer de Saint BarthÃ©lemy, France", "Cayman Islands [UK]", 
                          "Pays d'outre-mer de French Polynesia, France", "RÃ©gion d'outre-mer de Mayotte, France", 
                          "RÃ©gion d'outre-mer de French Guiana, France", "RÃ©gion d'outre-mer de RÃ©union, France", 
                          "RÃ©gion d'outre-mer de Guadeloupe, France", "RÃ©gion d'outre-mer de Martinique, France", 
                          "American Samoa [USA]", "Northern Mariana Islands [United States]", "Guam [USA]")

dataHM <- filter(dataHM, !(country %in% overseas_territories | place_name %in% overseas_territories))
dataHM$country <- factor(dataHM$country)


HM_byweek <- dataHM %>% group_by(date = floor_date(load_date, "week", week_start = 1), country = country, .drop = FALSE) %>%
  summarize(counts=n()) %>% as.data.frame()
HM_byweek <- HM_byweek[order(HM_byweek$country), ]

HM_byweek$count_smooth <- NA
for (i in seq_along(country_list)) { 
  HM_temp <- filter(HM_byweek, HM_byweek$country==country_list[i])
  count_smooth_temp <- loess(HM_temp$counts ~ as.numeric(HM_temp$date), span = 0.07, degree = 2)
  HM_byweek$count_smooth[HM_byweek$country==country_list[i]] <- count_smooth_temp$fitted
}

country_list <- levels(HM_byweek$country)


##### functions #####

epi_start <- function(df, col, cutoff = 0.5){# col in number
  for(j in 2:13){ # first few rows without looking back for previous outbreak
    if(df[j, col] >= cutoff & df[(j-1), col] < cutoff &
       df$count_smooth[j] > df$count_smooth[j-1]){
      bcp_start[j] <- df$date[j]
    } else {
      bcp_start[j] <- NA
    }
  }
  for(j in 10:nrow(df)){
    if(df[j, col] >= cutoff & df[(j-1), col] < cutoff & # transition from non-epidemic to epidemic
       df$count_smooth[j] > df$count_smooth[j-1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(bcp_start[(j-10):(j-1)])) == 0){ # No outbreak flagged during the previous 10 weeks
      bcp_start[j] <- df$date[j]
    } else {
      bcp_start[j] <- NA
    }
  }
  return(bcp_start)
  
}

epi_end <- function(df, col, cutoff = 0.5){ 
  for(j in (nrow(df)-2):2){# run in reverse because otherwise, third criterion cannot be recognized
    if(df[j-1, col] >= cutoff & df[j, col] < cutoff & # transition from non-epidemic to epidemic
       df$count_smooth[j] > df$count_smooth[j+1] & # running mean to ensure that curve is rising (beware of spikes!)
       sum(!is.na(bcp_end[(j+1):(j+15)])) == 0){ # No outbreak flagged during the previous 15 weeks
      bcp_end[j] <- df$date[j]
    } else {
      bcp_end[j] <- NA
    }
  }
  return(bcp_end)
}

##### varying bcp p0 and w0 values ##### 
### start and end of epidemics in HealthMap

p_prior <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)
w_prior <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)
p_prior_col <- paste("p_prior", p_prior, sep = "")
w_prior_col <- paste("w_prior", w_prior, sep = "")
HM_byweek <- cbind(HM_byweek, setNames(lapply(p_prior_col, function(x) x=NA), p_prior_col))
HM_byweek <- cbind(HM_byweek, setNames(lapply(w_prior_col, function(x) x=NA), w_prior_col))

bcp_list_p <- vector(mode = "list")
bcp_list_w <- vector(mode = "list")
bcp_postprob_list_p <- vector(mode = "list")
bcp_postprob_list_w <- vector(mode = "list")

for (i in seq_along(country_list)) { 
  HM_temp <- filter(HM_byweek, country == country_list[i])
  
  bcp_list_p <- lapply(p_prior, function(x) bcp(HM_temp$counts, p0 = x, burnin = 100, mcmc = 500))
  bcp_postprob_list_p <- lapply(bcp_list_p, '[[', "posterior.prob")
  bcp_postprob_df_p <- data.frame(matrix(unlist(bcp_postprob_list_p), nrow = 341, byrow = FALSE))
  HM_byweek[HM_byweek$country == country_list[i], 5:10] <- bcp_postprob_df_p
  
  bcp_list_w <- lapply(w_prior, function(x) bcp(HM_temp$counts, w0 = x, burnin = 100, mcmc = 500))
  bcp_postprob_list_w <- lapply(bcp_list_w, '[[', "posterior.prob")
  bcp_postprob_df_w <- data.frame(matrix(unlist(bcp_postprob_list_w), nrow = 341, byrow = FALSE))
  HM_byweek[HM_byweek$country == country_list[i], 11:16] <- bcp_postprob_df_w
}


start_col <- paste("bcp_start", rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)), sep = "_")
end_col <- paste("bcp_end", rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)), sep = "_")
HM_byweek <- cbind(HM_byweek, setNames(lapply(start_col, function(x) x=NA), start_col))
HM_byweek <- cbind(HM_byweek, setNames(lapply(end_col, function(x) x=NA), end_col))

bcp_start_list <- vector(mode = "list")
bcp_end_list <- vector(mode = "list")

for (i in seq_along(country_list)) {
  HM_temp <- filter(HM_byweek, country == country_list[i] & is.na(HM_byweek$p_prior0.05) == FALSE)
  
  # start of epidemics
  bcp_start_list <- lapply(5:16, function(x) epi_start(HM_temp, col = x))
  bcp_start_df <- data.frame(matrix(unlist(bcp_start_list), nrow = 341, byrow = FALSE))
  
  HM_byweek[HM_byweek$country == country_list[i], 17:28] <- bcp_start_df
  
  # end of epidemics
  bcp_end_list <- lapply(5:16, function(x) epi_end(HM_temp, col = x))
  bcp_end_df <- data.frame(matrix(unlist(bcp_end_list), nrow = 341, byrow = FALSE))
  
  HM_byweek[HM_byweek$country == country_list[i], 29:40] <- bcp_end_df
}

# replace all 0 with NA
HM_byweek[17:40] <- sapply(HM_byweek[17:40], na_if, 0)

# start end indicator column
patterns <- paste(rep(c("p", "w"), each = 6), rep(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8), 2), sep = "_")

start_end_list <- vector(mode = "list")
for(pattern in patterns){
  HM_temp <- HM_byweek[, grep(names(HM_byweek), pattern = pattern)]
  start_end_list[[pattern]] <- ifelse(is.na(HM_temp[1]) == FALSE, "start", 
                                      ifelse(is.na(HM_temp[2]) == FALSE, "end", NA))
}

start_end_col <- paste("start_end", patterns, sep = "_")
HM_byweek <- cbind(HM_byweek, setNames(lapply(start_end_col, function(x) x=NA), start_end_col))
HM_byweek[41:52] <- data.frame(matrix(unlist(start_end_list), nrow = 8184, byrow = FALSE), stringsAsFactors = FALSE)

# for spikes: if 'start' is not followed by an 'end' within 30 weeks, add an 'end' directly after - this period should be optimized
# for(i in 1:nrow(HM_byweek)){
#   if(HM_byweek$date[i] < "2019-09-01"){
#     if(is.na(HM_byweek[i, 40:56]) == FALSE){
#       if(HM_byweek[i, 40:56] == "start" & sum(HM_byweek[((i+1):(i+30)), 40:56] == "end", na.rm = TRUE) == 0){
#         HM_byweek[(i+1), 40:56] <- "end"
#       }
#     }
#   }
# }


# epidemic indicator

HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.05))) %>%
  mutate(epidemic_p_0.05 = replace(start_end_p_0.05, first(start_end_p_0.05) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.1))) %>%
  mutate(epidemic_p_0.1 = replace(start_end_p_0.1, first(start_end_p_0.1) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.2))) %>%
  mutate(epidemic_p_0.2 = replace(start_end_p_0.2, first(start_end_p_0.2) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.3))) %>%
  mutate(epidemic_p_0.3 = replace(start_end_p_0.3, first(start_end_p_0.3) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.5))) %>%
  mutate(epidemic_p_0.5 = replace(start_end_p_0.5, first(start_end_p_0.5) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_p_0.8))) %>%
  mutate(epidemic_p_0.8 = replace(start_end_p_0.8, first(start_end_p_0.8) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)

HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.05))) %>%
  mutate(epidemic_w_0.05 = replace(start_end_w_0.05, first(start_end_w_0.05) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.1))) %>%
  mutate(epidemic_w_0.1 = replace(start_end_w_0.1, first(start_end_w_0.1) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.2))) %>%
  mutate(epidemic_w_0.2 = replace(start_end_w_0.2, first(start_end_w_0.2) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.3))) %>%
  mutate(epidemic_w_0.3 = replace(start_end_w_0.3, first(start_end_w_0.3) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.5))) %>%
  mutate(epidemic_w_0.5 = replace(start_end_w_0.5, first(start_end_w_0.5) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)
HM_byweek <- HM_byweek %>% 
  group_by(country, grp = cumsum(!is.na(start_end_w_0.8))) %>%
  mutate(epidemic_w_0.8 = replace(start_end_w_0.8, first(start_end_w_0.8) == 'start', TRUE)) %>% 
  ungroup() %>% 
  select(-grp)

HM_byweek<- HM_byweek %>% mutate_at(vars(53:64), ~replace(., is.na(.), FALSE))
HM_byweek<- HM_byweek %>% mutate_at(vars(53:64), ~replace(., . == "end", TRUE))


write.csv(HM_byweek, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HM_epidemic_bcp_p0_w0.csv")

