library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(binom)

# load epidemic datasets with outbreak indicators
FluNet_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/FluNet_epidemic.csv")
FluNet_epidemic$SDATE <- as.POSIXct(FluNet_epidemic$SDATE)
HM_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HealthMap_epidemic_all_bcp.csv")
HM_epidemic$date <- as.POSIXct(HM_epidemic$date)
EIOS_epidemic <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_epidemic_all_bcp.csv")
EIOS_epidemic$date <- as.POSIXct(EIOS_epidemic$date)

country_list <- levels(FluNet_epidemic$Country)

# revise startend columns of datasets
HM_epidemic$startend <- NA
for(i in 2:nrow(HM_epidemic)){
  if(HM_epidemic$epidemic[i] == TRUE & HM_epidemic$epidemic[i-1] == FALSE){
    HM_epidemic$startend[i] <- "start"
  }else if(HM_epidemic$epidemic[i] == TRUE & HM_epidemic$epidemic[i+1] == FALSE){
    HM_epidemic$startend[i] <- "end"
  }
}

EIOS_epidemic$startend <- NA
for(i in 2:nrow(EIOS_epidemic)){
  if(EIOS_epidemic$epidemic[i] == TRUE & EIOS_epidemic$epidemic[i-1] == FALSE){
    EIOS_epidemic$startend[i] <- "start"
  }else if(EIOS_epidemic$epidemic[i] == TRUE & EIOS_epidemic$epidemic[i+1] == FALSE){
    EIOS_epidemic$startend[i] <- "end"
  }
}

FluNet_epidemic$startend <- NA
for(i in 2:nrow(FluNet_epidemic)){
  if(FluNet_epidemic$epidemic[i] == TRUE & FluNet_epidemic$epidemic[i-1] == FALSE){
    FluNet_epidemic$startend[i] <- "start"
  }else if(FluNet_epidemic$epidemic[i] == TRUE & FluNet_epidemic$epidemic[i+1] == FALSE){
    FluNet_epidemic$startend[i] <- "end"
  }
}

#some countries end with an ongoing epidemic, so automatically, an "end" is assigned to the start of the next country, which messes up the whole data
# therefore assign an NA to the first entry of every country
# for all countries except Saudi Arabia
FluNet_epidemic$startend[which(FluNet_epidemic$SDATE == min(FluNet_epidemic$SDATE))] <- NA
# for Saudi Arabia (data start later)
sa <- which(FluNet_epidemic$Country == "Saudi Arabia")
FluNet_epidemic$startend[sa[1]] <- NA


# new epidemic plots for all countries
for (i in seq_along(country_list)) {
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  
  plot <- ggplot(HM_temp, aes(x = date, y = counts, col = epidemic)) + 
    geom_line(aes(group = 1), size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "HealthMap event counts", 
         title = paste("HealthMap data for", country_list[i], "with start and end of epidemics", sep = " ")) + 
    geom_vline(xintercept = na.omit(HM_temp$bcp_start[HM_temp$startend == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(HM_temp$bcp_end[HM_temp$startend == "end"]), lty = 2, col = "darkgreen")  + 
    scale_color_manual(values = c("#6e6868", "#e64040"))
  
  print(plot)
  #ggsave(plot = plot, file = paste("HealthMap outbreak", country_list[i], ".jpeg", sep=' '))
}


for (i in seq_along(country_list)) {
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  plot <- ggplot(EIOS_temp, aes(x = date, y = counts, col = epidemic)) + 
    geom_line(aes(group = 1), size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") + 
    labs(x = "", y = "EIOS event counts", 
         title = paste("EIOS data for", country_list[i], "with start and end of epidemics", sep = " ")) + 
    geom_vline(xintercept = na.omit(EIOS_temp$bcp_start[EIOS_temp$startend == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(EIOS_temp$bcp_end[EIOS_temp$startend == "end"]), lty = 2, col = "darkgreen") + 
    scale_color_manual(values = c("#6e6868", "#e64040"))
  
  print(plot)
  #ggsave(plot = plot, file = paste("EIOS outbreak", country_list[i], ".jpeg", sep=' '))
}


for (i in seq_along(country_list)) {
  FluNet_temp <- filter(FluNet_epidemic, FluNet_epidemic$Country==country_list[i])
  
  plot <- ggplot(FluNet_temp, aes(x = SDATE, y = ALL_INF, col = epidemic)) +
    geom_line(aes(group = 1), size = 0.75) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 year") +
    labs(x = "", y = "influenza case counts",
         title = paste("WHO FluNet data for", country_list[i], "with epidemics", sep = " ")) +
    geom_vline(xintercept = na.omit(FluNet_temp$bcp_start[FluNet_temp$startend == "start"]), lty = 2, col = "red") +
    geom_vline(xintercept = na.omit(FluNet_temp$bcp_end[FluNet_temp$startend == "end"]), lty = 2, col = "darkgreen") + 
    scale_color_manual(values = c("#6e6868", "#e64040"))
  
  print(plot)
  #ggsave(plot = plot, file = paste("FluNet outbreak", country_list[i], ".jpeg", sep=' '))
}


FluNet_plot <- filter(FluNet_epidemic, Country %in% c("Nigeria", "United States", "Thailand", "Vietnam"))
FluNet_plot$Country <- factor(FluNet_plot$Country, levels = c("United States", "Nigeria", "Thailand", "Vietnam"))

ggplot(FluNet_plot, aes(x = SDATE, y = ALL_INF, col = epidemic)) +
  geom_line(aes(group = 1), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 years") +
  labs(x = "", y = "influenza case counts",
       title = paste("WHO FluNet data with epidemic periods", sep = " ")) +
  scale_color_manual(values = c("#6e6868", "#e64040")) + 
  facet_wrap(facets = ~Country, scales = "free_y")
#ggsave(filename = "FluNet epidemics four countries.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet")


HM_plot <- filter(HM_epidemic, country %in% c("Bulgaria", "India", "Nigeria", "United States"))
ggplot(HM_plot, aes(x = date, y = counts, col = epidemic)) +
  geom_line(aes(group = 1), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 years") +
  labs(x = "", y = "HealthMap event counts",
       title = "HealthMap data with epidemic periods") +
  scale_color_manual(values = c("#6e6868", "#e64040")) + 
  facet_wrap(facets = ~country, scales = "free_y")
#ggsave(filename = "HM epidemics four countries.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet")

EIOS_plot <- filter(EIOS_epidemic, country %in% c("Bulgaria", "India", "Nigeria", "United States"))
ggplot(EIOS_plot, aes(x = date, y = counts, col = epidemic)) +
  geom_line(aes(group = 1), size = 0.75) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 years") +
  labs(x = "", y = "EIOS event counts",
       title = "EIOS data with epidemic periods") +
  scale_color_manual(values = c("#6e6868", "#e64040")) + 
  facet_wrap(facets = ~country, scales = "free_y")
#ggsave(filename = "EIOS epidemics four countries.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/FluNet")



# outbreak periods (from one start to next start)
HM_epidemic <- HM_epidemic %>% group_by(country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
EIOS_epidemic <- EIOS_epidemic %>% group_by(country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
FluNet_epidemic <- FluNet_epidemic %>% group_by(Country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)



# calculate summary information on outbreaks
FluNet_outbreak_length <- FluNet_epidemic %>% filter(epidemic == TRUE) %>% group_by(Country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n(), start_date = min(SDATE), end_date = max(SDATE), height = max(ALL_INF)) %>% 
  mutate(outbreak_no = row_number()) %>% select(-grp)
FluNet_outbreak_length$outbreak_interval <- interval(FluNet_outbreak_length$start_date, FluNet_outbreak_length$end_date)

HM_outbreak_length <- HM_epidemic %>% filter(epidemic == TRUE) %>% group_by(country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n(), start_date = min(date), end_date = max(date), height = max(counts)) %>% 
  mutate(outbreak_no = row_number()) %>% select(-grp)
HM_outbreak_length$outbreak_interval <- interval(HM_outbreak_length$start_date, HM_outbreak_length$end_date)

EIOS_outbreak_length <- EIOS_epidemic %>% filter(epidemic == TRUE) %>% group_by(country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n(), start_date = min(date), end_date = max(date), height = max(counts)) %>% 
  mutate(outbreak_no = row_number()) %>% select(-grp)
EIOS_outbreak_length$outbreak_interval <- interval(EIOS_outbreak_length$start_date, EIOS_outbreak_length$end_date)


FluNet_outbreaks_for_HM <- filter(FluNet_outbreak_length, start_date %within% interval(min(HM_epidemic$date), max(HM_epidemic$date))) %>% 
  group_by(Country) %>% summarize(no_FluNet_for_HM = n())

FluNet_outbreaks_for_EIOS <- filter(FluNet_outbreak_length, start_date %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date))) %>% 
  group_by(Country) %>% summarize(no_FluNet_for_EIOS = n())



##### calculate sensitivity per outbreak for HM and EIOS #####
# input truth: outbreak intervals from FluNet data for every country
# input source: start and end dates from EIOS or HM data for every country

outbreaks_detected <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                               start_date %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                             start_date %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp_ol <- filter(HM_outbreak_length, country == country_list[i])
  EIOS_temp_ol <- filter(EIOS_outbreak_length, country == country_list[i])
 
  sum_HM <- NA
  for(j in 1:nrow(HM_temp_ol)){
    sum_HM[j] <- sum(int_overlaps(HM_temp_ol$outbreak_interval[j], FluNet_HM_temp_ol$outbreak_interval))
    no_detected_outbreaks_HM <- sum(sum_HM)
  }
  sum_EIOS <- NA
  for(k in 1:nrow(EIOS_temp_ol)){
    sum_EIOS[k] <- sum(int_overlaps(EIOS_temp_ol$outbreak_interval[k], FluNet_EIOS_temp_ol$outbreak_interval))
    no_detected_outbreaks_EIOS <- sum(sum_EIOS)
  } 
  
  outbreaks_detected$country[i] <- country_list[i]
  outbreaks_detected$no_HM[i] <- no_detected_outbreaks_HM
  outbreaks_detected$no_EIOS[i] <- no_detected_outbreaks_EIOS
}



sens_per_outbreak <- data.frame(outbreaks_detected$country, FluNet_outbreaks_for_HM$no_FluNet_for_HM, outbreaks_detected$no_HM,
                           FluNet_outbreaks_for_EIOS$no_FluNet_for_EIOS, outbreaks_detected$no_EIOS)
names(sens_per_outbreak) <- c("country", "no_FluNet_for_HM", "no_HM", "no_FluNet_for_EIOS", "no_EIOS")

sens_per_outbreak <- sens_per_outbreak %>% mutate(sens_HM = no_HM/no_FluNet_for_HM, sens_EIOS = no_EIOS/no_FluNet_for_EIOS)

#write.csv(no_outbreaks, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/sensitivity.csv")


# prepare data for plotting and plot
sens_per_outbreak_long <- pivot_longer(sens_per_outbreak, cols = c(sens_HM, sens_EIOS), names_to = "source", values_to = "sensitivity")
sens_per_outbreak_long$source <- factor(sens_per_outbreak_long$source, levels = c("sens_HM", "sens_EIOS"), labels = c("HealthMap", "EIOS"))
sens_per_outbreak_long$sensitivity <- ifelse(sens_per_outbreak_long$sensitivity > 1, 1, sens_per_outbreak_long$sensitivity)

ggplot(sens_per_outbreak_long, aes(x = country, y = sensitivity * 100, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Sensitivity per outbreak of HealthMap and EIOS systems", y = "Sensitivity per outbreak (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(sens_per_outbreak_long$country))) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set1")
#ggsave(filename = "all_countries_sens_per_outbreak.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

binom.confint(x = sens_per_outbreak$no_HM, n = sens_per_outbreak$no_FluNet_for_HM, method = "exact")
binom.confint(x = sens_per_outbreak$no_EIOS, n = sens_per_outbreak$no_FluNet_for_EIOS, method = "exact")


##### calculate sensitivity per week #####
sens_per_week_list <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  sens_HM <- sum(alarm_HM[state_HM == TRUE], na.rm = TRUE)
  sens_HM <- sens_HM / length(alarm_HM[state_HM == TRUE])
  sens_HM_CI <- binom.confint(sum(alarm_HM[state_HM == TRUE], na.rm = TRUE), length(alarm_HM[state_HM == TRUE]), method = "exact")
  
  state_EIOS <- FluNet_EIOS_temp$epidemic
  alarm_EIOS <- EIOS_temp$epidemic
  sens_EIOS <- sum(alarm_EIOS[state_EIOS == TRUE], na.rm = TRUE)
  sens_EIOS <- sens_EIOS / length(alarm_EIOS[state_EIOS == TRUE])
  sens_EIOS_CI <- binom.confint(sum(alarm_EIOS[state_EIOS == TRUE], na.rm = TRUE), length(alarm_EIOS[state_EIOS == TRUE]), method = "exact")
  
  sens_per_week_list$country[i] <- country_list[i]
  sens_per_week_list$sens_HM[i] <- sens_HM
  sens_per_week_list$sens_EIOS[i] <- sens_EIOS
  sens_per_week_list$CI_HM_lower[i] <- sens_HM_CI$lower
  sens_per_week_list$CI_HM_upper[i] <- sens_HM_CI$upper
  sens_per_week_list$CI_EIOS_lower[i] <- sens_EIOS_CI$lower
  sens_per_week_list$CI_EIOS_upper[i] <- sens_EIOS_CI$upper
}

sens_per_week <- data.frame(sens_per_week_list$country, sens_per_week_list$sens_HM, sens_per_week_list$sens_EIOS)
names(sens_per_week) <- c("country", "sens_HM", "sens_EIOS")

sens_per_week_long <- pivot_longer(sens_per_week, cols = c(sens_HM, sens_EIOS), names_to = "source", 
                                   values_to = "sensitivity_per_week")
sens_per_week_long$source <- factor(sens_per_week_long$source, levels = c("sens_HM", "sens_EIOS"), labels = c("HealthMap", "EIOS"))

ggplot(sens_per_week_long, aes(x = country, y = sensitivity_per_week * 100, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Sensitivity per week of HealthMap and EIOS systems", y = "Sensitivity per week (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(sens_per_week_long$country))) +
  scale_y_continuous(limits = c(0, 100)) + 
  scale_fill_brewer(palette = "Set1")
ggsave(filename = "all_countries_sens_per_week.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

##### calculate exact sensitivity (+/- 1 week) #####
detected_exact <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                                  start_date %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                                start_date %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp_ol <- filter(HM_outbreak_length, country == country_list[i])
  EIOS_temp_ol <- filter(EIOS_outbreak_length, country == country_list[i])
  
  sum_HM <- NA
  for(j in 1:nrow(HM_temp_ol)){
    sum_HM[j] <- sum(HM_temp_ol$start_date[j] %within%
                                  interval(FluNet_HM_temp_ol$start_date[j] - weeks(2), FluNet_HM_temp_ol$start_date[j] + weeks(2)),
                     na.rm = TRUE)
    no_detected_outbreaks_HM <- sum(sum_HM)
  }
  sum_EIOS <- NA
  for(k in 1:nrow(EIOS_temp_ol)){
    sum_EIOS[k] <- sum(EIOS_temp_ol$start_date[k]  %within%
                                    interval(FluNet_EIOS_temp_ol$start_date[k] - weeks(2), FluNet_EIOS_temp_ol$start_date[k] + weeks(2)),
                       na.rm = TRUE)
    no_detected_outbreaks_EIOS <- sum(sum_EIOS)
  } 
  
  detected_exact$country[i] <- country_list[i]
  detected_exact$no_HM[i] <- no_detected_outbreaks_HM
  detected_exact$no_EIOS[i] <- no_detected_outbreaks_EIOS
}

sens_exact <- data.frame(detected_exact$country, FluNet_outbreaks_for_HM$no_FluNet_for_HM, detected_exact$no_HM,
                                FluNet_outbreaks_for_EIOS$no_FluNet_for_EIOS, detected_exact$no_EIOS)
names(sens_exact) <- c("country", "no_FluNet_for_HM", "no_HM", "no_FluNet_for_EIOS", "no_EIOS")
sens_exact <- sens_exact %>% mutate(sens_HM = no_HM/no_FluNet_for_HM, sens_EIOS = no_EIOS/no_FluNet_for_EIOS)

sens_exact_long <- pivot_longer(sens_exact, cols = c(sens_HM, sens_EIOS), names_to = "source", values_to = "exact_sensitivity")
sens_exact_long$source <- factor(sens_exact_long$source, levels = c("sens_HM", "sens_EIOS"), labels = c("HealthMap", "EIOS"))

ggplot(sens_exact_long, aes(x = country, y = exact_sensitivity * 100, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Exact outbreak sensitivity (detection within +/- two weeks)", y = "Exact outbreak sensitivity (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(sens_exact_long$country))) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set1")
ggsave(filename = "all_countries_sens_exact.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

binom.confint(x = sens_exact$no_HM, n = sens_exact$no_FluNet_for_HM, method = "exact")
binom.confint(x = sens_exact$no_EIOS, n = sens_exact$no_FluNet_for_EIOS, method = "exact")

##### positive predictive value per week #####
PPV_list <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  TP_HM <- sum(alarm_HM[state_HM == TRUE], na.rm = TRUE)
  PPV_HM <- TP_HM / sum(alarm_HM, na.rm = TRUE)
  PPV_HM_CI <- binom.confint(TP_HM, sum(alarm_HM, na.rm = TRUE), methods = "exact")
  
  state_EIOS <- FluNet_EIOS_temp$epidemic
  alarm_EIOS <- EIOS_temp$epidemic
  TP_EIOS <- sum(alarm_EIOS[state_EIOS == TRUE], na.rm = TRUE)
  PPV_EIOS <- TP_EIOS / sum(alarm_EIOS, na.rm = TRUE)
  PPV_EIOS_CI <- binom.confint(TP_EIOS, sum(alarm_EIOS, na.rm = TRUE), methods = "exact")
  
  PPV_list$country[i] <- country_list[i]
  PPV_list$PPV_HM[i] <- PPV_HM
  PPV_list$PPV_EIOS[i] <- PPV_EIOS
  PPV_list$PPV_HM_CI_lower[i] <- PPV_HM_CI$lower
  PPV_list$PPV_EIOS_CI_lower[i] <- PPV_EIOS_CI$lower
  PPV_list$PPV_HM_CI_upper[i] <- PPV_HM_CI$upper
  PPV_list$PPV_EIOS_CI_upper[i] <- PPV_EIOS_CI$upper
}

PPV <- data.frame(PPV_list$country, PPV_list$PPV_HM, PPV_list$PPV_EIOS)
names(PPV) <- c("country", "PPV_HM", "PPV_EIOS")

PPV_long <- pivot_longer(PPV, cols = c(PPV_HM, PPV_EIOS), names_to = "source", values_to = "PPV")
PPV_long$source <- factor(PPV_long$source, levels = c("PPV_HM", "PPV_EIOS"), labels = c("HealthMap", "EIOS"))

ggplot(PPV_long, aes(x = country, y = PPV, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Positive predictive value of HealthMap and EIOS systems", y = "Positive predictive value", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(PPV_long$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
ggsave(filename = "all_countries_PPV.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")



NPV_list <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  TN_HM <- sum(alarm_HM[state_HM == FALSE] == FALSE, na.rm = TRUE)
  NPV_HM <- TN_HM / sum(alarm_HM == FALSE, na.rm = TRUE)
  NPV_HM_CI <- binom.confint(TN_HM, sum(alarm_HM == FALSE, na.rm = TRUE), methods = "exact")
  
  state_EIOS <- FluNet_EIOS_temp$epidemic
  alarm_EIOS <- EIOS_temp$epidemic
  TN_EIOS <- sum(alarm_EIOS[state_EIOS == FALSE] == FALSE, na.rm = TRUE)
  NPV_EIOS <- TN_EIOS / sum(alarm_EIOS == FALSE, na.rm = TRUE)
  NPV_EIOS_CI <- binom.confint(TN_EIOS, sum(alarm_EIOS == FALSE, na.rm = TRUE), methods = "exact")
  
  
  NPV_list$country[i] <- country_list[i]
  NPV_list$NPV_HM[i] <- NPV_HM
  NPV_list$NPV_EIOS[i] <- NPV_EIOS
  NPV_list$NPV_HM_CI_lower[i] <- NPV_HM_CI$lower
  NPV_list$NPV_EIOS_CI_lower[i] <- NPV_EIOS_CI$lower
  NPV_list$NPV_HM_CI_upper[i] <- NPV_HM_CI$upper
  NPV_list$NPV_EIOS_CI_upper[i] <- NPV_EIOS_CI$upper
}


NPV <- data.frame(NPV_list$country, NPV_list$NPV_HM, NPV_list$NPV_EIOS)
names(NPV) <- c("country", "NPV_HM", "NPV_EIOS")

NPV_long <- pivot_longer(NPV, cols = c(NPV_HM, NPV_EIOS), names_to = "source", values_to = "NPV")
NPV_long$source <- factor(NPV_long$source, levels = c("NPV_HM", "NPV_EIOS"), labels = c("HealthMap", "EIOS"))

ggplot(NPV_long, aes(x = country, y = NPV, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Negative predictive value of HealthMap and EIOS systems", y = "Negative predictive value", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(NPV_long$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
ggsave(filename = "all_countries_PPV.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")



##### calculate false alarm rate #####
# calculate false alert rate during non-outbreak intervals

false_alarm_rate <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  FA_HM <- sum(alarm_HM[state_HM == FALSE], na.rm = TRUE)
  FAR_HM <- FA_HM / length(alarm_HM[state_HM == FALSE])
  FAR_HM_CI <- binom.confint(FAR_HM, length(alarm_HM[state_HM == FALSE]), methods = "exact")
  
  state_EIOS <- FluNet_EIOS_temp$epidemic
  alarm_EIOS <- EIOS_temp$epidemic
  FA_EIOS <- sum(alarm_EIOS[state_EIOS == FALSE], na.rm = TRUE)
  FAR_EIOS <- FA_EIOS / length(alarm_EIOS[state_EIOS == FALSE])
  FAR_EIOS_CI <- binom.confint(FAR_EIOS, length(alarm_EIOS[state_EIOS == FALSE]), methods = "exact")
  
  false_alarm_rate$country[i] <- country_list[i]
  false_alarm_rate$FAR_HM[i] <- FAR_HM
  false_alarm_rate$FAR_EIOS[i] <- FAR_EIOS
  false_alarm_rate$FAR_HM_CI_lower[i] <- FAR_HM_CI$lower
  false_alarm_rate$FAR_HM_CI_upper[i] <- FAR_HM_CI$upper
  false_alarm_rate$FAR_EIOS_CI_lower[i] <- FAR_EIOS_CI$lower
  false_alarm_rate$FAR_EIOS_CI_upper[i] <- FAR_EIOS_CI$upper
}

FAR <- data.frame(false_alarm_rate$country, false_alarm_rate$FAR_HM, false_alarm_rate$FAR_EIOS)
names(FAR) <- c("country", "FAR_HM", "FAR_EIOS")

FAR_long <- pivot_longer(FAR, cols = c(FAR_HM, FAR_EIOS), names_to = "source", values_to = "false_alarm_rate")
FAR_long$source <- factor(FAR_long$source, levels = c("FAR_HM", "FAR_EIOS"), labels = c("HealthMap", "EIOS"))
FAR_long$specificity <- 1 - FAR_long$false_alarm_rate

ggplot(FAR_long, aes(x = country, y = false_alarm_rate, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "False alarm rate of HealthMap and EIOS systems", y = "False alarm rate", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(FAR_long$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
ggsave(filename = "all_countries_false_alarm_rate.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

ggplot(FAR_long, aes(x = country, y = specificity * 100, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Specificity of HealthMap and EIOS systems", y = "Specificity (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(FAR_long$country))) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set1")
ggsave(filename = "all_countries_specificity.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")


##### calculate timeliness#####
#  input df is data frame with 3 col: date, baseline, baseline + outbreak
#  input truth is list with element: state, which is binary vector indicating presence of outbreak
#  output is proportion of outbreak "prevented" (from 1 for first day to 0 for last day of outbreak, or later)


timeliness_HM <- vector(mode = "list")
timeliness_EIOS <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  prevented <- 0
  for(j in 1:max(FluNet_HM_temp$outbreak_period)){
    FluNet_HM_temp_ob <- filter(FluNet_HM_temp, outbreak_period == j)
    HM_temp <- filter(HM_epidemic, country == country_list[i] & 
                        date %within% interval(min(FluNet_HM_temp_ob$SDATE), max(FluNet_HM_temp_ob$SDATE)))
    
    state <- FluNet_HM_temp_ob$epidemic
    alarm <- HM_temp$epidemic
    length <- sum(FluNet_HM_temp_ob$epidemic)
    
    detect <- ifelse(sum(alarm[state == TRUE], na.rm = TRUE) > 0, TRUE, FALSE)
    if (detect) {
      first.alarm <- min(which(alarm[state == TRUE] == TRUE))
      prevented[j] <- (length - first.alarm) / length
    }else{
      prevented[j] <- 0
    }
  }
  timeliness_HM$country[i] <- country_list[i]
  timeliness_HM$timeliness[i] <- mean(prevented, na.rm = TRUE)
  timeliness_HM$time_CI_lower[i] <- mean(prevented, na.rm = TRUE) - 1.96*sd(prevented)/sqrt(length(prevented))
  timeliness_HM$time_CI_upper[i] <- mean(prevented, na.rm = TRUE) + 1.96*sd(prevented)/sqrt(length(prevented))
  
}  

for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_EIOS_temp <- FluNet_EIOS_temp %>% group_by(Country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
  
  prevented <- 0
  for(j in 1:max(FluNet_EIOS_temp$outbreak_period)){
    FluNet_EIOS_temp_ob <- filter(FluNet_EIOS_temp, outbreak_period == j)
    EIOS_temp <- filter(EIOS_epidemic, country == country_list[i] & 
                        date %within% interval(min(FluNet_EIOS_temp_ob$SDATE), max(FluNet_EIOS_temp_ob$SDATE)))
    
    state <- FluNet_EIOS_temp_ob$epidemic
    alarm <- EIOS_temp$epidemic
    length <- sum(FluNet_EIOS_temp_ob$epidemic)
    
    detect <- ifelse(sum(alarm[state == TRUE], na.rm = TRUE) > 0, TRUE, FALSE)
    if (detect == TRUE){
      first.alarm <- min(which(alarm[state == TRUE] == TRUE))
      prevented[j] <- (length - first.alarm) / length
    }else{
      prevented[j] <- 0
    }
  }
  timeliness_EIOS$country[i] <- country_list[i]
  timeliness_EIOS$timeliness[i] <- mean(prevented, na.rm = TRUE)
  timeliness_EIOS$time_CI_lower <- mean(prevented, na.rm = TRUE) - 1.96*sd(prevented)/sqrt(length(prevented))
  timeliness_EIOS$time_CI_upper <- mean(prevented, na.rm = TRUE) + 1.96*sd(prevented)/sqrt(length(prevented))
  
} 

timeliness <- data.frame(timeliness_HM$country, timeliness_HM$timeliness, timeliness_EIOS$timeliness)
names(timeliness) <- c("country", "frac_prevented_HM", "frac_prevented_EIOS")

timeliness_long <- pivot_longer(timeliness, cols = c(frac_prevented_HM, frac_prevented_EIOS), 
                                names_to = "source", values_to = "frac_prevented")
timeliness_long$source <- factor(timeliness_long$source, levels = c("frac_prevented_HM", "frac_prevented_EIOS"), 
                                 labels = c("HealthMap", "EIOS"))


ggplot(timeliness_long, aes(x = country, y = frac_prevented, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Timeliness (prevented fraction) of HealthMap and EIOS systems", y = "prevented fraction", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(timeliness_long$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
ggsave(filename = "all_countries_prevented_frac.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")


mean_time_HM <- vector(mode = "list") 
for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  time_to_detect <- 0
  for(j in 1:max(FluNet_HM_temp$outbreak_period)){
    FluNet_HM_temp_ob <- filter(FluNet_HM_temp, outbreak_period == j)
    HM_temp <- filter(HM_epidemic, country == country_list[i] & 
                        date %within% interval(min(FluNet_HM_temp_ob$SDATE), max(FluNet_HM_temp_ob$SDATE)))
    
    state <- FluNet_HM_temp_ob$epidemic
    alarm <- HM_temp$epidemic
    length <- sum(FluNet_HM_temp_ob$epidemic)
    
    detect <- ifelse(sum(alarm[state == TRUE], na.rm = TRUE) > 0, TRUE, FALSE)
    if (detect) {
      first.alarm_source <- min(which(alarm[state == TRUE] == TRUE))
      first.alarm_truth <- min(which(state == TRUE))
      time_to_detect[j] <- first.alarm_source - first.alarm_truth
    }else{
      time_to_detect[j] <- NA
    }
  }
  mean_time_HM$country[i] <- country_list[i]
  mean_time_HM$time_to_detect[i] <- mean(time_to_detect, na.rm = TRUE)
}

mean_time_EIOS <- vector(mode = "list") 
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  time_to_detect <- 0
  for(j in 1:max(FluNet_EIOS_temp$outbreak_period)){
    FluNet_EIOS_temp_ob <- filter(FluNet_EIOS_temp, outbreak_period == j)
    EIOS_temp <- filter(EIOS_epidemic, country == country_list[i] & 
                        date %within% interval(min(FluNet_EIOS_temp_ob$SDATE), max(FluNet_EIOS_temp_ob$SDATE)))
    
    state <- FluNet_EIOS_temp_ob$epidemic
    alarm <- EIOS_temp$epidemic
    length <- sum(FluNet_EIOS_temp_ob$epidemic)
    
    detect <- ifelse(sum(alarm[state == TRUE], na.rm = TRUE) > 0, TRUE, FALSE)
    if (detect) {
      first.alarm_source <- min(which(alarm[state == TRUE] == TRUE))
      first.alarm_truth <- min(which(state == TRUE))
      time_to_detect[j] <- first.alarm_source - first.alarm_truth
    }else{
      time_to_detect[j] <- NA
    }
  }
  mean_time_EIOS$country[i] <- country_list[i]
  mean_time_EIOS$time_to_detect[i] <- mean(time_to_detect, na.rm = TRUE)
}

mean_time_detect <- data.frame(mean_time_HM$country, mean_time_HM$time_to_detect, mean_time_EIOS$time_to_detect)
names(mean_time_detect) <- c("country", "mean_time_HM", "mean_time_EIOS")

mean_time_detect_long <- pivot_longer(mean_time_detect, cols = c(mean_time_HM, mean_time_EIOS), 
                                names_to = "source", values_to = "mean_time")
mean_time_detect_long$source <- factor(mean_time_detect_long$source, levels = c("mean_time_HM", "mean_time_EIOS"), 
                                 labels = c("HealthMap", "EIOS"))
mean_time_detect_long$mean_time <- ifelse(is.na(mean_time_detect_long$mean_time), 0, mean_time_detect_long$mean_time)


ggplot(mean_time_detect_long, aes(x = country, y = mean_time, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Mean time to outbreak detection of HealthMap and EIOS systems", y = "mean time to outbreak detection (weeks)", x = "", 
       caption = "WHO FluNet data were used as gold standard \n
       Non-detected outbreaks are not counted") +
  scale_x_discrete(limits = rev(levels(mean_time_detect_long$country))) +
  scale_fill_brewer(palette = "Set1")
#ggsave(filename = "all_countries_mean_time.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")


##### accuracy #####
accuracy_list <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  TP_HM <- sum(alarm_HM[state_HM == TRUE], na.rm = TRUE)
  TN_HM <- sum(alarm_HM[state_HM == FALSE] == FALSE, na.rm = TRUE)
  accuracy_HM <- (TP_HM + TN_HM)/nrow(HM_temp)
  accuracy_HM_CI <- binom.confint(TP_HM + TN_HM, nrow(HM_temp), methods = "exact")
  
  state_EIOS <- FluNet_EIOS_temp$epidemic
  alarm_EIOS <- EIOS_temp$epidemic
  TP_EIOS <- sum(alarm_EIOS[state_EIOS == TRUE], na.rm = TRUE)
  TN_EIOS <- sum(alarm_EIOS[state_EIOS == FALSE] == FALSE, na.rm = TRUE)
  accuracy_EIOS <- (TP_EIOS + TN_EIOS)/nrow(EIOS_temp)
  accuracy_EIOS_CI <- binom.confint(TP_EIOS + TN_EIOS, nrow(EIOS_temp), methods = "exact")
  
  accuracy_list$country[i] <- country_list[i]
  accuracy_list$accuracy_HM[i] <- accuracy_HM
  accuracy_list$accuracy_EIOS[i] <- accuracy_EIOS
  accuracy_list$accuracy_HM_CI_lower[i] <- accuracy_HM_CI$lower
  accuracy_list$accuracy_HM_CI_upper[i] <- accuracy_HM_CI$upper
  accuracy_list$accuracy_EIOS_CI_lower[i] <- accuracy_EIOS_CI$lower
  accuracy_list$accuracy_EIOS_CI_upper[i] <- accuracy_EIOS_CI$upper
}

accuracy <- data.frame(accuracy_list$country, accuracy_list$accuracy_HM, accuracy_list$accuracy_EIOS)
names(accuracy) <- c("country", "accuracy_HM", "accuracy_EIOS")

accuracy_long <- pivot_longer(accuracy, cols = c(accuracy_HM, accuracy_EIOS), 
                                      names_to = "source", values_to = "accuracy")
accuracy_long$source <- factor(accuracy_long$source, levels = c("accuracy_HM", "accuracy_EIOS"), 
                                       labels = c("HealthMap", "EIOS"))

ggplot(accuracy_long, aes(x = country, y = accuracy, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Accuracy of HealthMap and EIOS systems", y = "Accuracy", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(accuracy_long$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
ggsave(filename = "all_countries_accuracy.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

##### combine all metrics into one dataframe ##### 
metrics <- data.frame(sens_per_outbreak_long$country, sens_per_outbreak_long$source, sens_per_outbreak_long$sensitivity, 
                      sens_exact_long$exact_sensitivity, sens_per_week_long$sensitivity_per_week, PPV_long$PPV, FAR_long$specificity, timeliness_long$frac_prevented)
names(metrics) <- c("country", "source", "sens_per_outbreak", "sens_exact", "sens_per_week", "PPV", "specificity", "frac_prevented")
#write.csv(metrics, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics.csv")

metrics_long <- pivot_longer(metrics, cols = c(sens_per_outbreak, sens_per_week, sens_exact, PPV, specificity, frac_prevented), 
                             names_to = "metric", values_to = "values")


for(i in seq_along(country_list)){
  metrics_temp <- filter(metrics_long, country == country_list[i])

  plot <- ggplot(metrics_temp, aes(x = metric, y = values, fill = source)) +
    geom_col(position = "dodge") +
    labs(title = paste(country_list[i], ": Evaluation of HealthMap and EIOS systems", sep = ""), y = "", x = "", 
         caption = "WHO FluNet data were used as gold standard") +
    scale_x_discrete(labels = c('Prevented fraction \n (Timeliness)', 'Positive \n predictive value', 'Exact sensitivity \n (+/- 2 weeks)', 
                               'Sensitivity per \n outbreak', 'Sensitivity per \n week', 'Specificity')) +
    scale_y_continuous(limits = c(0, 1)) + 
    scale_fill_brewer(palette = "Set1") 
  print(plot)
  ggsave(plot = plot, file = paste("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/Evaluation metrics", country_list[i], ".jpeg", sep=' '))
}



##### repeat for all_bcp dfs ##### - obsolete 

HM_epidemic_all_bcp <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HealthMap_epidemic_all_bcp.csv")
HM_epidemic_all_bcp$date <- as.POSIXct(HM_epidemic_all_bcp$date)
EIOS_epidemic_all_bcp <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_epidemic_all_bcp.csv")
EIOS_epidemic_all_bcp$date <- as.POSIXct(EIOS_epidemic_all_bcp$date)


HM_epidemic_all_bcp$startend <- NA
for(i in 2:nrow(HM_epidemic_all_bcp)){
  if(HM_epidemic_all_bcp$epidemic[i] == TRUE & HM_epidemic_all_bcp$epidemic[i-1] == FALSE){
    HM_epidemic_all_bcp$startend[i] <- "start"
  }else if(HM_epidemic_all_bcp$epidemic[i] == FALSE & HM_epidemic_all_bcp$epidemic[i-1] == TRUE){
    HM_epidemic_all_bcp$startend[i] <- "end"
  }
}

EIOS_epidemic_all_bcp$startend <- NA
for(i in 2:nrow(EIOS_epidemic_all_bcp)){
  if(EIOS_epidemic_all_bcp$epidemic[i] == TRUE & EIOS_epidemic_all_bcp$epidemic[i-1] == FALSE){
    EIOS_epidemic_all_bcp$startend[i] <- "start"
  }else if(EIOS_epidemic_all_bcp$epidemic[i] == FALSE & EIOS_epidemic_all_bcp$epidemic[i-1] == TRUE){
    EIOS_epidemic_all_bcp$startend[i] <- "end"
  }
}

HM_epidemic_all_bcp <- HM_epidemic_all_bcp %>% group_by(country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
EIOS_epidemic_all_bcp <- EIOS_epidemic_all_bcp %>% group_by(country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)


# calculate summary information on outbreaks
HM_outbreak_length_all_bcp <- HM_epidemic_all_bcp %>% filter(epidemic == TRUE) %>% group_by(country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n(), start_date = min(date), end_date = max(date), height = max(counts)) %>% 
  mutate(outbreak_no = row_number()) %>% select(-grp)
HM_outbreak_length_all_bcp$outbreak_interval <- interval(HM_outbreak_length_all_bcp$start_date, HM_outbreak_length_all_bcp$end_date)

EIOS_outbreak_length_all_bcp <- EIOS_epidemic_all_bcp %>% filter(epidemic == TRUE) %>% group_by(country, grp = cumsum(!is.na(startend))) %>% 
  summarize(length = n(), start_date = min(date), end_date = max(date), height = max(counts)) %>% 
  mutate(outbreak_no = row_number()) %>% select(-grp)
EIOS_outbreak_length_all_bcp$outbreak_interval <- interval(EIOS_outbreak_length_all_bcp$start_date, EIOS_outbreak_length_all_bcp$end_date)


### sensitivity per outbreak
outbreaks_detected_all_bcp <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                                  start_date %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                                start_date %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp_ol <- filter(HM_outbreak_length_all_bcp, country == country_list[i])
  EIOS_temp_ol <- filter(EIOS_outbreak_length_all_bcp, country == country_list[i])
  
  sum_HM <- NA
  for(j in 1:nrow(HM_temp_ol)){
    sum_HM[j] <- sum(int_overlaps(HM_temp_ol$outbreak_interval[j], FluNet_HM_temp_ol$outbreak_interval))
    no_detected_outbreaks_HM <- sum(sum_HM)
  }
  sum_EIOS <- NA
  for(k in 1:nrow(EIOS_temp_ol)){
    sum_EIOS[k] <- sum(int_overlaps(EIOS_temp_ol$outbreak_interval[k], FluNet_EIOS_temp_ol$outbreak_interval))
    no_detected_outbreaks_EIOS <- sum(sum_EIOS)
  } 
  
  outbreaks_detected_all_bcp$country[i] <- country_list[i]
  outbreaks_detected_all_bcp$no_HM[i] <- no_detected_outbreaks_HM
  outbreaks_detected_all_bcp$no_EIOS[i] <- no_detected_outbreaks_EIOS
}


FluNet_outbreaks_for_HM <- filter(FluNet_outbreak_length, start_date %within% interval(min(HM_epidemic$date), max(HM_epidemic$date))) %>% 
  group_by(Country) %>% summarize(no_FluNet_for_HM = n())

FluNet_outbreaks_for_EIOS <- filter(FluNet_outbreak_length, start_date %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date))) %>% 
  group_by(Country) %>% summarize(no_FluNet_for_EIOS = n())


sens_per_outbreak_all_bcp <- data.frame(outbreaks_detected_all_bcp$country, FluNet_outbreaks_for_HM$no_FluNet_for_HM, 
                                        outbreaks_detected_all_bcp$no_HM,
                                FluNet_outbreaks_for_EIOS$no_FluNet_for_EIOS, outbreaks_detected_all_bcp$no_EIOS)
names(sens_per_outbreak_all_bcp) <- c("country", "no_FluNet_for_HM", "no_HM", "no_FluNet_for_EIOS", "no_EIOS")

sens_per_outbreak_all_bcp <- sens_per_outbreak_all_bcp %>% 
  mutate(sens_HM_all_bcp = no_HM/no_FluNet_for_HM, sens_EIOS_all_bcp = no_EIOS/no_FluNet_for_EIOS)

# prepare data for plotting and plot
sens_per_outbreak_long_all_bcp <- pivot_longer(sens_per_outbreak_all_bcp, cols = c(sens_HM_all_bcp, sens_EIOS_all_bcp), names_to = "source", values_to = "sensitivity")
sens_per_outbreak_long_all_bcp$source <- factor(sens_per_outbreak_long_all_bcp$source, levels = c("sens_HM_all_bcp", "sens_EIOS_all_bcp"), labels = c("HealthMap", "EIOS"))
sens_per_outbreak_long_all_bcp$sensitivity <- ifelse(sens_per_outbreak_long_all_bcp$sensitivity > 1, 1, sens_per_outbreak_long_all_bcp$sensitivity)

ggplot(sens_per_outbreak_long_all_bcp, aes(x = country, y = sensitivity * 100, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Sensitivity per outbreak of HealthMap and EIOS systems (all bcp)", y = "Sensitivity per outbreak (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(sens_per_outbreak_long_all_bcp$country))) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set1")
#ggsave(filename = "all_countries_sens_per_outbreak_all_bcp.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")


##### calculate sensitivity per week #####
sens_per_week_list_all_bcp <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic_all_bcp, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic_all_bcp, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  sens_HM <- sum(alarm_HM[state_HM == TRUE], na.rm = TRUE)
  sens_HM <- sens_HM / length(alarm_HM[state_HM == TRUE])
  
  state_EIOS <- FluNet_EIOS_temp$epidemic
  alarm_EIOS <- EIOS_temp$epidemic
  sens_EIOS <- sum(alarm_EIOS[state_EIOS == TRUE], na.rm = TRUE)
  sens_EIOS <- sens_EIOS / length(alarm_EIOS[state_EIOS == TRUE])
  
  sens_per_week_list_all_bcp$country[i] <- country_list[i]
  sens_per_week_list_all_bcp$sens_HM[i] <- sens_HM
  sens_per_week_list_all_bcp$sens_EIOS[i] <- sens_EIOS
}

sens_per_week_all_bcp <- data.frame(sens_per_week_list_all_bcp$country, sens_per_week_list_all_bcp$sens_HM, sens_per_week_list_all_bcp$sens_EIOS)
names(sens_per_week_all_bcp) <- c("country", "sens_HM", "sens_EIOS")

sens_per_week_long_all_bcp <- pivot_longer(sens_per_week_all_bcp, cols = c(sens_HM, sens_EIOS), names_to = "source", 
                                   values_to = "sensitivity_per_week")
sens_per_week_long_all_bcp$source <- factor(sens_per_week_long_all_bcp$source, levels = c("sens_HM", "sens_EIOS"), labels = c("HealthMap", "EIOS"))

ggplot(sens_per_week_long_all_bcp, aes(x = country, y = sensitivity_per_week * 100, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Sensitivity per week of HealthMap and EIOS systems (all bcp)", y = "Sensitivity per week (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(sens_per_week_long_all_bcp$country))) +
  scale_y_continuous(limits = c(0, 100)) + 
  scale_fill_brewer(palette = "Set1")
#ggsave(filename = "all_countries_sens_per_week_all_bcp.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

##### calculate exact sensitivity (+/- one week) #####
detected_exact_all_bcp <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                                  start_date %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp_ol <- filter(FluNet_outbreak_length, Country == country_list[i] & 
                                start_date %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp_ol <- filter(HM_outbreak_length_all_bcp, country == country_list[i])
  EIOS_temp_ol <- filter(EIOS_outbreak_length_all_bcp, country == country_list[i])
  
  sum_HM <- NA
  for(j in 1:nrow(HM_temp_ol)){
    sum_HM[j] <- sum(HM_temp_ol$start_date[j] %within%
                       interval(FluNet_HM_temp_ol$start_date - weeks(1), FluNet_HM_temp_ol$start_date + weeks(1)),
                     na.rm = TRUE)
    no_detected_outbreaks_HM <- sum(sum_HM)
  }
  sum_EIOS <- NA
  for(k in 1:nrow(EIOS_temp_ol)){
    sum_EIOS[k] <- sum(EIOS_temp_ol$start_date[k]  %within%
                         interval(FluNet_EIOS_temp_ol$start_date[k] - weeks(1), FluNet_EIOS_temp_ol$start_date[k] + weeks(1)),
                       na.rm = TRUE)
    no_detected_outbreaks_EIOS <- sum(sum_EIOS)
  } 
  
  detected_exact_all_bcp$country[i] <- country_list[i]
  detected_exact_all_bcp$no_HM[i] <- no_detected_outbreaks_HM
  detected_exact_all_bcp$no_EIOS[i] <- no_detected_outbreaks_EIOS
}

sens_exact_all_bcp <- data.frame(detected_exact_all_bcp$country, FluNet_outbreaks_for_HM$no_FluNet_for_HM, detected_exact_all_bcp$no_HM,
                         FluNet_outbreaks_for_EIOS$no_FluNet_for_EIOS, detected_exact_all_bcp$no_EIOS)
names(sens_exact_all_bcp) <- c("country", "no_FluNet_for_HM", "no_HM", "no_FluNet_for_EIOS", "no_EIOS")
sens_exact_all_bcp <- sens_exact_all_bcp %>% mutate(sens_HM = no_HM/no_FluNet_for_HM, sens_EIOS = no_EIOS/no_FluNet_for_EIOS)

sens_exact_long_all_bcp <- pivot_longer(sens_exact_all_bcp, cols = c(sens_HM, sens_EIOS), names_to = "source", values_to = "exact_sensitivity")
sens_exact_long_all_bcp$source <- factor(sens_exact_long_all_bcp$source, levels = c("sens_HM", "sens_EIOS"), labels = c("HealthMap", "EIOS"))

ggplot(sens_exact_long_all_bcp, aes(x = country, y = exact_sensitivity * 100, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Exact outbreak sensitivity (detection within +/- one week, all bcp)", y = "Exact outbreak sensitivity (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(sens_exact_long_all_bcp$country))) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set1")
#ggsave(filename = "all_countries_sens_exact_all_bcp.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")


##### positive predictive value per week #####
PPV_list_all_bcp <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic_all_bcp, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic_all_bcp, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  TP_HM <- sum(alarm_HM[state_HM == TRUE], na.rm = TRUE)
  PPV_HM <- TP_HM / sum(alarm_HM, na.rm = TRUE)
  
  state_EIOS <- FluNet_EIOS_temp$epidemic
  alarm_EIOS <- EIOS_temp$epidemic
  TP_EIOS <- sum(alarm_EIOS[state_EIOS == TRUE], na.rm = TRUE)
  PPV_EIOS <- TP_EIOS / sum(alarm_EIOS, na.rm = TRUE)
  
  PPV_list_all_bcp$country[i] <- country_list[i]
  PPV_list_all_bcp$PPV_HM[i] <- PPV_HM
  PPV_list_all_bcp$PPV_EIOS[i] <- PPV_EIOS
}

PPV_all_bcp <- data.frame(PPV_list_all_bcp$country, PPV_list_all_bcp$PPV_HM, PPV_list_all_bcp$PPV_EIOS)
names(PPV_all_bcp) <- c("country", "PPV_HM", "PPV_EIOS")

PPV_long_all_bcp <- pivot_longer(PPV_all_bcp, cols = c(PPV_HM, PPV_EIOS), names_to = "source", values_to = "PPV")
PPV_long_all_bcp$source <- factor(PPV_long_all_bcp$source, levels = c("PPV_HM", "PPV_EIOS"), labels = c("HealthMap", "EIOS"))

ggplot(PPV_long_all_bcp, aes(x = country, y = PPV, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Positive predictive value of HealthMap and EIOS systems (all bcp)", y = "Positive predictive value", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(PPV_long_all_bcp$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
#ggsave(filename = "all_countries_PPV_all_bcp.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")



##### calculate false alarm rate #####
# calculate false alert rate during non-outbreak intervals

false_alarm_rate_all_bcp <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  HM_temp <- filter(HM_epidemic_all_bcp, country == country_list[i])
  EIOS_temp <- filter(EIOS_epidemic_all_bcp, country == country_list[i])
  
  state_HM <- FluNet_HM_temp$epidemic
  alarm_HM <- HM_temp$epidemic
  FA_HM <- sum(alarm_HM[state_HM == FALSE], na.rm = TRUE)
  FAR_HM <- FA_HM / length(alarm_HM[state_HM == FALSE])
  
  state_EIOS <- FluNet_EIOS_temp$epidemic
  alarm_EIOS <- EIOS_temp$epidemic
  FA_EIOS <- sum(alarm_EIOS[state_EIOS == FALSE], na.rm = TRUE)
  FAR_EIOS <- FA_EIOS / length(alarm_EIOS[state_EIOS == FALSE])
  
  false_alarm_rate_all_bcp$country[i] <- country_list[i]
  false_alarm_rate_all_bcp$FAR_HM[i] <- FAR_HM
  false_alarm_rate_all_bcp$FAR_EIOS[i] <- FAR_EIOS
}

FAR_all_bcp <- data.frame(false_alarm_rate_all_bcp$country, false_alarm_rate_all_bcp$FAR_HM, false_alarm_rate_all_bcp$FAR_EIOS)
names(FAR_all_bcp) <- c("country", "FAR_HM", "FAR_EIOS")

FAR_long_all_bcp <- pivot_longer(FAR_all_bcp, cols = c(FAR_HM, FAR_EIOS), names_to = "source", values_to = "false_alarm_rate")
FAR_long_all_bcp$source <- factor(FAR_long_all_bcp$source, levels = c("FAR_HM", "FAR_EIOS"), labels = c("HealthMap", "EIOS"))
FAR_long_all_bcp$specificity <- 1 - FAR_long_all_bcp$false_alarm_rate

ggplot(FAR_long_all_bcp, aes(x = country, y = false_alarm_rate, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "False alarm rate of HealthMap and EIOS systems (all bcp)", y = "False alarm rate", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(FAR_long_all_bcp$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
#ggsave(filename = "all_countries_false_alarm_rate_all_bcp.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")

ggplot(FAR_long_all_bcp, aes(x = country, y = specificity * 100, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Specificity of HealthMap and EIOS systems (all bcp)", y = "Specificity (%)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(FAR_long_all_bcp$country))) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set1")
#ggsave(filename = "all_countries_specificity_all_bcp.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")


##### calculate timeliness#####
#  input df is data frame with 3 col: date, baseline, baseline + outbreak
#  input truth is list with element: state, which is binary vector indicating presence of outbreak
#  output is proportion of outbreak "prevented" (from 1 for first day to 0 for last day of outbreak, or later)


timeliness_HM_all_bcp <- vector(mode = "list")
timeliness_EIOS_all_bcp <- vector(mode = "list")
for(i in seq_along(country_list)){
  FluNet_HM_temp <- filter(FluNet_epidemic, Country == country_list[i] &
                             SDATE %within% interval(min(HM_epidemic$date), max(HM_epidemic$date)))
  prevented <- 0
  for(j in 1:max(FluNet_HM_temp$outbreak_period)){
    FluNet_HM_temp_ob <- filter(FluNet_HM_temp, outbreak_period == j)
    HM_temp <- filter(HM_epidemic_all_bcp, country == country_list[i] & 
                        date %within% interval(min(FluNet_HM_temp_ob$SDATE), max(FluNet_HM_temp_ob$SDATE)))
    
    state <- FluNet_HM_temp_ob$epidemic
    alarm <- HM_temp$epidemic
    length <- sum(FluNet_HM_temp_ob$epidemic)
    
    detect <- ifelse(sum(alarm[state == TRUE], na.rm = TRUE) > 0, TRUE, FALSE)
    if (detect) {
      first.alarm <- min(which(alarm[state == TRUE] == TRUE))
      prevented[j] <- (length - first.alarm) / length
    }else{
      prevented[j] <- 0
    }
  }
  timeliness_HM_all_bcp$country[i] <- country_list[i]
  timeliness_HM_all_bcp$timeliness[i] <- mean(prevented, na.rm = TRUE)
}  

for(i in seq_along(country_list)){
  FluNet_EIOS_temp <- filter(FluNet_epidemic, Country == country_list[i] & 
                               SDATE %within% interval(min(EIOS_epidemic$date), max(EIOS_epidemic$date)))
  FluNet_EIOS_temp <- FluNet_EIOS_temp %>% group_by(Country) %>% mutate(outbreak_period = ((cumsum(!is.na(startend))-1) %/% 2) + 1)
  
  prevented <- 0
  for(j in 1:max(FluNet_EIOS_temp$outbreak_period)){
    FluNet_EIOS_temp_ob <- filter(FluNet_EIOS_temp, outbreak_period == j)
    EIOS_temp <- filter(EIOS_epidemic_all_bcp, country == country_list[i] & 
                          date %within% interval(min(FluNet_EIOS_temp_ob$SDATE), max(FluNet_EIOS_temp_ob$SDATE)))
    
    state <- FluNet_EIOS_temp_ob$epidemic
    alarm <- EIOS_temp$epidemic
    length <- sum(FluNet_EIOS_temp_ob$epidemic)
    
    detect <- ifelse(sum(alarm[state == TRUE], na.rm = TRUE) > 0, TRUE, FALSE)
    if (detect == TRUE){
      first.alarm <- min(which(alarm[state == TRUE] == TRUE))
      prevented[j] <- (length - first.alarm) / length
    }else{
      prevented[j] <- 0
    }
  }
  timeliness_EIOS_all_bcp$country[i] <- country_list[i]
  timeliness_EIOS_all_bcp$timeliness[i] <- mean(prevented, na.rm = TRUE)
} 

timeliness_all_bcp <- data.frame(timeliness_HM_all_bcp$country, timeliness_HM_all_bcp$timeliness, timeliness_EIOS_all_bcp$timeliness)
names(timeliness_all_bcp) <- c("country", "frac_prevented_HM", "frac_prevented_EIOS")

timeliness_long_all_bcp <- pivot_longer(timeliness_all_bcp, cols = c(frac_prevented_HM, frac_prevented_EIOS), 
                                names_to = "source", values_to = "frac_prevented")
timeliness_long_all_bcp$source <- factor(timeliness_long_all_bcp$source, levels = c("frac_prevented_HM", "frac_prevented_EIOS"), 
                                 labels = c("HealthMap", "EIOS"))


ggplot(timeliness_long_all_bcp, aes(x = country, y = frac_prevented, fill = source)) +
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Timeliness (prevented fraction) of HealthMap and EIOS systems (all bcp)", y = "prevented fraction", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(timeliness_long_all_bcp$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
#ggsave(filename = "all_countries_prevented_frac_all_bcp.jpeg", path = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic")


### combine all metrics into one dataframe ###
metrics_all_bcp <- data.frame(sens_per_outbreak_long_all_bcp$country, sens_per_outbreak_long_all_bcp$source, 
                              sens_per_outbreak_long_all_bcp$sensitivity, sens_exact_long_all_bcp$exact_sensitivity, 
                              sens_per_week_long_all_bcp$sensitivity_per_week, PPV_long_all_bcp$PPV, 
                              FAR_long_all_bcp$specificity, timeliness_long_all_bcp$frac_prevented)
names(metrics_all_bcp) <- c("country", "source", "sens_per_outbreak", "sens_exact", "sens_per_week", "PPV", "specificity", "frac_prevented")
#write.csv(metrics_all_bcp, file = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics_all_bcp.csv")

metrics_long_all_bcp <- pivot_longer(metrics_all_bcp, cols = c(sens_per_outbreak, sens_per_week, sens_exact, PPV, specificity, frac_prevented), 
                             names_to = "metric", values_to = "values")



##### compare all_bcp and multiple methods evaluation metrics#####
metrics <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics.csv", stringsAsFactors = TRUE)
metrics$method <- "multiple"
metrics <- metrics %>% unite("source_method", c(source, method), remove = FALSE)
metrics <- metrics %>% select(-X)

metrics_all_bcp <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics_all_bcp.csv")
metrics_all_bcp$method <- "bcp"
metrics_all_bcp <- metrics_all_bcp %>% unite("source_method", c(source, method), remove = FALSE)
metrics_all_bcp <- metrics_all_bcp %>% select(-X)

 
comparison_df <- rbind(metrics, metrics_all_bcp)

ggplot(filter(comparison_df, source == "HealthMap"), aes(x = country, y = frac_prevented, fill = method)) + 
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Comparison of all bcp and multiple methods for outbreak detection in HealthMap", y = "prevented fraction", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(comparison_df$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")

ggplot(filter(comparison_df, source == "EIOS"), aes(x = country, y = frac_prevented, fill = method)) + 
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Comparison of all bcp and multiple methods for outbreak detection in EIOS", y = "prevented fraction", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(comparison_df$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")

ggplot(comparison_df, aes(x = country, y = frac_prevented, fill = source_method)) + 
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Comparison of all bcp and multiple methods for outbreak detection", y = "prevented fraction", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(comparison_df$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
ggsave("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/all_countries_frac_prevented_algo_comp.jpeg", scale = 1.2)


ggplot(comparison_df, aes(x = country, y = sens_per_outbreak, fill = source_method)) + 
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Comparison of all bcp and multiple methods for outbreak detection", y = "Sensitivity per outbreak", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(comparison_df$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
ggsave("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/all_countries_sens_per_outbreak_algo_comp.jpeg", scale = 1.2)


ggplot(comparison_df, aes(x = country, y = sens_exact, fill = source_method)) + 
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Comparison of all bcp and multiple methods for outbreak detection", y = "Exact sensitivity (detection within +/- 1 week)", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(comparison_df$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
ggsave("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/all_countries_sens_exact_algo_comp.jpeg", scale = 1.2)


ggplot(comparison_df, aes(x = country, y = PPV, fill = source_method)) + 
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Comparison of all bcp and multiple methods for outbreak detection", y = "Positive predictive value", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(comparison_df$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
ggsave("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/all_countries_PPV_algo_comp.jpeg", scale = 1.2)


ggplot(comparison_df, aes(x = country, y = specificity, fill = source_method)) + 
  geom_col(position = "dodge") + 
  coord_flip() + 
  labs(title = "Comparison of all bcp and multiple methods for outbreak detection", y = "Specificity", x = "", 
       caption = "WHO FluNet data were used as gold standard") +
  scale_x_discrete(limits = rev(levels(comparison_df$country))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set1")
ggsave("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/all_countries_specificity_algo_comp.jpeg", scale = 1.2)


mean_metrics <- comparison_df %>% group_by(source_method) %>% summarize_if(is.numeric, mean)
median_metrics <- comparison_df %>% group_by(source_method) %>% summarize_if(is.numeric, median)

mean_metrics_long <- pivot_longer(mean_metrics, cols = 2:7, names_to = "metric", values_to = "value")
median_metrics_long <- pivot_longer(median_metrics, cols = 2:7, names_to = "metric", values_to = "value")


ggplot(mean_metrics_long, aes(y = value, x = metric, fill = source_method)) + 
  geom_col(position = "dodge") + 
  labs(title = "Mean evaluation metrics", y ="", x = "", fill = "source and outbreak \ndetection method") + 
  scale_x_discrete(labels = c("prevented fraction", "positive \npredictive value", "exact sensitivity", "sensitivity \nper outbreak",
                              "sensitivity \nper week", "specificity")) + 
  scale_fill_brewer(palette = "Paired") + 
  scale_y_continuous(limits = c(0, 0.99))
ggsave("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/algo_comp_mean.jpeg")

ggplot(median_metrics_long, aes(y = value, x = metric, fill = source_method)) + 
  geom_col(position = "dodge") + 
  labs(title = "Median evaluation metrics", y ="", x = "", fill = "source and outbreak \ndetection method") + 
  scale_x_discrete(labels = c("prevented fraction", "positive \npredictive value", "exact sensitivity", "sensitivity \nper outbreak",
                              "sensitivity \nper week", "specificity")) + 
  scale_fill_brewer(palette = "Paired") + 
  scale_y_continuous(limits = c(0, 0.99))
ggsave("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/algo_comp_median.jpeg")

