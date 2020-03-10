library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)



### load metrics and indicator data
metrics <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics.csv", stringsAsFactors = TRUE)
metrics$false_alarm_rate <- 1 - metrics$specificity

metrics_long <- pivot_longer(metrics, cols = c(sens_per_outbreak, sens_per_week, sens_exact, PPV, specificity, frac_prevented), 
                             names_to = "metric", values_to = "values")

indicators <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/country_indicators.csv")


met_ind <- full_join(metrics, indicators, by = "country")
met_ind$HM_total_cat <- ordered(met_ind$HM_total_cat, levels = c("low", "medium","high"))
met_ind$EIOS_total_cat <- ordered(met_ind$EIOS_total_cat, levels = c("low", "medium","high"))

##### ROC & AMOC curves #####
### HealthMap
# outbreak sensitivity, false alarm rate
metrics_HM <- filter(met_ind, source == "HealthMap")
ggplot(metrics_HM, aes(x = false_alarm_rate, y = sens_per_outbreak)) + 
  geom_point(aes(shape = problematic_FluNet), size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = metrics_HM$country, col = HM_total_cat)) + 
  labs(title = "ROC curve HealthMap", y = "sensitivity per outbreak", x = "false alarm rate", 
       col = "Total event counts", shape = "FluNet data problems") +
  scale_y_continuous(limits = c(0, 1)) + 
  scale_color_brewer(palette = "Set1") 

# exact sensitivity, false alarm rate
ggplot(metrics_HM, aes(x = false_alarm_rate, y = sens_exact)) + 
  geom_point(aes(shape = problematic_FluNet), size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = metrics_HM$country, col = HM_total_cat)) + 
  labs(title = "ROC curve HealthMap", y = "sensitivity (detection within +/- 1 week)", x = "false alarm rate", 
       col = "Total event counts", shape = "FluNet data problems") +
  scale_y_continuous(limits = c(0, 0.5)) + 
  scale_color_brewer(palette = "Set1")

# outbreak timeliness, false alarm rate
ggplot(metrics_HM, aes(x = false_alarm_rate, y = frac_prevented)) + 
  geom_point(aes(shape = problematic_FluNet), size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = metrics_HM$country, col = HM_total_cat)) + 
  labs(title = "AMOC curve HealthMap", y = "timeliness (prevented fraction)", x = "false alarm rate", 
       col = "Total event counts", shape = "FluNet data problems") +
  scale_y_continuous(limits = c(0, 1)) + 
  scale_color_brewer(palette = "Set1") 

# outbreak sensitivity, PPV
ggplot(metrics_HM, aes(x = PPV, y = sens_per_outbreak)) + 
  geom_point(aes(shape = problematic_FluNet), size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = metrics_HM$country, col = HM_total_cat)) + 
  labs(title = "ROC curve HealthMap", y = "sensitivity per outbreak", x = "positive predictive value", 
       col = "Total event counts", shape = "FluNet data problems") +
  scale_y_continuous(limits = c(0, 1)) + 
  scale_color_brewer(palette = "Set1") 
  

### EIOS
metrics_EIOS <- filter(met_ind, source == "EIOS")
ggplot(metrics_EIOS, aes(x = false_alarm_rate, y = sens_per_outbreak)) + 
  geom_point(aes(shape = problematic_FluNet), size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = metrics_EIOS$country, col = EIOS_total_cat)) + 
  labs(title = "ROC curve EIOS", y = "sensitivity per outbreak", x = "false alarm rate", 
       col = "Total event counts", shape = "FluNet data problems") +
  scale_y_continuous(limits = c(0, 1)) + 
  scale_color_brewer(palette = "Set1") 

# exact sensitivity, false alarm rate
ggplot(metrics_EIOS, aes(x = false_alarm_rate, y = sens_exact)) + 
  geom_point(aes(shape = problematic_FluNet), size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = metrics_EIOS$country, col = EIOS_total_cat)) + 
  labs(title = "ROC curve EIOS", y = "sensitivity (detection within +/- 1 week)", x = "false alarm rate", 
       col = "Total event counts", shape = "FluNet data problems") +
  scale_y_continuous(limits = c(0, 0.5)) + 
  scale_color_brewer(palette = "Set1")

# outbreak timeliness, false alarm rate
ggplot(metrics_EIOS, aes(x = false_alarm_rate, y = frac_prevented)) + 
  geom_point(aes(shape = problematic_FluNet), size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = metrics_EIOS$country, col = EIOS_total_cat)) + 
  labs(title = "AMOC curve EIOS", y = "timeliness (prevented fraction)", x = "false alarm rate", 
       col = "Total event counts", shape = "FluNet data problems") +
  scale_y_continuous(limits = c(0, 1)) + 
  scale_color_brewer(palette = "Set1") 

# outbreak sensitivity, PPV
ggplot(metrics_EIOS, aes(x = PPV, y = sens_per_outbreak)) + 
  geom_point(aes(shape = problematic_FluNet), size = 3, alpha = 0.8) + 
  geom_text_repel(aes(label = metrics_EIOS$country, col = EIOS_total_cat)) + 
  labs(title = "ROC curve EIOS", y = "sensitivity per outbreak", x = "positive predictive value", 
       col = "Total event counts", shape = "FluNet data problems") +
  scale_y_continuous(limits = c(0, 1)) + 
  scale_color_brewer(palette = "Set1") 


##### individual predictor plots #####
### HealthMap
# y is prevented fraction, x is data abundance, language, region
p1 <- ggplot(metrics_HM, aes(x = HM_total_cat, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts", y = "timeliness (prevented fraction)") + 
  scale_y_continuous(limits = c(0, 1))

p2 <- ggplot(metrics_HM, aes(x = problematic_FluNet, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "FluNet data problems", y = "") + 
  scale_y_continuous(limits = c(0, 1))  

p3 <- ggplot(metrics_HM, aes(x = english, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "English", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p4 <- ggplot(metrics_HM, aes(x = global_region, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "Geography", y = "timeliness (prevented fraction)") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = c("temp. Northern \n hemisphere", "temp. Southern \n hemisphere", "tropical"))

p5 <- ggplot(metrics_HM, aes(x = HM_total, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts over 6 years", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p6 <- ggplot(metrics_HM, aes(x = HM_max, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "Maximum event counts per week", y = "") + 
  scale_y_continuous(limits = c(0, 1))

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "HealthMap timeliness predictors")
ggsave(plot = grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "HealthMap timeliness predictors"),
       filename = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HM_timeliness_predictors.jpeg", scale = 1.3)

# y is outbreak sensitivity, x is data abundance, language, region
p1 <- ggplot(metrics_HM, aes(x = HM_total_cat, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts", y = "outbreak sensitivity") + 
  scale_y_continuous(limits = c(0, 1))

p2 <- ggplot(metrics_HM, aes(x = problematic_FluNet, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "FluNet data problems", y = "") + 
  scale_y_continuous(limits = c(0, 1))  

p3 <- ggplot(metrics_HM, aes(x = english, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "English", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p4 <- ggplot(metrics_HM, aes(x = global_region, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "Geography", y = "outbreak sensitivity") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = c("temp. Northern \n hemisphere", "temp. Southern \n hemisphere", "tropical"))

p5 <- ggplot(metrics_HM, aes(x = HM_total, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts over 6 years", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p6 <- ggplot(metrics_HM, aes(x = HM_max, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "Maximum event counts per week", y = "") + 
  scale_y_continuous(limits = c(0, 1))

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "HealthMap outbreak sensitivity predictors")
ggsave(plot = grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "HealthMap outbreak sensitivity predictors"),
       filename = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HM_outbreak_sens_predictors.jpeg", scale = 1.3)

# y is specificity, x is data abundance, language, region
p1 <- ggplot(metrics_HM, aes(x = HM_total_cat, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts", y = "specificity") + 
  scale_y_continuous(limits = c(0, 1))

p2 <- ggplot(metrics_HM, aes(x = problematic_FluNet, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "FluNet data problems", y = "") + 
  scale_y_continuous(limits = c(0, 1))  

p3 <- ggplot(metrics_HM, aes(x = english, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "English", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p4 <- ggplot(metrics_HM, aes(x = global_region, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "Geography", y = "specificity") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = c("temp. Northern \n hemisphere", "temp. Southern \n hemisphere", "tropical"))

p5 <- ggplot(metrics_HM, aes(x = HM_total, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts over 6 years", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p6 <- ggplot(metrics_HM, aes(x = HM_max, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "Maximum event counts per week", y = "") + 
  scale_y_continuous(limits = c(0, 1))

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "HealthMap specificity predictors")
ggsave(plot = grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "HealthMap specificity predictors"),
       filename = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HM_specificity_predictors.jpeg", scale = 1.3)


# y is PPV, x is data abundance, language, region
p1 <- ggplot(metrics_HM, aes(x = HM_total_cat, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts", y = "positive predictive value") + 
  scale_y_continuous(limits = c(0, 1))

p2 <- ggplot(metrics_HM, aes(x = problematic_FluNet, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "FluNet data problems", y = "") + 
  scale_y_continuous(limits = c(0, 1))  

p3 <- ggplot(metrics_HM, aes(x = english, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "English", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p4 <- ggplot(metrics_HM, aes(x = global_region, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "Geography", y = "positive predictive value") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = c("temp. Northern \n hemisphere", "temp. Southern \n hemisphere", "tropical"))

p5 <- ggplot(metrics_HM, aes(x = HM_total, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts over 6 years", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p6 <- ggplot(metrics_HM, aes(x = HM_max, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "Maximum event counts per week", y = "") + 
  scale_y_continuous(limits = c(0, 1))

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "HealthMap positive predictive value predictors")
ggsave(plot = grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "HealthMap positive predictive value predictors"),
       filename = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/HM_PPV_predictors.jpeg", scale = 1.3)



### EIOS
p1 <- ggplot(metrics_EIOS, aes(x = EIOS_total_cat, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts", y = "timeliness (prevented fraction)") + 
  scale_y_continuous(limits = c(0, 1))

p2 <- ggplot(metrics_EIOS, aes(x = problematic_FluNet, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "FluNet data problems", y = "") + 
  scale_y_continuous(limits = c(0, 1))  

p3 <- ggplot(metrics_EIOS, aes(x = english, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "English", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p4 <- ggplot(metrics_EIOS, aes(x = global_region, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "Geography", y = "timeliness (prevented fraction)") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = c("temp. Northern \n hemisphere", "temp. Southern \n hemisphere", "tropical"))

p5 <- ggplot(metrics_EIOS, aes(x = EIOS_total, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts over 2 years", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p6 <- ggplot(metrics_EIOS, aes(x = EIOS_max, y = frac_prevented)) + 
  geom_point(size = 2) + 
  labs(x = "Maximum event counts per week", y = "") + 
  scale_y_continuous(limits = c(0, 1))

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "EIOS timeliness predictors")
ggsave(plot = grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "EIOS timeliness predictors"),
       filename = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_timeliness_predictors.jpeg", scale = 1.3)


# y is outbreak sensitivity, x is data abundance, language, region
p1 <- ggplot(metrics_EIOS, aes(x = EIOS_total_cat, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts", y = "outbreak sensitivity") + 
  scale_y_continuous(limits = c(0, 1))

p2 <- ggplot(metrics_EIOS, aes(x = problematic_FluNet, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "FluNet data problems", y = "") + 
  scale_y_continuous(limits = c(0, 1))  

p3 <- ggplot(metrics_EIOS, aes(x = english, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "English", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p4 <- ggplot(metrics_EIOS, aes(x = global_region, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "Geography", y = "outbreak sensitivity") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = c("temp. Northern \n hemisphere", "temp. Southern \n hemisphere", "tropical"))

p5 <- ggplot(metrics_EIOS, aes(x = EIOS_total, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts over 6 years", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p6 <- ggplot(metrics_EIOS, aes(x = EIOS_max, y = sens_per_outbreak)) + 
  geom_point(size = 2) + 
  labs(x = "Maximum event counts per week", y = "") + 
  scale_y_continuous(limits = c(0, 1))

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "EIOS outbreak sensitivity predictors")
ggsave(plot = grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "EIOS outbreak sensitivity predictors"),
       filename = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_outbreak_sens_predictors.jpeg", scale = 1.3)


# y is specificity, x is data abundance, language, region
p1 <- ggplot(metrics_EIOS, aes(x = EIOS_total_cat, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts", y = "specificity") + 
  scale_y_continuous(limits = c(0, 1))

p2 <- ggplot(metrics_EIOS, aes(x = problematic_FluNet, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "FluNet data problems", y = "") + 
  scale_y_continuous(limits = c(0, 1))  

p3 <- ggplot(metrics_EIOS, aes(x = english, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "English", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p4 <- ggplot(metrics_EIOS, aes(x = global_region, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "Geography", y = "specificity") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = c("temp. Northern \n hemisphere", "temp. Southern \n hemisphere", "tropical"))

p5 <- ggplot(metrics_EIOS, aes(x = EIOS_total, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts over 6 years", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p6 <- ggplot(metrics_EIOS, aes(x = EIOS_max, y = specificity)) + 
  geom_point(size = 2) + 
  labs(x = "Maximum event counts per week", y = "") + 
  scale_y_continuous(limits = c(0, 1))

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "EIOS specificity predictors")
ggsave(plot = grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "EIOS specificity predictors"),
       filename = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_specificity_predictors.jpeg", scale = 1.3)

# y is PPV, x is data abundance, language, region
p1 <- ggplot(metrics_EIOS, aes(x = EIOS_total_cat, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts", y = "positive predictive value") + 
  scale_y_continuous(limits = c(0, 1))

p2 <- ggplot(metrics_EIOS, aes(x = problematic_FluNet, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "FluNet data problems", y = "") + 
  scale_y_continuous(limits = c(0, 1))  

p3 <- ggplot(metrics_EIOS, aes(x = english, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "English", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p4 <- ggplot(metrics_EIOS, aes(x = global_region, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "Geography", y = "positive predictive value") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = c("temp. Northern \n hemisphere", "temp. Southern \n hemisphere", "tropical"))

p5 <- ggplot(metrics_EIOS, aes(x = EIOS_total, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "Total event counts over 6 years", y = "") + 
  scale_y_continuous(limits = c(0, 1)) 

p6 <- ggplot(metrics_EIOS, aes(x = EIOS_max, y = PPV)) + 
  geom_point(size = 2) + 
  labs(x = "Maximum event counts per week", y = "") + 
  scale_y_continuous(limits = c(0, 1))

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "EIOS positive predictive value predictors")
ggsave(plot = grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = "EIOS positive predictive value predictors"),
       filename = "D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/EIOS_PPV_predictors.jpeg", scale = 1.3)

