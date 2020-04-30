library(ggplot2)
library(lubridate)
library(gamlss)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(readxl)


metrics <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics.csv", stringsAsFactors = TRUE)
metrics$false_alarm_rate <- 1 - metrics$specificity

metrics_long <- pivot_longer(metrics, cols = c(sens_per_outbreak, sens_per_week, sens_exact, PPV, specificity, frac_prevented), 
                             names_to = "metric", values_to = "values")

indicators <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/country_indicators.csv")


met_ind <- full_join(metrics, indicators, by = "country")
met_ind$HM_total_cat <- ordered(met_ind$HM_total_cat, levels = c("low", "medium","high"))
met_ind$EIOS_total_cat <- ordered(met_ind$EIOS_total_cat, levels = c("low", "medium","high"))

# select important predictors and log total and max counts
metrics_HM <- filter(met_ind, source == "EIOS") %>% 
  select(-c(X.x, X.y, FluNet_total_cat, FluNet_total, FluNet_max, EIOS_total_cat, EIOS_total, EIOS_max, ISO, 
            influenza_transmission_zone, ISO, problematic_FluNet, problematic_EIOS, problematic_HM, language)) %>% 
  mutate(latitude = abs(latitude), HM_total = log(HM_total), HM_max = log(HM_max))
# set counts for USA to NA so that they are disregarded in regressions, because they are outliers
# metrics_HM[metrics_HM$country == "United States", 11:12] <- NA
metrics_EIOS <- filter(met_ind, source == "EIOS") %>% 
  select(-c(X.x, X.y, FluNet_total_cat, FluNet_total, FluNet_max, HM_total_cat, HM_total, HM_max, ISO, 
            influenza_transmission_zone, ISO, problematic_FluNet, problematic_EIOS, problematic_HM, language)) %>% 
  mutate(latitude = abs(latitude), EIOS_total = log(EIOS_total), EIOS_max = log(EIOS_max))
# set counts for USA to NA so that they are disregarded in regressions, because they are outliers
# metrics_EIOS[metrics_EIOS$country == "United States", 11:12] <- NA

# list of all predictors
pred_EIOS <- names(metrics_EIOS)[10:19]
pred_HM <- names(metrics_HM)[10:19]

##### univariable regressions#####
### EIOS
coef <- c(1, 2, 5, 7, 9, 10, seq(from = 13, to = 23, by = 2))
pval <- c(3, 4, 6, 8, 11, 12, seq(from = 14, to = 24, by = 2))
# sensitivity per outbreak
lm_EIOS_sensOB <- lapply(10:19, function(x) lm(metrics_EIOS$sens_per_outbreak ~ metrics_EIOS[,x]))
names(lm_EIOS_sensOB) <- names(metrics_EIOS)[10:19]
lm_EIOS_sensOB_summary <- lapply(lm_EIOS_sensOB, summary)
lm_EIOS_sensOB_coef_list <- lapply(lm_EIOS_sensOB_summary, function(x) x$coefficients[-1, c(1,4)])
lm_EIOS_sensOB_coef_vec <- unlist(lm_EIOS_sensOB_coef_list)
lm_EIOS_sensOB_coef <- data.frame("coefficient" = lm_EIOS_sensOB_coef_vec[coef], "p-value" = lm_EIOS_sensOB_coef_vec[pval])
row.names(lm_EIOS_sensOB_coef) <- c("EIOS_total_cat.L", "EIOS_total_cat.Q", "EIOS_total", "EIOS_max", "global_region.temp_South", 
                                  "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_EIOS_sensOB_coef, 3)

lm_EIOS_sensOB_R2 <- unlist(lapply(lm_EIOS_sensOB_summary, function(x) x$r.squared))
lm_EIOS_sensOB_R2

par(mfrow = c(2, 2))
sapply(lm_EIOS_sensOB, plot)

# exact sensitivity 
lm_EIOS_sens_exact <- lapply(10:19, function(x) lm(metrics_EIOS$sens_exact ~ metrics_EIOS[,x]))
lm_EIOS_sens_exact_summary <- lapply(lm_EIOS_sens_exact, summary)
lm_EIOS_sens_exact_coef_list <- lapply(lm_EIOS_sens_exact_summary, function(x) x$coefficients[-1, c(1,4)])
lm_EIOS_sens_exact_coef_vec <- unlist(lm_EIOS_sens_exact_coef_list)
lm_EIOS_sens_exact_coef <- data.frame("coefficient" = lm_EIOS_sens_exact_coef_vec[coef], "p-value" = lm_EIOS_sens_exact_coef_vec[pval])
row.names(lm_EIOS_sens_exact_coef) <- c("EIOS_total_cat.L", "EIOS_total_cat.Q", "EIOS_total", "EIOS_max", "global_region.temp_South", 
                                  "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_EIOS_sens_exact_coef, 3)

lm_EIOS_sens_exact_R2 <- unlist(lapply(lm_EIOS_sens_exact_summary, function(x) x$r.squared))
names(lm_EIOS_sens_exact_R2) <- pred_EIOS
lm_EIOS_sens_exact_R2

par(mfrow = c(2, 2))
sapply(lm_EIOS_sens_exact, plot)  


zinf_EIOS_sens_exact <- lapply(10:19, function(x) gamlss(metrics_EIOS$sens_exact ~ metrics_EIOS[,x], family = "BEZI",
                                                         trace = F))
zinf_EIOS_sens_exact_summary <- lapply(zinf_EIOS_sens_exact, summary)
zinf_EIOS_sens_exact_coef <- sapply(zinf_EIOS_sens_exact_summary, '[', ,1)
zinf_EIOS_sens_exact_coef <- unlist(zinf_EIOS_sens_exact_coef)
zinf_EIOS_sens_exact_coef <- zinf_EIOS_sens_exact_coef[!grepl(x = names(zinf_EIOS_sens_exact_coef), pattern = "Intercept")]
zinf_EIOS_sens_exact_pval <- sapply(zinf_EIOS_sens_exact_summary, '[', ,4)
zinf_EIOS_sens_exact_pval <- unlist(zinf_EIOS_sens_exact_pval)
zinf_EIOS_sens_exact_pval <- zinf_EIOS_sens_exact_pval[!grepl(x = names(zinf_EIOS_sens_exact_pval), pattern = "Intercept")]

zinf_EIOS_sens_exact_df <- data.frame("coefficient" = zinf_EIOS_sens_exact_coef, "p-value" = zinf_EIOS_sens_exact_pval)
row.names(zinf_EIOS_sens_exact_df) <- c("EIOS_total_cat.L", "EIOS_total_cat.Q", "EIOS_total", "EIOS_max", "global_region.temp_South", 
                                        "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")


# sensitivity per week
lm_EIOS_sens_week <- lapply(10:19, function(x) lm(metrics_EIOS$sens_per_week ~ metrics_EIOS[,x]))
lm_EIOS_sens_week_summary <- lapply(lm_EIOS_sens_week, summary)
lm_EIOS_sens_week_coef_list <- lapply(lm_EIOS_sens_week_summary, function(x) x$coefficients[-1, c(1,4)])
lm_EIOS_sens_week_coef_vec <- unlist(lm_EIOS_sens_week_coef_list)
lm_EIOS_sens_week_coef <- data.frame("coefficient" = lm_EIOS_sens_week_coef_vec[coef], "p-value" = lm_EIOS_sens_week_coef_vec[pval])
row.names(lm_EIOS_sens_week_coef) <- c("EIOS_total_cat.L", "EIOS_total_cat.Q", "EIOS_total", "EIOS_max", "global_region.temp_South", 
                                      "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_EIOS_sens_week_coef, 3)

lm_EIOS_sens_week_R2 <- unlist(lapply(lm_EIOS_sens_week_summary, function(x) x$r.squared))
names(lm_EIOS_sens_week_R2) <- pred_EIOS
lm_EIOS_sens_week_R2

par(mfrow = c(2, 2))
sapply(lm_EIOS_sens_week, plot) 


# specificity 
lm_EIOS_spec <- lapply(10:19, function(x) lm(metrics_EIOS$specificity ~ metrics_EIOS[,x]))
lm_EIOS_spec_summary <- lapply(lm_EIOS_spec, summary)
lm_EIOS_spec_coef_list <- lapply(lm_EIOS_spec_summary, function(x) x$coefficients[-1, c(1,4)])
lm_EIOS_spec_coef_vec <- unlist(lm_EIOS_spec_coef_list)
lm_EIOS_spec_coef <- data.frame("coefficient" = lm_EIOS_spec_coef_vec[coef], "p-value" = lm_EIOS_spec_coef_vec[pval])
row.names(lm_EIOS_spec_coef) <- c("EIOS_total_cat.L", "EIOS_total_cat.Q", "EIOS_total", "EIOS_max", "global_region.temp_South", 
                                     "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_EIOS_spec_coef, 5)

lm_EIOS_spec_R2 <- unlist(lapply(lm_EIOS_spec_summary, function(x) x$r.squared))
names(lm_EIOS_spec_R2) <- pred_EIOS
lm_EIOS_spec_R2

par(mfrow = c(2, 2))
sapply(lm_EIOS_spec, plot) 


# timeliness
lm_EIOS_frac_prev <- lapply(10:19, function(x) lm(metrics_EIOS$frac_prevented ~ metrics_EIOS[,x]))
lm_EIOS_frac_prev_summary <- lapply(lm_EIOS_frac_prev, summary)
lm_EIOS_frac_prev_coef_list <- lapply(lm_EIOS_frac_prev_summary, function(x) x$coefficients[-1, c(1,4)])
lm_EIOS_frac_prev_coef_vec <- unlist(lm_EIOS_frac_prev_coef_list)
lm_EIOS_frac_prev_coef <- data.frame("coefficient" = lm_EIOS_frac_prev_coef_vec[coef], "p-value" = lm_EIOS_frac_prev_coef_vec[pval])
row.names(lm_EIOS_frac_prev_coef) <- c("EIOS_total_cat.L", "EIOS_total_cat.Q", "EIOS_total", "EIOS_max", "global_region.temp_South", 
                                "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_EIOS_frac_prev_coef, 3)

lm_EIOS_frac_prev_R2 <- unlist(lapply(lm_EIOS_frac_prev_summary, function(x) x$r.squared))
names(lm_EIOS_frac_prev_R2) <- pred_EIOS
lm_EIOS_frac_prev_R2

par(mfrow = c(2, 2))
sapply(lm_EIOS_frac_prev, plot) 


# PPV
lm_EIOS_PPV <- lapply(10:19, function(x) lm(metrics_EIOS$PPV ~ metrics_EIOS[,x]))
lm_EIOS_PPV_summary <- lapply(lm_EIOS_PPV, summary)
lm_EIOS_PPV_coef_list <- lapply(lm_EIOS_PPV_summary, function(x) x$coefficients[-1, c(1,4)])
lm_EIOS_PPV_coef_vec <- unlist(lm_EIOS_PPV_coef_list)
lm_EIOS_PPV_coef <- data.frame("coefficient" = lm_EIOS_PPV_coef_vec[coef], "p-value" = lm_EIOS_PPV_coef_vec[pval])
row.names(lm_EIOS_PPV_coef) <- c("EIOS_total_cat.L", "EIOS_total_cat.Q", "EIOS_total", "EIOS_max", "global_region.temp_South", 
                                "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_EIOS_PPV_coef, 3)

lm_EIOS_PPV_R2 <- unlist(lapply(lm_EIOS_PPV_summary, function(x) x$r.squared))
names(lm_EIOS_PPV_R2) <- pred_EIOS
lm_EIOS_PPV_R2

par(mfrow = c(2, 2))
sapply(lm_EIOS_PPV, plot) 



### HealthMap
# sensitivity per outbreak
lm_HM_sensOB <- lapply(10:19, function(x) lm(metrics_HM$sens_per_outbreak ~ metrics_HM[,x]))
lm_HM_sensOB_summary <- lapply(lm_HM_sensOB, summary)
lm_HM_sensOB_coef_list <- lapply(lm_HM_sensOB_summary, function(x) x$coefficients[-1, c(1,4)])
lm_HM_sensOB_coef_vec <- unlist(lm_HM_sensOB_coef_list)
lm_HM_sensOB_coef <- data.frame("coefficient" = lm_HM_sensOB_coef_vec[coef], "p-value" = lm_HM_sensOB_coef_vec[pval])
row.names(lm_HM_sensOB_coef) <- c("HM_total_cat.L", "HM_total_cat.Q", "HM_total", "HM_max", "global_region.temp_South", 
                                  "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_HM_sensOB_coef, 3)

lm_HM_sensOB_R2 <- unlist(lapply(lm_HM_sensOB_summary, function(x) x$r.squared))
names(lm_HM_sensOB_R2) <- pred_HM
lm_HM_sensOB_R2

par(mfrow = c(2, 2))
sapply(lm_HM_sensOB, plot)

# exact sensitivity 
lm_HM_sens_exact <- lapply(10:19, function(x) lm(metrics_HM$sens_exact ~ metrics_HM[,x]))
lm_HM_sens_exact_summary <- lapply(lm_HM_sens_exact, summary)
lm_HM_sens_exact_coef_list <- lapply(lm_HM_sens_exact_summary, function(x) x$coefficients[-1, c(1,4)])
lm_HM_sens_exact_coef_vec <- unlist(lm_HM_sens_exact_coef_list)
lm_HM_sens_exact_coef <- data.frame("coefficient" = lm_HM_sens_exact_coef_vec[coef], "p-value" = lm_HM_sens_exact_coef_vec[pval])
row.names(lm_HM_sens_exact_coef) <- c("HM_total_cat.L", "HM_total_cat.Q", "HM_total", "HM_max", "global_region.temp_South", 
                                      "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_HM_sens_exact_coef, 3)

lm_HM_sens_exact_R2 <- unlist(lapply(lm_HM_sens_exact_summary, function(x) x$r.squared))
names(lm_HM_sens_exact_R2) <- pred_HM
lm_HM_sens_exact_R2

par(mfrow = c(2, 2))
sapply(lm_HM_sens_exact, plot)  


# sensitivity per week
lm_HM_sens_week <- lapply(10:19, function(x) lm(metrics_HM$sens_per_week ~ metrics_HM[,x]))
lm_HM_sens_week_summary <- lapply(lm_HM_sens_week, summary)
lm_HM_sens_week_coef_list <- lapply(lm_HM_sens_week_summary, function(x) x$coefficients[-1, c(1,4)])
lm_HM_sens_week_coef_vec <- unlist(lm_HM_sens_week_coef_list)
lm_HM_sens_week_coef <- data.frame("coefficient" = lm_HM_sens_week_coef_vec[coef], "p-value" = lm_HM_sens_week_coef_vec[pval])
row.names(lm_HM_sens_week_coef) <- c("HM_total_cat.L", "HM_total_cat.Q", "HM_total", "HM_max", "global_region.temp_South", 
                                     "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_HM_sens_week_coef, 3)

lm_HM_sens_week_R2 <- unlist(lapply(lm_HM_sens_week_summary, function(x) x$r.squared))
names(lm_HM_sens_week_R2) <- pred_HM
lm_HM_sens_week_R2

par(mfrow = c(2, 2))
sapply(lm_HM_sens_week, plot) 


# specificity 
lm_HM_spec <- lapply(10:19, function(x) lm(log(metrics_HM$specificity) ~ metrics_HM[,x]))
lm_HM_spec_summary <- lapply(lm_HM_spec, summary)
lm_HM_spec_coef_list <- lapply(lm_HM_spec_summary, function(x) x$coefficients[-1, c(1,4)])
lm_HM_spec_coef_vec <- unlist(lm_HM_spec_coef_list)
lm_HM_spec_coef <- data.frame("coefficient" = lm_HM_spec_coef_vec[coef], "p-value" = lm_HM_spec_coef_vec[pval])
row.names(lm_HM_spec_coef) <- c("HM_total_cat.L", "HM_total_cat.Q", "HM_total", "HM_max", "global_region.temp_South", 
                                "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_HM_spec_coef, 3)

lm_HM_spec_R2 <- unlist(lapply(lm_HM_spec_summary, function(x) x$r.squared))
names(lm_HM_spec_R2) <- pred_HM
lm_HM_spec_R2

par(mfrow = c(2, 2))
sapply(lm_HM_spec, plot) 


# timeliness
lm_HM_frac_prev <- lapply(10:19, function(x) lm(metrics_HM$frac_prevented ~ metrics_HM[,x]))
lm_HM_frac_prev_summary <- lapply(lm_HM_frac_prev, summary)
lm_HM_frac_prev_coef_list <- lapply(lm_HM_frac_prev_summary, function(x) x$coefficients[-1, c(1,4)])
lm_HM_frac_prev_coef_vec <- unlist(lm_HM_frac_prev_coef_list)
lm_HM_frac_prev_coef <- data.frame("coefficient" = lm_HM_frac_prev_coef_vec[coef], "p-value" = lm_HM_frac_prev_coef_vec[pval])
row.names(lm_HM_frac_prev_coef) <- c("HM_total_cat.L", "HM_total_cat.Q", "HM_total", "HM_max", "global_region.temp_South", 
                                     "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_HM_frac_prev_coef, 3)

lm_HM_frac_prev_R2 <- unlist(lapply(lm_HM_frac_prev_summary, function(x) x$r.squared))
names(lm_HM_frac_prev_R2) <- pred_HM
lm_HM_frac_prev_R2

par(mfrow = c(2, 2))
sapply(lm_HM_frac_prev, plot) 


# PPV
lm_HM_PPV <- lapply(10:19, function(x) lm(metrics_HM$PPV ~ metrics_HM[,x]))
lm_HM_PPV_summary <- lapply(lm_HM_PPV, summary)
lm_HM_PPV_coef_list <- lapply(lm_HM_PPV_summary, function(x) x$coefficients[-1, c(1,4)])
lm_HM_PPV_coef_vec <- unlist(lm_HM_PPV_coef_list)
lm_HM_PPV_coef <- data.frame("coefficient" = lm_HM_PPV_coef_vec[coef], "p-value" = lm_HM_PPV_coef_vec[pval])
row.names(lm_HM_PPV_coef) <- c("HM_total_cat.L", "HM_total_cat.Q", "HM_total", "HM_max", "global_region.temp_South", 
                               "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_HM_PPV_coef, 3)

lm_HM_PPV_R2 <- unlist(lapply(lm_HM_PPV_summary, function(x) x$r.squared))
names(lm_HM_PPV_R2) <- pred_HM
lm_HM_PPV_R2

par(mfrow = c(2, 2))
sapply(lm_HM_PPV, plot) 


#### regressions with interactions #####

### HealthMap
predictors <- colnames(metrics_HM)[10:19]
combos <- combn(predictors, 2, simplify = FALSE)
combos <- sapply(combos, FUN = paste, collapse = "*")
# sensitivity per outbreak
lm_HM_sensOB_int <- lapply(combos, FUN = function(x) {
  frm <- as.formula(paste("sens_per_outbreak ~ ", x))
  summary(lm(frm, data = metrics_HM))
})
names(lm_HM_sensOB_int) <- combos

lm_HM_sensOB_int_coef <- lapply(lm_HM_sensOB_int, '[[', "coefficients")
lm_HM_sensOB_int_pval <- sapply(lm_HM_sensOB_int_coef, '[', ,4)
lm_HM_sensOB_int_pval <- unlist(lapply(lm_HM_sensOB_int_pval, last))
names(lm_HM_sensOB_int_pval)[lm_HM_sensOB_int_pval < 0.05]
# no significant interactions 

# sensitivity per week
lm_HM_sens_week_int <- lapply(combos, FUN = function(x) {
  frm <- as.formula(paste("sens_per_week ~ ", x))
  summary(lm(frm, data = metrics_HM))
})
names(lm_HM_sens_week_int) <- combos

lm_HM_sens_week_int_coef <- lapply(lm_HM_sens_week_int, '[[', "coefficients")
lm_HM_sens_week_int_pval <- sapply(lm_HM_sens_week_int_coef, '[', ,4)
lm_HM_sens_week_int_pval <- unlist(lapply(lm_HM_sens_week_int_pval, last))
names(lm_HM_sens_week_int_pval)[lm_HM_sens_week_int_pval < 0.05]


# specificity
lm_HM_spec_int <- lapply(combos, FUN = function(x) {
  frm <- as.formula(paste("specificity ~ ", x))
  summary(lm(frm, data = metrics_HM))
})
names(lm_HM_spec_int) <- combos

lm_HM_spec_int_coef <- lapply(lm_HM_spec_int, '[[', "coefficients")
lm_HM_spec_int_pval <- sapply(lm_HM_spec_int_coef, '[', ,4)
lm_HM_spec_int_pval <- unlist(lapply(lm_HM_spec_int_pval, last))
names(lm_HM_spec_int_pval)[lm_HM_spec_int_pval < 0.05]


# timeliness
lm_HM_frac_prevented_int <- lapply(combos, FUN = function(x) {
  frm <- as.formula(paste("frac_prevented ~ ", x))
  summary(lm(frm, data = metrics_HM))
})
names(lm_HM_frac_prevented_int) <- combos

lm_HM_frac_prevented_int_coef <- lapply(lm_HM_frac_prevented_int, '[[', "coefficients")
lm_HM_frac_prevented_int_pval <- sapply(lm_HM_frac_prevented_int_coef, '[', ,4)
lm_HM_frac_prevented_int_pval <- unlist(lapply(lm_HM_frac_prevented_int_pval, last))
names(lm_HM_frac_prevented_int_pval)[lm_HM_frac_prevented_int_pval < 0.05]


# PPV
lm_HM_PPV_int <- lapply(combos, FUN = function(x) {
  frm <- as.formula(paste("PPV ~ ", x))
  summary(lm(frm, data = metrics_HM))
})
names(lm_HM_PPV_int) <- combos

lm_HM_PPV_int_coef <- lapply(lm_HM_PPV_int, '[[', "coefficients")
lm_HM_PPV_int_pval <- sapply(lm_HM_PPV_int_coef, '[', ,4)
lm_HM_PPV_int_pval <- unlist(lapply(lm_HM_PPV_int_pval, last))
names(lm_HM_PPV_int_pval)[lm_HM_PPV_int_pval < 0.05]



### EIOS
predictors <- colnames(metrics_EIOS)[10:19]
combos <- combn(predictors, 2, simplify = FALSE)
combos <- sapply(combos, FUN = paste, collapse = "*")
# sensitivity per outbreak
lm_EIOS_sensOB_int <- lapply(combos, FUN = function(x) {
  frm <- as.formula(paste("sens_per_outbreak ~ ", x))
  summary(lm(frm, data = metrics_EIOS))
})
names(lm_EIOS_sensOB_int) <- combos

lm_EIOS_sensOB_int_coef <- lapply(lm_EIOS_sensOB_int, '[[', "coefficients")
lm_EIOS_sensOB_int_pval <- sapply(lm_EIOS_sensOB_int_coef, '[', ,4)
lm_EIOS_sensOB_int_pval <- unlist(lapply(lm_EIOS_sensOB_int_pval, last))
names(lm_EIOS_sensOB_int_pval)[lm_EIOS_sensOB_int_pval < 0.05]


# sensitivity per week
lm_EIOS_sens_week_int <- lapply(combos, FUN = function(x) {
  frm <- as.formula(paste("sens_per_week ~ ", x))
  summary(lm(frm, data = metrics_EIOS))
})
names(lm_EIOS_sens_week_int) <- combos

lm_EIOS_sens_week_int_coef <- lapply(lm_EIOS_sens_week_int, '[[', "coefficients")
lm_EIOS_sens_week_int_pval <- sapply(lm_EIOS_sens_week_int_coef, '[', ,4)
lm_EIOS_sens_week_int_pval <- unlist(lapply(lm_EIOS_sens_week_int_pval, last))
names(lm_EIOS_sens_week_int_pval)[lm_EIOS_sens_week_int_pval < 0.05]


# specificity
lm_EIOS_spec_int <- lapply(combos, FUN = function(x) {
  frm <- as.formula(paste("specificity ~ ", x))
  summary(lm(frm, data = metrics_EIOS))
})
names(lm_EIOS_spec_int) <- combos

lm_EIOS_spec_int_coef <- lapply(lm_EIOS_spec_int, '[[', "coefficients")
lm_EIOS_spec_int_pval <- sapply(lm_EIOS_spec_int_coef, '[', ,4)
lm_EIOS_spec_int_pval <- unlist(lapply(lm_EIOS_spec_int_pval, last))
names(lm_EIOS_spec_int_pval)[lm_EIOS_spec_int_pval < 0.05]


# timeliness
lm_EIOS_frac_prevented_int <- lapply(combos, FUN = function(x) {
  frm <- as.formula(paste("frac_prevented ~ ", x))
  summary(lm(frm, data = metrics_EIOS))
})
names(lm_EIOS_frac_prevented_int) <- combos

lm_EIOS_frac_prevented_int_coef <- lapply(lm_EIOS_frac_prevented_int, '[[', "coefficients")
lm_EIOS_frac_prevented_int_pval <- sapply(lm_EIOS_frac_prevented_int_coef, '[', ,4)
lm_EIOS_frac_prevented_int_pval <- unlist(lapply(lm_EIOS_frac_prevented_int_pval, last))
names(lm_EIOS_frac_prevented_int_pval)[lm_EIOS_frac_prevented_int_pval < 0.05]


# PPV
lm_EIOS_PPV_int <- lapply(combos, FUN = function(x) {
  frm <- as.formula(paste("PPV ~ ", x))
  summary(lm(frm, data = metrics_EIOS))
})
names(lm_EIOS_PPV_int) <- combos

lm_EIOS_PPV_int_coef <- lapply(lm_EIOS_PPV_int, '[[', "coefficients")
lm_EIOS_PPV_int_pval <- sapply(lm_EIOS_PPV_int_coef, '[', ,4)
lm_EIOS_PPV_int_pval <- unlist(lapply(lm_EIOS_PPV_int_pval, last))
names(lm_EIOS_PPV_int_pval)[lm_EIOS_PPV_int_pval < 0.05]




##### multivariable regressions #####

### HealthMap
lm_HM_multi <- lapply(3:8, function(x) lm(metrics_HM[,x] ~ metrics_HM$HM_total_cat + metrics_HM$global_region + metrics_HM$english + 
                                          metrics_HM$HDI.2018 + metrics_HM$latitude))
lm_HM_multi_summary <- lapply(lm_HM_multi, summary)
names(lm_HM_multi_summary) <- names(metrics_HM)[3:8]
lm_HM_multi_coef <- lapply(lm_HM_multi_summary, '[[', "coefficients")

# step function based on AIC to select best model
lm_HM_multi_step <- lapply(lapply(lm_HM_multi, step), summary)
lm_HM_multi_step_coef <- lapply(lm_HM_multi_step, '[[', "coefficients")
names(lm_HM_multi_step_coef) <- names(metrics_HM)[3:8]
lm_HM_multi_step_coef

lm_HM_multi_step_R2 <- sapply(lm_HM_multi_step, '[[', "r.squared")
names(lm_HM_multi_step_R2) <- names(metrics_HM)[3:8]
lm_HM_multi_step_R2

lm_HM_multi_step_adjR2 <- sapply(lm_HM_multi_step, '[[', "adj.r.squared")
names(lm_HM_multi_step_adjR2) <- names(metrics_HM)[3:8]
lm_HM_multi_step_adjR2


### EIOS
lm_EIOS_multi <- lapply(3:8, function(x) lm(metrics_EIOS[,x] ~ metrics_EIOS$EIOS_total_cat + metrics_EIOS$global_region + metrics_EIOS$english + 
                                            metrics_EIOS$HDI.2018 + metrics_EIOS$latitude))
lm_EIOS_multi_summary <- lapply(lm_EIOS_multi, summary)
names(lm_EIOS_multi_summary) <- names(metrics_EIOS)[3:8]
lm_EIOS_multi_coef <- lapply(lm_EIOS_multi_summary, '[[', "coefficients")

# step function based on AIC to select best model
lm_EIOS_multi_step <- lapply(lapply(lm_EIOS_multi, step), summary)
lm_EIOS_multi_step_coef <- lapply(lm_EIOS_multi_step, '[[', "coefficients")
names(lm_EIOS_multi_step_coef) <- names(metrics_EIOS)[3:8]
lm_EIOS_multi_step_coef

lm_EIOS_multi_step_R2 <- sapply(lm_EIOS_multi_step, '[[', "r.squared")
names(lm_EIOS_multi_step_R2) <- names(metrics_EIOS)[3:8]
lm_EIOS_multi_step_R2

lm_EIOS_multi_step_adjR2 <- sapply(lm_EIOS_multi_step, '[[', "adj.r.squared")
names(lm_EIOS_multi_step_adjR2) <- names(metrics_EIOS)[3:8]
lm_EIOS_multi_step_adjR2


#### advanced variable selection ####
library(caret)
library(leaps)
library(MASS)
library(relaimpo)
library(car)
library(gvlma)

full_model_HM_sensOB <- lm(sens_per_outbreak ~ HM_total_cat + HM_total + HM_max + global_region + english + 
                             HDI.2018 + latitude + longitude + PFI.2018 + TIU.2017 + HM_filter_lang, 
                           data = metrics_HM)
vif(full_model_HM_sensOB)
reduced_model_HM_sensOB <- lm(sens_per_outbreak ~ HM_total + global_region + english + 
                             HDI.2018 + latitude + PFI.2018 + HM_filter_lang, 
                           data = metrics_HM)
vif(reduced_model_HM_sensOB) # all these predictors don't have collinearity (VIF < 3)
summary(reduced_model_HM_sensOB)

lm_HM_multi <- lapply(3:8, function(x) lm(metrics_HM[,x] ~ metrics_HM$HM_total + metrics_HM$global_region + metrics_HM$english + 
                                            metrics_HM$HDI.2018 + metrics_HM$latitude + metrics_HM$PFI.2018 + metrics_HM$HM_filter_lang))
AIC_HM_multi <- lapply(lm_HM_multi, stepAIC, direction = "both", trace = FALSE)
AIC_HM_multi_summary <- lapply(AIC_HM_multi, summary)
AIC_HM_multi_coef <- lapply(AIC_HM_multi_summary, '[[', "coefficients")
names(AIC_HM_multi_coef) <- names(metrics_HM)[3:8]

AIC_HM_multi_AIC <- sapply(AIC_HM_multi, AIC)
names(AIC_HM_multi_AIC) <- names(metrics_HM)[3:8]

AIC_HM_multi_radj <- sapply(AIC_HM_multi_summary, '[[', "adj.r.squared")
names(AIC_HM_multi_radj) <- names(metrics_HM)[3:8]


AIC_HM_multi_outlier <- lapply(AIC_HM_multi, outlierTest)
AIC_HM_multi_leverageplot <- lapply(AIC_HM_multi, leveragePlots)

AIC_HM_multi_spreadplot <- lapply(AIC_HM_multi, spreadLevelPlot)
AIC_HM_multi_ncvtest <- lapply(AIC_HM_multi, ncvTest)
AIC_HM_multi_ncvtest_p <- sapply(AIC_HM_multi_ncvtest, '[[', "p")

AIC_HM_multi_res <- data.frame(sapply(AIC_HM_multi, '[[', "residuals"))
colnames(AIC_HM_multi_res) <- names(metrics_HM)[3:8]

par(mfrow = c(1, 1))
AIC_HM_multi_hist_res <- apply(AIC_HM_multi_res, 2, hist, main = "")
AIC_HM_multi_qqplot <- lapply(AIC_HM_multi, qqPlot)

par(mfrow= c(2, 2))
lapply(AIC_HM_multi, plot)

AIC_HM_multi_crplot <- lapply(AIC_HM_multi, crPlots)

AIC_HM_multi_gvlma <- lapply(AIC_HM_multi, gvlma)



leaps_HM_sensOB <- regsubsets(sens_per_outbreak ~ HM_total + global_region + english + 
                               HDI.2018 + latitude + PFI.2018 + HM_filter_lang, 
                             data = metrics_HM, nvmax = 5, method = "seqrep", all.best = T)
summary(leaps_HM_sensOB)



set.seed(123)
train.control <- trainControl(method = "cv", number = 5)
train_HM_sensOB <- train(sens_per_outbreak ~ HM_total + global_region + english + 
                           HDI.2018 + latitude + PFI.2018 + HM_filter_lang, 
                         data = metrics_HM,
                        method = "leapSeq", 
                        tuneGrid = data.frame(nvmax = 1:5),
                        metric = "Rsquared",
                        trControl = train.control
)
train_HM_sensOB$results
train_HM_sensOB$bestTune
summary(train_HM_sensOB$finalModel)
# cross-validation gives weird results, r-squared is higher with 1 variable than with multiple


train_HM_sensOB <- train(sens_per_outbreak ~ HM_total + global_region + english + 
                           HDI.2018 + latitude + PFI.2018 + HM_filter_lang, 
                         data = metrics_HM,
                         method = "lmStepAIC", 
                         trace = FALSE,
                         trControl = train.control,
                         metric = "Rsquared"
)
train_HM_sensOB$results
train_HM_sensOB$bestTune
summary(train_HM_sensOB$finalModel)


calc.relimp(reduced_model_HM_sensOB, type = c("lmg","last","first"),
            rela=TRUE)


outlierTest(full_model_HM_sensOB)
qqPlot(full_model_HM_sensOB, main="QQ Plot")
leveragePlots(full_model_HM_sensOB)
influencePlot(full_model_HM_sensOB, id.method="identify")
ncvTest(full_model_HM_sensOB)
spreadLevelPlot(full_model_HM_sensOB)

# variance inflation factors
vif(full_model_HM_sensOB)
sqrt(vif((full_model_HM_sensOB)))

crPlots(full_model_HM_sensOB)

