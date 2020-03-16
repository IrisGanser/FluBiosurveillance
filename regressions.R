library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)
library(GGally)
library(readxl)

metrics <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics_all_bcp.csv", stringsAsFactors = TRUE)
metrics$false_alarm_rate <- 1 - metrics$specificity

metrics_long <- pivot_longer(metrics, cols = c(sens_per_outbreak, sens_per_week, sens_exact, PPV, specificity, frac_prevented), 
                             names_to = "metric", values_to = "values")

indicators <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/country_indicators.csv")


met_ind <- full_join(metrics, indicators, by = "country")
met_ind$HM_total_cat <- ordered(met_ind$HM_total_cat, levels = c("low", "medium","high"))
met_ind$EIOS_total_cat <- ordered(met_ind$EIOS_total_cat, levels = c("low", "medium","high"))

metrics_HM <- filter(met_ind, source == "HealthMap") %>% 
  select(-c(X.x, X.y, FluNet_total_cat, FluNet_total, FluNet_max, EIOS_total_cat, EIOS_total, EIOS_max, ISO, 
            influenza_transmission_zone, ISO, problematic_FluNet, problematic_EIOS, problematic_HM, language)) %>% 
  mutate(latitude = abs(latitude))
metrics_EIOS <- filter(met_ind, source == "EIOS") %>% select(-c(X.x, X.y))

##### univariable regressions#####
### HealthMap
coef <- c(1, 2, 5, 7, 9, 10, seq(from = 13, to = 23, by = 2))
pval <- c(3, 4, 6, 8, 11, 12, seq(from = 14, to = 24, by = 2))
# sensitivity per outbreak
lm_HM_sensOB <- lapply(10:19, function(x) lm(metrics_HM$sens_per_outbreak ~ metrics_HM[,x]))
lm_HM_sensOB_summary <- lapply(lm_HM_sensOB, summary)
lm_HM_sensOB_coef_list <- lapply(lm_HM_sensOB_summary, function(x) x$coefficients[-1, c(1,4)])
lm_HM_sensOB_coef_vec <- unlist(lm_HM_sensOB_coef_list)
lm_HM_sensOB_coef <- data.frame("coefficient" = lm_HM_sensOB_coef_vec[coef], "p-value" = lm_HM_sensOB_coef_vec[pval])
row.names(lm_HM_sensOB_coef) <- c("HM_total_cat.L", "HM_total_cat.Q", "HM_total", "HM_max", "global_region.temp_South", 
                                  "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_HM_sensOB_coef, 3)

lm_HM_sensOB_res <- data.frame(sapply(lm_HM_sensOB_summary, residuals))
names(lm_HM_sensOB_res) <- names(metrics_HM)[10:19]

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

lm_HM_sens_exact_res <- data.frame(sapply(lm_HM_sens_exact_summary, residuals))
names(lm_HM_sens_exact_res) <- names(metrics_HM)[10:19]

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

lm_HM_sens_week_res <- data.frame(sapply(lm_HM_sens_week_summary, residuals))
names(lm_HM_sens_week_res) <- names(metrics_HM)[10:19]

par(mfrow = c(2, 2))
sapply(lm_HM_sens_week, plot) 


# specificity 
lm_HM_spec <- lapply(10:19, function(x) lm(metrics_HM$specificity ~ metrics_HM[,x]))
lm_HM_spec_summary <- lapply(lm_HM_spec, summary)
lm_HM_spec_coef_list <- lapply(lm_HM_spec_summary, function(x) x$coefficients[-1, c(1,4)])
lm_HM_spec_coef_vec <- unlist(lm_HM_spec_coef_list)
lm_HM_spec_coef <- data.frame("coefficient" = lm_HM_spec_coef_vec[coef], "p-value" = lm_HM_spec_coef_vec[pval])
row.names(lm_HM_spec_coef) <- c("HM_total_cat.L", "HM_total_cat.Q", "HM_total", "HM_max", "global_region.temp_South", 
                                     "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", "TIU.2017")
round(lm_HM_spec_coef, 3)

lm_HM_spec_res <- data.frame(sapply(lm_HM_spec_summary, residuals))
names(lm_HM_spec_res) <- names(metrics_HM)[10:19]

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

lm_HM_frac_prev_res <- data.frame(sapply(lm_HM_frac_prev_summary, residuals))
names(lm_HM_frac_prev_res) <- names(metrics_HM)[10:19]

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

lm_HM_PPV_res <- data.frame(sapply(lm_HM_PPV_summary, residuals))
names(lm_HM_PPV_res) <- names(metrics_HM)[10:19]

par(mfrow = c(2, 2))
sapply(lm_HM_PPV, plot) 




#### regressions with interactions ##### 
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
# no significant interactions 
