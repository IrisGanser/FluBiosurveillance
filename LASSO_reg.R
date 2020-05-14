library(tidyverse)
library(dplyr)
library(caret)
library(glmnet)

metrics <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics.csv", stringsAsFactors = TRUE)
metrics$false_alarm_rate <- 1 - metrics$specificity

indicators <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/country_indicators.csv")


met_ind <- full_join(metrics, indicators, by = "country")
met_ind$HM_total_cat <- ordered(met_ind$HM_total_cat, levels = c("low", "medium","high"))
met_ind$EIOS_total_cat <- ordered(met_ind$EIOS_total_cat, levels = c("low", "medium","high"))

# select important predictors and log total and max counts
metrics_HM <- filter(met_ind, source == "HealthMap") %>% 
  dplyr::select(-c(X.x, X.y, FluNet_total_cat, FluNet_total, FluNet_max, EIOS_total_cat, EIOS_total, EIOS_max, ISO, 
            influenza_transmission_zone, ISO, problematic_FluNet, problematic_EIOS, problematic_HM, language)) %>% 
  mutate(latitude = abs(latitude), HM_total = log(HM_total), HM_max = log(HM_max))
# set Nigeria to NA for sens_per_outbreak and frac_prevented because it is an outlier and distorts the regressions
# metrics_HM[metrics_HM$country == "Nigeria", c(3, 8)] <- NA

metrics_EIOS <- filter(met_ind, source == "EIOS") %>% 
  dplyr::select(-c(X.x, X.y, FluNet_total_cat, FluNet_total, FluNet_max, HM_total_cat, HM_total, HM_max, ISO, 
            influenza_transmission_zone, ISO, problematic_FluNet, problematic_EIOS, problematic_HM, language)) %>% 
  mutate(latitude = abs(latitude), EIOS_total = log(EIOS_total), EIOS_max = log(EIOS_max))

country_list <- levels(met_ind$country)



model_matrix_HM <- model.matrix(sens_per_outbreak~HM_total + global_region + english + 
                    HDI.2018 + latitude + PFI.2018 + HM_filter_lang, 
                  data = metrics_HM)[,-1]
model_matrix_EIOS <- model.matrix(sens_per_outbreak ~ EIOS_total + global_region + english + 
                                    HDI.2018 + latitude + PFI.2018, 
                                  data = metrics_EIOS)[,-1]


y <- NA
cv.lasso <- NA

lasso_HM <- lapply(3:8, function(x){
  y <- metrics_HM[, x]
  cv.lasso <- cv.glmnet(model_matrix_HM, y, alpha = 1, family = "gaussian")
  coef <- coef(cv.lasso, cv.lasso$lambda.min)
  return(coef)
})
names(lasso_HM) <- names(metrics_HM)[3:8]

lasso_EIOS <- lapply(3:8, function(x){
  y <- metrics_EIOS[, x]
  cv.lasso <- cv.glmnet(model_matrix_EIOS, y, alpha = 1, family = "gaussian")
  coef <- coef(cv.lasso, cv.lasso$lambda.min)
  return(coef)
})
names(lasso_EIOS) <- names(metrics_EIOS)[3:8]
