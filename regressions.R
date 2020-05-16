library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(readxl)
library(gvlma)

metrics <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics.csv", stringsAsFactors = TRUE)
metrics$false_alarm_rate <- 1 - metrics$specificity

metrics_long <- pivot_longer(metrics, cols = c(sens_per_outbreak, sens_per_week, sens_exact, PPV, specificity, frac_prevented), 
                             names_to = "metric", values_to = "values")

indicators <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/country_indicators.csv")


met_ind <- full_join(metrics, indicators, by = "country")
met_ind$HM_total_cat <- ordered(met_ind$HM_total_cat, levels = c("low", "medium","high"))
met_ind$EIOS_total_cat <- ordered(met_ind$EIOS_total_cat, levels = c("low", "medium","high"))

# select important predictors and log total and max counts
metrics_HM <- filter(met_ind, source == "HealthMap") %>% 
  select(-c(X.x, X.y, FluNet_total_cat, FluNet_total, FluNet_max, EIOS_total_cat, EIOS_total, EIOS_max, ISO, 
            influenza_transmission_zone, ISO, problematic_FluNet, problematic_EIOS, problematic_HM, language)) %>% 
  mutate(latitude = abs(latitude), HM_total = log(HM_total), HM_max = log(HM_max))
# set Nigeria to NA for sens_per_outbreak and frac_prevented because it is an outlier and distorts the regressions
# metrics_HM[metrics_HM$country == "Nigeria", c(3, 8)] <- NA

metrics_EIOS <- filter(met_ind, source == "EIOS") %>% 
  select(-c(X.x, X.y, FluNet_total_cat, FluNet_total, FluNet_max, HM_total_cat, HM_total, HM_max, ISO, 
            influenza_transmission_zone, ISO, problematic_FluNet, problematic_EIOS, problematic_HM, language)) %>% 
  mutate(latitude = abs(latitude), EIOS_total = log(EIOS_total), EIOS_max = log(EIOS_max))
# set counts for USA to NA so that they are disregarded in regressions, because they are outliers
# metrics_EIOS[metrics_EIOS$country == "United States", 11:12] <- NA

country_list <- levels(met_ind$country)


##### univariable regressions #####

reg.fun <- function(x, y, df){
  reg <- vector(mode = "list", length = 6)
  .reg.fun<- function(x, y, df){
    reg <- lm(df[,y] ~ df[,x])
    return(reg)
  }
  
  for (i in seq_along(y)){
    reg[[i]] <- mapply(.reg.fun, x, y[i], MoreArgs = list(df = df))
  }
  return(reg)
}

reg.summary <- function(x, y, df){
  reg <- vector(mode = "list", length = 6)
  .reg.fun<- function(x, y, df){
    reg <- lm(df[,y] ~ df[,x])
    summary <- summary(reg)
    return(summary)
  }
  
  for (i in seq_along(y)){
    reg[[i]] <- mapply(.reg.fun, x, y[i], MoreArgs = list(df = df))
  }
  return(reg)
}

predictors <- names(metrics_HM[10:20])
outcomes <- names(metrics_HM[3:8])
reg_list_HM <- reg.fun(x = predictors, y = outcomes, df = metrics_HM)
names(reg_list_HM) <- names(metrics_HM[3:8])

reg_HM_confint_list <- vector(mode = "list")
for(i in 1:6){
  reg_HM_confint_list[[i]] <- lapply(reg_list_HM[[i]], confint)
  reg_HM_confint_list[[i]] <- do.call("rbind", lapply(reg_HM_confint_list[[i]], as.data.frame))
}
reg_HM_confint_df <- do.call("rbind", lapply(reg_HM_confint_list, as.data.frame))
reg_HM_confint_df <- reg_HM_confint_df[-grep("Intercept", rownames(reg_HM_confint_df)), ]
reg_HM_confint_df <- round(reg_HM_confint_df, 4)
reg_HM_confint_df <- reg_HM_confint_df %>% unite(col = "confint_95", sep = " - ")


reg_summary_HM <- reg.summary(x = predictors, y = outcomes, df = metrics_HM)

reg_HM_coef_list <- vector(mode = "list")
reg_HM_rsqu_list <- vector(mode = "list")
for(i in 1:6){
  reg_HM_coef_list[[i]] <- do.call("rbind", reg_summary_HM[[i]][4, 1:11])
  reg_HM_rsqu_list[[i]] <- do.call("rbind", reg_summary_HM[[i]][8, 1:11])
}

reg_HM_coef_df <- do.call("rbind", lapply(reg_HM_coef_list, as.data.frame))
reg_HM_coef_df <- reg_HM_coef_df[-grep("Intercept", rownames(reg_HM_coef_df)), ]
reg_HM_coef_df$confint_95 <- reg_HM_confint_df$confint_95 
rownames(reg_HM_coef_df) <- paste(rep(c("HM_total_cat.L", "HM_total_cat.Q", "HM_total", "HM_max", "global_region.temp_South", 
                                        "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", 
                                        "TIU.2017", "HM_filter_lang"), 6), rep(names(metrics_HM)[3:8], each = 13), sep = ".")

write.csv(reg_HM_coef_df, "univariable regressions HM.csv")
write.csv(unlist(reg_HM_rsqu_list), "univariabe regressions HM rsquared.csv", row.names = FALSE)

lapply(reg_list_HM[[1]], plot)
lapply(reg_list_HM[[1]], gvlma)


## EIOS ## 
predictors <- names(metrics_EIOS[10:19])
outcomes <- names(metrics_EIOS[3:8])
reg_list_EIOS <- reg.fun(x = predictors, y = outcomes, df = metrics_EIOS)
names(reg_list_EIOS) <- names(metrics_EIOS[3:8])

reg_EIOS_confint_list <- vector(mode = "list")
for(i in 1:6){
  reg_EIOS_confint_list[[i]] <- lapply(reg_list_EIOS[[i]], confint)
  reg_EIOS_confint_list[[i]] <- do.call("rbind", lapply(reg_EIOS_confint_list[[i]], as.data.frame))
}
reg_EIOS_confint_df <- do.call("rbind", lapply(reg_EIOS_confint_list, as.data.frame))
reg_EIOS_confint_df <- reg_EIOS_confint_df[-grep("Intercept", rownames(reg_EIOS_confint_df)), ]
reg_EIOS_confint_df <- round(reg_EIOS_confint_df, 4)
reg_EIOS_confint_df <- reg_EIOS_confint_df %>% unite(col = "confint_95", sep = " - ")


reg_summary_EIOS <- reg.summary(x = predictors, y = outcomes, df = metrics_EIOS)

reg_EIOS_coef_list <- vector(mode = "list")
reg_EIOS_rsqu_list <- vector(mode = "list")
for(i in 1:6){
  reg_EIOS_coef_list[[i]] <- do.call("rbind", reg_summary_EIOS[[i]][4, 1:10])
  reg_EIOS_rsqu_list[[i]] <- do.call("rbind", reg_summary_EIOS[[i]][8, 1:10])
}

reg_EIOS_coef_df <- do.call("rbind", lapply(reg_EIOS_coef_list, as.data.frame))
reg_EIOS_coef_df <- reg_EIOS_coef_df[-grep("Intercept", rownames(reg_EIOS_coef_df)), ]
reg_EIOS_coef_df$confint_95 <- reg_EIOS_confint_df$confint_95 
rownames(reg_EIOS_coef_df) <- paste(rep(c("EIOS_total_cat.L", "EIOS_total_cat.Q", "EIOS_total", "EIOS_max", "global_region.temp_South", 
                                        "global_region.tropical", "English", "HDI.2018", "abs.latitude", "longitude", "PFI.2018", 
                                        "TIU.2017"), 6), rep(names(metrics_EIOS)[3:8], each = 12), sep = ".")

write.csv(reg_EIOS_coef_df, "univariable regressions EIOS.csv")
write.csv(unlist(reg_EIOS_rsqu_list), "univariabe regressions EIOS rsquared.csv", row.names = FALSE)

lapply(reg_list_EIOS[[1]], plot)
lapply(reg_list_EIOS[[1]], gvlma)


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




#### multivariable models: advanced variable selection ####
library(caret)
library(leaps)
library(MASS)
library(relaimpo)
library(car)


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
AIC_HM_multi_confint <- lapply(AIC_HM_multi, confint)
names(AIC_HM_multi_confint) <- names(metrics_HM)[3:8]


AIC_HM_multi_coef_df <- do.call("rbind", lapply(AIC_HM_multi_coef, as.data.frame))
#write.csv(AIC_HM_multi_coef_df, "regression_coefficients_HM.csv")
AIC_HM_multi_confint_df <- do.call("rbind", lapply(AIC_HM_multi_confint, as.data.frame))
AIC_HM_multi_confint_df <- round(AIC_HM_multi_confint_df, 4)
AIC_HM_multi_confint_df <- AIC_HM_multi_confint_df %>% unite(col = "95%_confint", sep = " - ", remove = TRUE)
#write.csv(AIC_HM_multi_confint_df, "regression_coefficients_confint_HM.csv", quote = FALSE)

AIC_HM_multi_AIC <- sapply(AIC_HM_multi, AIC)
names(AIC_HM_multi_AIC) <- names(metrics_HM)[3:8]

AIC_HM_multi_radj <- sapply(AIC_HM_multi_summary, '[[', "adj.r.squared")
names(AIC_HM_multi_radj) <- names(metrics_HM)[3:8]


AIC_HM_multi_outlier <- lapply(AIC_HM_multi, outlierTest)
AIC_HM_multi_leverageplot <- lapply(AIC_HM_multi, leveragePlots)

AIC_HM_multi_spreadplot <- lapply(AIC_HM_multi, spreadLevelPlot)
AIC_HM_multi_ncvtest <- lapply(AIC_HM_multi, ncvTest)
AIC_HM_multi_ncvtest_p <- sapply(AIC_HM_multi_ncvtest, '[[', "p")

AIC_HM_multi_res <- lapply(AIC_HM_multi, '[[', "residuals")
names(AIC_HM_multi_res) <- names(metrics_HM)[3:8]

par(mfrow = c(1, 1))
AIC_HM_multi_hist_res <- lapply(AIC_HM_multi_res, hist, main = "")
AIC_HM_multi_qqplot <- lapply(AIC_HM_multi, qqPlot)

par(mfrow= c(2, 2))
lapply(AIC_HM_multi, plot)

AIC_HM_multi_crplot <- lapply(AIC_HM_multi, crPlots)

AIC_HM_multi_gvlma <- lapply(AIC_HM_multi, gvlma)

# have a closer look at PPV
model_HM_PPV <- lm(PPV ~ HM_total + global_region + english + 
                                HDI.2018 + latitude + PFI.2018 + HM_filter_lang, 
                              data = metrics_HM)
summary(model_HM_PPV)

### EIOS ###

reduced_model_EIOS_sensOB <- lm(sens_per_outbreak ~ EIOS_total + global_region + english + 
                                HDI.2018 + latitude + PFI.2018, 
                              data = metrics_EIOS)
vif(reduced_model_EIOS_sensOB) # all these predictors don't have collinearity (VIF < 3)
summary(reduced_model_EIOS_sensOB)

lm_EIOS_multi <- lapply(3:8, function(x) lm(metrics_EIOS[,x] ~ metrics_EIOS$EIOS_total + metrics_EIOS$global_region + metrics_EIOS$english + 
                                            metrics_EIOS$HDI.2018 + metrics_EIOS$latitude + metrics_EIOS$PFI.2018))
AIC_EIOS_multi <- lapply(lm_EIOS_multi, stepAIC, direction = "both", trace = FALSE)
AIC_EIOS_multi_summary <- lapply(AIC_EIOS_multi, summary)
AIC_EIOS_multi_coef <- lapply(AIC_EIOS_multi_summary, '[[', "coefficients")
names(AIC_EIOS_multi_coef) <- names(metrics_EIOS)[3:8]
AIC_EIOS_multi_confint <- lapply(AIC_EIOS_multi, confint)
names(AIC_EIOS_multi_confint) <- names(metrics_EIOS)[3:8]


AIC_EIOS_multi_coef_df <- do.call("rbind", lapply(AIC_EIOS_multi_coef, as.data.frame))
#write.csv(AIC_EIOS_multi_coef_df, "regression_coefficients_EIOS.csv")
AIC_EIOS_multi_confint_df <- do.call("rbind", lapply(AIC_EIOS_multi_confint, as.data.frame))
AIC_EIOS_multi_confint_df <- round(AIC_EIOS_multi_confint_df, 4)
AIC_EIOS_multi_confint_df <- AIC_EIOS_multi_confint_df %>% unite(col = "95%_confint", sep = " - ", remove = TRUE)
write.csv(AIC_EIOS_multi_confint_df, "regression_coefficients_confint_EIOS.csv", quote = FALSE)

AIC_EIOS_multi_AIC <- sapply(AIC_EIOS_multi, AIC)
names(AIC_EIOS_multi_AIC) <- names(metrics_EIOS)[3:8]

AIC_EIOS_multi_radj <- sapply(AIC_EIOS_multi_summary, '[[', "adj.r.squared")
names(AIC_EIOS_multi_radj) <- names(metrics_EIOS)[3:8]


AIC_EIOS_multi_outlier <- lapply(AIC_EIOS_multi, outlierTest)
AIC_EIOS_multi_leverageplot <- lapply(AIC_EIOS_multi, leveragePlots)

AIC_EIOS_multi_spreadplot <- lapply(AIC_EIOS_multi, spreadLevelPlot)
AIC_EIOS_multi_ncvtest <- lapply(AIC_EIOS_multi, ncvTest)
AIC_EIOS_multi_ncvtest_p <- sapply(AIC_EIOS_multi_ncvtest, '[[', "p")

AIC_EIOS_multi_res <- data.frame(sapply(AIC_EIOS_multi, '[[', "residuals"))
colnames(AIC_EIOS_multi_res) <- names(metrics_EIOS)[3:8]

par(mfrow = c(1, 1))
AIC_EIOS_multi_hist_res <- apply(AIC_EIOS_multi_res, 2, hist, main = "")
AIC_EIOS_multi_qqplot <- lapply(AIC_EIOS_multi, qqPlot)

par(mfrow= c(2, 2))
lapply(AIC_EIOS_multi, plot)

AIC_EIOS_multi_crplot <- lapply(AIC_EIOS_multi, crPlots)

AIC_EIOS_multi_gvlma <- lapply(AIC_EIOS_multi, gvlma)




## try some stuff

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



##### sens_exact:logistic regression ####
#HealthMap
metrics_HM <- metrics_HM %>% mutate(sens_exact_cat = ifelse(sens_exact == 0, 0 , 1))
glm_HM_sens_exact <- lapply(10:20, function(x) glm(metrics_HM$sens_exact_cat ~ metrics_HM[, x], 
                                                   family = binomial(link = "logit")))
glm_HM_sens_exact_summary <- lapply(glm_HM_sens_exact, summary)
glm_HM_sens_exact_coefs <- lapply(glm_HM_sens_exact_summary, '[[', "coefficients")
names(glm_HM_sens_exact_coefs) <- names(metrics_HM)[10:20]
glm_HM_sens_exact_coef_df <- do.call("rbind", lapply(glm_HM_sens_exact_coefs, as.data.frame))
glm_HM_sens_exact_coef_df <- glm_HM_sens_exact_coef_df[-grep("Intercept", rownames(glm_HM_sens_exact_coef_df)), ]
glm_HM_sens_exact_coef_df$exp.Estimate <- exp(glm_HM_sens_exact_coef_df$Estimate)

glm_HM_sens_exact_confint <- lapply(glm_HM_sens_exact, confint)
names(glm_HM_sens_exact_confint) <- names(metrics_HM)[10:20]
glm_HM_sens_exact_confint_df <- do.call("rbind", lapply(glm_HM_sens_exact_confint, as.data.frame))
glm_HM_sens_exact_confint_df <- glm_HM_sens_exact_confint_df[-grep("Intercept", rownames(glm_HM_sens_exact_confint_df)), ]
glm_HM_sens_exact_confint_df <- apply(glm_HM_sens_exact_confint_df, 2, exp) %>% as.data.frame()
glm_HM_sens_exact_confint_df <- round(glm_HM_sens_exact_confint_df, 4) %>% unite(col = "95%_confint", sep = " - ", remove = TRUE)

glm_HM_sens_exact_coef_df$exp.confint <- glm_HM_sens_exact_confint_df$`95%_confint`
write.csv(glm_HM_sens_exact_coef_df, "reg_coefs_sens_exact_HM_univariable.csv", quote = FALSE)


# multivariable
glm_HM_sens_exact_multi <- glm(sens_exact_cat ~ HM_total + global_region + english + HDI.2018 + latitude + 
                                 PFI.2018 + HM_filter_lang, 
                               data = metrics_HM,
                               family = binomial(link = "logit"))
summary(glm_HM_sens_exact_multi)
glm_HM_sens_exact_multi_AIC <- stepAIC(glm_HM_sens_exact_multi, direction = "both", trace = FALSE)
glm_HM_sens_exact_multi_AIC_summary <- summary(glm_HM_sens_exact_multi_AIC)
glm_HM_sens_exact_multi_AIC_coefs <- data.frame(glm_HM_sens_exact_multi_AIC_summary$coefficients)
glm_HM_sens_exact_multi_AIC_coefs$exp.Estimate <- exp(glm_HM_sens_exact_multi_AIC_coefs$Estimate)

glm_HM_sens_exact_multi_AIC_confint <- confint(glm_HM_sens_exact_multi_AIC)
glm_HM_sens_exact_multi_AIC_confint <- apply(glm_HM_sens_exact_multi_AIC_confint, 2, exp) %>% as.data.frame()
glm_HM_sens_exact_multi_AIC_confint <- round(glm_HM_sens_exact_multi_AIC_confint, 4) %>% unite(col = "95%_confint", sep = " - ", remove = TRUE)

glm_HM_sens_exact_multi_AIC_coefs$confint <- glm_HM_sens_exact_multi_AIC_confint
glm_HM_sens_exact_multi_AIC_coefs <- glm_HM_sens_exact_multi_AIC_coefs[-grep("Intercept", rownames(glm_HM_sens_exact_multi_AIC_coefs)), ]
write.csv(glm_HM_sens_exact_multi_AIC_coefs, "reg_coefs_sens_exact_HM_multivariable.csv", quote = FALSE)


#EIOS
metrics_EIOS <- metrics_EIOS %>% mutate(sens_exact_cat = ifelse(sens_exact == 0, 0 , 1))
glm_EIOS_sens_exact <- lapply(10:20, function(x) glm(metrics_EIOS$sens_exact_cat ~ metrics_EIOS[, x], 
                                                   family = binomial(link = "logit")))
glm_EIOS_sens_exact_summary <- lapply(glm_EIOS_sens_exact, summary)
glm_EIOS_sens_exact_coefs <- lapply(glm_EIOS_sens_exact_summary, '[[', "coefficients")
names(glm_EIOS_sens_exact_coefs) <- names(metrics_EIOS)[10:20]
glm_EIOS_sens_exact_coef_df <- do.call("rbind", lapply(glm_EIOS_sens_exact_coefs, as.data.frame))
glm_EIOS_sens_exact_coef_df <- glm_EIOS_sens_exact_coef_df[-grep("Intercept", rownames(glm_EIOS_sens_exact_coef_df)), ]
glm_EIOS_sens_exact_coef_df$exp.Estimate <- exp(glm_EIOS_sens_exact_coef_df$Estimate)

glm_EIOS_sens_exact_confint <- lapply(glm_EIOS_sens_exact, confint)
names(glm_EIOS_sens_exact_confint) <- names(metrics_EIOS)[10:20]
glm_EIOS_sens_exact_confint_df <- do.call("rbind", lapply(glm_EIOS_sens_exact_confint, as.data.frame))
glm_EIOS_sens_exact_confint_df <- glm_EIOS_sens_exact_confint_df[-grep("Intercept", rownames(glm_EIOS_sens_exact_confint_df)), ]
glm_EIOS_sens_exact_confint_df <- apply(glm_EIOS_sens_exact_confint_df, 2, exp) %>% as.data.frame()
glm_EIOS_sens_exact_confint_df <- round(glm_EIOS_sens_exact_confint_df, 4) %>% unite(col = "95%_confint", sep = " - ", remove = TRUE)

glm_EIOS_sens_exact_coef_df$exp.confint <- glm_EIOS_sens_exact_confint_df$`95%_confint`
write.csv(glm_EIOS_sens_exact_coef_df, "reg_coefs_sens_exact_EIOS_univariable.csv", quote = FALSE)


# multivariable
glm_EIOS_sens_exact_multi <- glm(sens_exact_cat ~ EIOS_total + global_region + english + HDI.2018 + latitude + PFI.2018, 
                               data = metrics_EIOS,
                               family = binomial(link = "logit"))
summary(glm_EIOS_sens_exact_multi)
glm_EIOS_sens_exact_multi_AIC <- stepAIC(glm_EIOS_sens_exact_multi, direction = "both", trace = FALSE)
glm_EIOS_sens_exact_multi_AIC_summary <- summary(glm_EIOS_sens_exact_multi_AIC)
glm_EIOS_sens_exact_multi_AIC_coefs <- data.frame(glm_EIOS_sens_exact_multi_AIC_summary$coefficients)
glm_EIOS_sens_exact_multi_AIC_coefs$exp.Estimate <- exp(glm_EIOS_sens_exact_multi_AIC_coefs$Estimate)

glm_EIOS_sens_exact_multi_AIC_confint <- confint(glm_EIOS_sens_exact_multi_AIC)
glm_EIOS_sens_exact_multi_AIC_confint <- apply(glm_EIOS_sens_exact_multi_AIC_confint, 2, exp) %>% as.data.frame()
glm_EIOS_sens_exact_multi_AIC_confint <- round(glm_EIOS_sens_exact_multi_AIC_confint, 4) %>% unite(col = "95%_confint", sep = " - ", remove = TRUE)

glm_EIOS_sens_exact_multi_AIC_coefs$confint <- glm_EIOS_sens_exact_multi_AIC_confint
glm_EIOS_sens_exact_multi_AIC_coefs <- glm_EIOS_sens_exact_multi_AIC_coefs[-grep("Intercept", rownames(glm_EIOS_sens_exact_multi_AIC_coefs)), ]
write.csv(glm_EIOS_sens_exact_multi_AIC_coefs, "reg_coefs_sens_exact_EIOS_multivariable.csv", quote = FALSE)




#### regressions of postprob cutoff #####
HM_single_cutoff_FAR <- c(1,1,1,0,1,0,0,0,1,0,0,1,1,1,0,1,0,0,0,0,1,1,0,0)
EIOS_single_cutoff_FAR <- c(1,1,0,0,0,1,0,0,1,1,0,0,1,1,1,0,0,1,0,0,0,1,1,0)

glm_HM_cutoff_FAR <- lapply(10:20, function(x) glm(HM_single_cutoff_FAR ~ metrics_HM[, x], 
                                                   family = binomial(link = "logit")))
glm_HM_cutoff_FAR_summary <- lapply(glm_HM_cutoff_FAR, summary)
glm_HM_cutoff_FAR_coefs <- lapply(glm_HM_cutoff_FAR_summary, '[[', "coefficients")
names(glm_HM_cutoff_FAR_coefs) <- names(metrics_HM)[10:20]
glm_HM_cutoff_FAR_coef_df <- do.call("rbind", lapply(glm_HM_cutoff_FAR_coefs, as.data.frame))
glm_HM_cutoff_FAR_coef_df <- glm_HM_cutoff_FAR_coef_df[-grep("Intercept", rownames(glm_HM_cutoff_FAR_coef_df)), ]
glm_HM_cutoff_FAR_coef_df$exp.Estimate <- exp(glm_HM_cutoff_FAR_coef_df$Estimate)


glm_EIOS_cutoff_FAR <- lapply(10:19, function(x) glm(EIOS_single_cutoff_FAR ~ metrics_EIOS[, x], 
                                                   family = binomial(link = "logit")))
glm_EIOS_cutoff_FAR_summary <- lapply(glm_EIOS_cutoff_FAR, summary)
glm_EIOS_cutoff_FAR_coefs <- lapply(glm_EIOS_cutoff_FAR_summary, '[[', "coefficients")
names(glm_EIOS_cutoff_FAR_coefs) <- names(metrics_EIOS)[10:19]
glm_EIOS_cutoff_FAR_coef_df <- do.call("rbind", lapply(glm_EIOS_cutoff_FAR_coefs, as.data.frame))
glm_EIOS_cutoff_FAR_coef_df <- glm_EIOS_cutoff_FAR_coef_df[-grep("Intercept", rownames(glm_EIOS_cutoff_FAR_coef_df)), ]
glm_EIOS_cutoff_FAR_coef_df$exp.Estimate <- exp(glm_EIOS_cutoff_FAR_coef_df$Estimate)

