library(ggplot2)
library(MASS)
library(lubridate)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(binom)
library(ggrepel)
library(gridExtra)
library(grid)
library(cowplot)

metrics <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/metrics.csv", stringsAsFactors = TRUE)
metrics$false_alarm_rate <- 1 - metrics$specificity

metrics_long <- pivot_longer(metrics, cols = c(sens_per_outbreak, sens_per_week, sens_exact, PPV, specificity, frac_prevented), 
                             names_to = "metric", values_to = "values")

indicators <- read.csv("D:/Dokumente (D)/McGill/Thesis/SurveillanceData/data_epidemic/country_indicators.csv")


met_ind <- full_join(metrics, indicators, by = "country")
met_ind$HM_total_cat <- ordered(met_ind$HM_total_cat, levels = c("low", "medium","high"))
met_ind$EIOS_total_cat <- ordered(met_ind$EIOS_total_cat, levels = c("low", "medium","high"))

# select important predictors and log total and max counts
metrics_HM <- filter(met_ind, source == "HealthMap", !country %in% c("Nigeria", "Thailand", "Vietnam")) %>% 
  select(-c(X.x, X.y, FluNet_total_cat, FluNet_total, FluNet_max, EIOS_total_cat, EIOS_total, EIOS_max, ISO, 
            influenza_transmission_zone, ISO, problematic_FluNet, problematic_EIOS, problematic_HM, language)) %>% 
  mutate(latitude = abs(latitude), HM_total = log(HM_total), HM_max = log(HM_max))
# set Nigeria to NA for sens_per_outbreak and frac_prevented because it is an outlier and distorts the regressions
# metrics_HM[metrics_HM$country == "Nigeria", c(3, 8)] <- NA

metrics_EIOS <- filter(met_ind, source == "EIOS", !country %in% c("Nigeria", "Thailand", "Vietnam")) %>% 
  select(-c(X.x, X.y, FluNet_total_cat, FluNet_total, FluNet_max, HM_total_cat, HM_total, HM_max, ISO, 
            influenza_transmission_zone, ISO, problematic_FluNet, problematic_EIOS, problematic_HM, language)) %>% 
  mutate(latitude = abs(latitude), EIOS_total = log(EIOS_total), EIOS_max = log(EIOS_max))
# set counts for USA to NA so that they are disregarded in regressions, because they are outliers
# metrics_EIOS[metrics_EIOS$country == "United States", 11:12] <- NA

country_list <- levels(met_ind$country)


## regressions
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

predictors <- names(metrics_HM)[10:20]
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



HM_comp_df <- read.csv("univariable regressions HM.csv", row.names = 1)

HM_abs_diff_estimate <- reg_HM_coef_df$Estimate - HM_comp_df$Estimate
names(HM_abs_diff_estimate) <- rownames(reg_HM_coef_df)
plot(HM_abs_diff_estimate)
HM_abs_diff_estimate[which(abs(HM_abs_diff_estimate) > 0.1)]

HM_rel_diff_estimate <- reg_HM_coef_df$Estimate/HM_comp_df$Estimate
names(HM_rel_diff_estimate) <- rownames(reg_HM_coef_df)
plot(HM_rel_diff_estimate)
HM_rel_diff_estimate[which(abs(HM_rel_diff_estimate) > 10)]



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


EIOS_comp_df <- read.csv("univariable regressions EIOS.csv", row.names = 1)

EIOS_abs_diff_estimate <- reg_EIOS_coef_df$Estimate - EIOS_comp_df$Estimate
names(EIOS_abs_diff_estimate) <- rownames(reg_EIOS_coef_df)
plot(EIOS_abs_diff_estimate)
EIOS_abs_diff_estimate[which(abs(EIOS_abs_diff_estimate) > 0.1)]

EIOS_rel_diff_estimate <- reg_EIOS_coef_df$Estimate/EIOS_comp_df$Estimate
names(EIOS_rel_diff_estimate) <- rownames(reg_EIOS_coef_df)
plot(EIOS_rel_diff_estimate)
EIOS_rel_diff_estimate[which(abs(EIOS_rel_diff_estimate) > 10)]



### multivariable 
lm_HM_multi <- lapply(3:8, function(x) lm(metrics_HM[,x] ~ metrics_HM$HM_total + metrics_HM$global_region + metrics_HM$english + 
                                            metrics_HM$HDI.2018 + metrics_HM$latitude + metrics_HM$PFI.2018 + metrics_HM$HM_filter_lang))
AIC_HM_multi <- lapply(lm_HM_multi, stepAIC, direction = "both", trace = FALSE)
AIC_HM_multi_summary <- lapply(AIC_HM_multi, summary)
AIC_HM_multi_coef <- lapply(AIC_HM_multi_summary, '[[', "coefficients")
names(AIC_HM_multi_coef) <- names(metrics_HM)[3:8]
AIC_HM_multi_confint <- lapply(AIC_HM_multi, confint)
names(AIC_HM_multi_confint) <- names(metrics_HM)[3:8]
AIC_HM_multi_coef_df <- do.call("rbind", lapply(AIC_HM_multi_coef, as.data.frame))
AIC_HM_multi_coef_df <- AIC_HM_multi_coef_df[-grep("Intercept", rownames(AIC_HM_multi_coef_df)), ]

AIC_HM_multi_coef_comp <- read.csv("regression_coefficients_HM.csv", row.names = 1)
AIC_HM_multi_coef_comp <- AIC_HM_multi_coef_comp[-grep("Intercept", rownames(AIC_HM_multi_coef_comp)), ]

rownames(AIC_HM_multi_coef_comp)
rownames(AIC_HM_multi_coef_df)
# not selected in sensitivity analysis: "sens_per_week.metrics_HM$HDI.2018", "sens_per_week.metrics_HM$HM_filter_langTRUE",
# selected in sensitivity analysis but not in original: "PPV.metrics_HM$englishTRUE", "PPV.metrics_HM$HDI.2018", "PPV.metrics_HM$HM_filter_langTRUE", "specificity.metrics_HM$HM_total","specificity.metrics_HM$englishTRUE"


lm_EIOS_multi <- lapply(3:8, function(x) lm(metrics_EIOS[,x] ~ metrics_EIOS$EIOS_total + metrics_EIOS$global_region + metrics_EIOS$english + 
                                              metrics_EIOS$HDI.2018 + metrics_EIOS$latitude + metrics_EIOS$PFI.2018))
AIC_EIOS_multi <- lapply(lm_EIOS_multi, stepAIC, direction = "both", trace = FALSE)
AIC_EIOS_multi_summary <- lapply(AIC_EIOS_multi, summary)
AIC_EIOS_multi_coef <- lapply(AIC_EIOS_multi_summary, '[[', "coefficients")
names(AIC_EIOS_multi_coef) <- names(metrics_EIOS)[3:8]
AIC_EIOS_multi_confint <- lapply(AIC_EIOS_multi, confint)
names(AIC_EIOS_multi_confint) <- names(metrics_EIOS)[3:8]
AIC_EIOS_multi_coef_df <- do.call("rbind", lapply(AIC_EIOS_multi_coef, as.data.frame))
AIC_EIOS_multi_coef_df <- AIC_EIOS_multi_coef_df[-grep("Intercept", rownames(AIC_EIOS_multi_coef_df)), ]


AIC_EIOS_multi_coef_comp <- read.csv("regression_coefficients_EIOS.csv", row.names = 1)
AIC_EIOS_multi_coef_comp <- AIC_EIOS_multi_coef_comp[-grep("Intercept", rownames(AIC_EIOS_multi_coef_comp)), ]

rownames(AIC_EIOS_multi_coef_comp)
rownames(AIC_EIOS_multi_coef_df)

# not selected in sensitivity analysis: "PPV.metrics_EIOS$EIOS_total" "specificity.metrics_EIOS$global_region"
# selected in sensitivity analysis but not in original: "specificity.metrics_EIOS$HDI.2018"