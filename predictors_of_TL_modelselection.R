"
Summary:

B. Significant predictors
  5) Multicollinearity
      A. LASSO model (with interaction terms): full model regularized with LASSO
          - lambda 1se predictors: level_of_precipitation*atm_pressure, pop_size_meancentered, temperature, NAO
      B. No-interaction model: family-wise inference
          - significant predictor(s): temperature only
      C. ANOVA to compare LASSO model and no-interaction model
          - p-value = 0.10 => temperature is the only significant predictor
  6) Diagnostic plots
      - QQ plot
      - residuals v. fitted
      - residuals v. temperature
      - Cook's
      - independence of residuals (plot against date, observation no, other 5 ommitted var)
  7) re-transform back into original units
      
"
#############################################################################################################
dev.off()
rm(list = ls())
# Libraries
library(corrplot)
library(glmnet)

# Functions
#############################################################################################################
# Accessing data
tl_df <- read.csv("./data_clean/predictors_of_tl_dataset.csv")
head(tl_df)
#############################################################################################################
# checking for multicollinearity
full_pred <- c("island_name",
                      "pop_size_meancentered",
                      "temperature*level_of_precipitation",
                      "level_of_precipitation * atm_pressure", "NAO"
)
full_formula <- as.formula(
  paste("bc_TL", 
        paste(full_pred, collapse = " + "),
        sep = " ~ ")
)
X <- model.matrix(lm(full_formula, data = tl_df), data = tl_df)[, -1]
corrplot(cor(X), type = "lower", tl.cex = 0.65)
# High multicollinearity among precipitation variables. Regularize with LASSO later.

# minus interactions
pred_nointer <- c("island_name",
                  "pop_size_meancentered",
                  "temperature", "level_of_precipitation",
                  "atm_pressure", "NAO"
)
nointer_formula <- as.formula(
  paste("bc_TL", 
        paste(pred_nointer, collapse = " + "),
        sep = " ~ ")
)
X_nointer <- model.matrix(lm(nointer_formula, data = tl_df), data = tl_df)[, -1]
corrplot(cor(X_nointer), type = "lower", tl.cex = 0.65)
# looks much better
####################################################################################################################
# LASSO Regularization on full model
X <- model.matrix(full_formula, data = tl_df)
cvfit <- cv.glmnet(x = X, y = tl_df$bc_TL, alpha = 1)
print(cvfit)
coef(glmnet(X, tl_df$bc_TL), s = 0.001070)
# level_of_precipitation*atm_pressure, pop_size_meancentered, temperature, NAO
lasso_model <- lm(bc_TL ~ level_of_precipitation*atm_pressure
                  + pop_size_meancentered + temperature + NAO,
                  data = tl_df)

# t test model
nointer_model <- lm(nointer_formula, data = tl_df)
bonferroni_a <- 0.05/6
bonferroni_a
# 0.008333333
summary(nointer_model)
# temperature is the only significant predictor
reduced_model <- lm(bc_TL ~ temperature, data = tl_df)
####################################################################################################################
# ANOVA comparison
anova(lasso_model, reduced_model)
# Does not appear to be significantly better
best_model <- reduced_model
####################################################################################################################
# Diagnostics
par(mfrow = c(2,2))
plot(best_model, which = c(1, 2, 4))
plot(tl_df$temperature, best_model$residuals)

# Independence of residuals
par(mfrow = c())
# v. hatching time
plot(best_model$residuals, tl_df$year)
plot(best_model$residuals, tl_df$yearday)
# v. telomere collection time

# v. ID

# v. ommitted predictors
omitted_pred <- c("island_name",
                  "pop_size_meancentered",
                  "level_of_precipitation",
                  "atm_pressure", 
                  "NAO"
)

for (pred in omitted_pred) {
  continuoustypeof(tl_df[[pred]])
  if (
  plot(tl_df[[pred]], best_model$residuals, main = pred)
  }
plot(best_model$residuals, tl_df$year)
independence of residuals (plot against day_of_the_month, year, observation no, other 5 ommitted var)
head(tl_df)
