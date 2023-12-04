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
library(ggplot2)

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
dev.off()

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

corrplot(cor(X_nointer), type = "lower", tl.cex = 0.8)
# looks much better
dev.off()
####################################################################################################################
# LASSO Regularization on full model
X <- model.matrix(full_formula, data = tl_df)
cvfit <- cv.glmnet(x = X, y = tl_df$bc_TL, alpha = 1)

plot(cvfit)
dev.off()

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
# temperature is the only significant predictor p-val = 7.45e-05
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
plot(x = tl_df$temperature, y = best_model$residuals, main = "Residuals v. Temperature", xlab = "Temperature")
dev.off()

# Independence of residuals
par(mfrow = c(3,3))
# v. hatching time
plot(y = best_model$residuals, x = tl_df$year, main="hatch year", xlab = "hatch year", ylab="residuals")
plot(y = best_model$residuals, x = tl_df$yearday, main="hatch day of year", xlab = "hatch day of year (1-365)", ylab="residuals")
# v. telomere collection time
tl_df <- tl_df %>% mutate(sample_year = as.integer(substr(sample_date, 7, 10))) %>%
  mutate(sample_month = as.integer(substr(sample_date, 4, 5)))

plot(y=best_model$residuals, x=tl_df$sample_year, main = "telomere collection year",xlab="year", ylab = "residuals")
plot(y=best_model$residuals, x=tl_df$sample_month, main = "telomere collection month", xlab="month",ylab="residuals")

# v. ommitted predictors
cts_omitted_pred <- c("pop_size_meancentered",
                  "atm_pressure", 
                  "NAO"
)
cat_omitted_pred <- c("island_name", "level_of_precipitation")

for (pred in cts_omitted_pred) {
  plot(tl_df[[pred]], best_model$residuals, main = pred, xlab = pred)
}

for (pred in cat_omitted_pred) {
  plot(as.factor(tl_df[[pred]]), best_model$residuals, main = pred)
}

dev.off()
###################################################################################################################
# plotting best model (box-cox transformed)
predictions <- data.frame(temperature = 270:300, bc_TL = predict(best_model, newdata = data.frame(temperature = 270:300)))

ggplot() +
  geom_point(aes(x=temperature, y = bc_TL, color = "blue"), data = tl_df) +
  geom_line(aes(x=temperature, y = bc_TL, col = "darkblue"), data = predictions) +
  ggtitle("Temperature v. Transformed Telomere Length") +
  scale_x_continuous("Temperature (K)") +
  scale_y_continuous("TL^(3/4)") 
###################################################################################################################
# plot out  back-transformed predictions
lambda <- .75
predictions <- predictions %>% 
  mutate(TL = (bc_TL)^(1/lambda))
ggplot() +
  geom_point(aes(x=temperature, y = TL, color = "blue"), data = tl_df) +
  geom_line(aes(x=temperature, y = TL, col = "darkblue"), data = predictions) +
  ggtitle("Temperature v. Telomere Length") +
  scale_x_continuous("Temperature (K)") +
  scale_y_continuous("TL") 

est <- lm(TL ~ temperature, data = predictions)
summary(est)