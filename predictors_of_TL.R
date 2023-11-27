############################################################################################################
# Libraries
library(dplyr)
library(MASS)

# Functions
plot_pvr <- function(model, df, pvr_pred) {     # plotting predictors v. residuals
  resid <- model$residuals
    for (pred in pvr_pred){
      predictor <- df[[pred]]
      plot(predictor, resid, xlab = pred)
  }
}

getDFFITS <- function(model, data){
  n <- nrow(data)
  p <- length(model$coefficients)
  SSE <- sum(model$residuals^2)
  X <- model.matrix(model)
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  
  e <- rep(0,n)
  t <- rep(0,n)
  DFFITS <- rep(0,n)
  
  for (i in 1:n) {
    e[i] <- model$residuals[i]
    t[i] <- e[i] * ( (n-p-1) / ( SSE*(1-H[i,i]) - e[i]^2 ) )^(1/2)
    DFFITS[i] <- t[i] * ( (H[i,i]) / (1 - H[i,i]) )^(1/2)
  }
  return(DFFITS)
}

getDFBETAS <- function(model,df){
  n <- nrow(df)
  p <- length(model$coefficients)
  X <- model.matrix(model)
  original_formula <- eval(model$call[[2]])
  
  DFBETAS <- as.data.frame(matrix(0, nrow=n, ncol=p, byrow=T))
  names(DFBETAS) <- names(model$coefficients)
  rownames(DFBETAS) <- 1:n 
  
  for_c <- solve(t(X)%*%X)
  for (i in 1:n){
    diff_fit <- lm(original_formula, data = df[-c(i),])
    MSE_diff <- mean(diff_fit$residuals^2)
    
    for (k in 1:p) {
      b_k <- model$coefficients[[k]]
      b_k_diff <- diff_fit$coefficients[[k]]
      c_kk <- for_c[k,k]
      DFBETAS[i,k] <- (b_k - b_k_diff) / sqrt(MSE_diff * c_kk)
    }
  }
  return(DFBETAS)
}
###########################################################################################################
# Accessing data
df_sp <- read.csv("./data_clean/sparrow_data_Pepke_etal_EcolEvol2022clean.csv") %>%
  dplyr::select("ID", "TL", "island_name", "pop_size_meancentered", "hatch_date")

df_wt <- read.csv("./data_clean/weather_data_Pepke_etal_EcolEvol2022clean.csv") %>%
  dplyr::select("date","temperature":"NAO")

# Merge data.frames by hatch_date/date
df_merged <- df_sp %>% merge(df_wt, by.x = "hatch_date", by.y = "date")
#############################################################################################################
# Full model
tl_formula <- TL ~ island_name +
  pop_size_meancentered + temperature +
  precipitation + atm_pressure + NAO
tl_fit <- lm(tl_formula, data = df_merged)

# Diagnostics plots 
par(mfrow = c(2,2))
plot(tl_fit)                         # 1) QQ: heavy upper tail 
plot_pvr('df_merged', 'tl_fit')      # 2) heteroskedastic: precipitation v. residuals

dev.off()
###############################################################################################################
# Non-Normality and Heteroskedasticity: Log v. Root v. Box-Cox on Response

# Box-Cox
bc <- boxcox(tl_fit)
lambda <- bc$x[which.max(bc$y)] # to choose value with maximum likelihood

# Merge transformed response
df_merged <- df_merged %>% mutate(logTL = log(TL), bc_TL = TL^lambda)

# Compare QQ-plots and residuals
responses <- c("TL","logTL", "bc_TL")

pvr_pred <- c("precipitation") # for heteroskedasticity 

par(mfrow = c(3,2))
compare_qq_pvr(model, df, pvr_pred, cts_pred, cat_pred, responses){
  for (response in responses) {
    predictors <- paste("island_name", "pop_size_meancentered",
                        "temperature", "precipitation", "atm_pressure", "NAO", sep = " + ")
    formula_str <- paste(response, predictors, sep = " ~ ")
    formula <- as.formula(formula_str)
    model <- lm(formula, data = df_merged)
    
    plot(model, 2, main = response) # QQ plot
    
    plot_pvr(model, df_merged, cts_pred, cat_pred) # precip vs resid plot
  }
}

plot_qq_pvr

dev.off()
###############################################################################################################
# 0 precipitation seems to have a different behavior
# Separate into 2 subsets based on precipitation 1) 0 and 2) non-0
noprecip_subset <- df_merged[df_merged$precipitation==0,]
precip_subset <- df_merged[df_merged$precipitation!=0,]


par(mfrow = c(3,2))
for (response in responses) {
  predictors <- paste("island_name", "pop_size_meancentered",
                      "temperature", "precipitation", "atm_pressure", "NAO", sep = " + ")
  formula_str <- paste(response, predictors, sep = " ~ ")
  formula <- as.formula(formula_str)
  model <- lm(formula, data = df_merged)
  
  plot(model, 2, main = response) # QQ plot
  
  plot_pvr(model, df_merged, cts_pred, cat_pred) # precip vs resid plot
}

dev.off()

###############################################################################################################
# Possible outliers
bc_formula <- TL^lambda ~ island_name +
  pop_size_meancentered + temperature +
  precipitation + atm_pressure + NAO
best_model <- lm(bc_formula, df_merged)
plot(best_model, 4) # Cook's distance

DFFITS_out <- getDFFITS(best_model, df_merged)
sort(abs(DFFITS_out), decreasing = T)[1:50] # magnitude > 0.3 appear to be abnormal
pos_extreme_DFFITS <- sort(DFFITS_out, decreasing = T)[1:5]
neg_extreme_DFFITS <- sort(DFFITS_out, decreasing = F)[1:5]
largestDFFITS <- c(pos_extreme_DFFITS[1:3], neg_extreme_DFFITS[1:2]) # magnitude > 0.3

DFFITS_i <- match(largestDFFITS, DFFITS_out) # 685 1693 2020  612 1911

df_merged[DFFITS_i, c("bc_TL", "island_name", "pop_size_meancentered", "temperature",
                      "precipitation", "atm_pressure", "NAO")]

DFBETAS_out <- getDFBETAS(best_model, df_merged)
sort(abs(unlist(DFBETAS_out, use.names=FALSE)), decreasing = T)[1:50]
pos_extreme_DFBETAS <- sort(unlist(DFBETAS_out, use.names=FALSE), decreasing = T)[1:5]
neg_extreme_DFBETAS <- sort(unlist(DFBETAS_out, use.names=FALSE), decreasing = F)[1:5]
largestDFBETAS <- c(pos_extreme_DFBETAS[1:3], neg_extreme_DFBETAS[1:2])

indexed_DFBETAS <- data.frame(DFBETAS_out, removed_obs = 1:nrow(df_merged))
indexed_DFBETAS %>%
  filter(if_any(.cols = everything(),
                .fns = ~ .x %in% largestDFBETAS)
  )
# 3 out of 4 of them are in the precipitation category
####################################################################################################################
# Heteroskedasticity: inverse weighting

