"
Summary: Wrangling data so that LINE assumptions are met

A. Model assumptions
  1) Full model finds nonlinearity (QQ) and nonconstant residuals (vs. ‘precipitation’)
  2) Box-cox transform on Yi produces best diagnostic plots (QQ/resid)
  3) Bin precipitation (none, low, high) to drop constant variance assumption of precipitation predictor v. resid
     - check changes in box-cox, QQ, and resid plots
     - Fixed: all continuous resid plots good
  4) remove outliers (DFFITS, DFBETAS, and Cook's leverage)         
     - check changes in box-cox (new lambda), QQ, and resid plots
        - still better than log transform
     - fixed: QQ plot looks normal now
     - write into csv: 'predictors_of_tl_dataset.csv'
"

############################################################################################################
dev.off()
rm(list = ls())
# Libraries
library(dplyr)
library(MASS)

# Functions
getVarName <- function(object, name = deparse(substitute(object))) {
  return(name)
}

plot_multi_pvr <- function(model, data, pvr_predictors,    # plotting predictors v. residuals (pvr) 
                           change_layout = T,                                # optional arg
                           model_name = deparse(substitute(model)),       # for jpeg naming
                           data_name = deparse(substitute(data))) {        # for jpeg naming
  if (change_layout == T){
    # jpeg part
    path = "./plots/"
    jpegtitle <- paste(path, getVarName(model, model_name), getVarName(data, data_name), "pvr_plot.jpg", sep = '_')
    jpeg(jpegtitle)
    print(jpegtitle) # sanity check
    
    par(mfrow = c(2, ceiling( length(pvr_predictors ) / 2 ) ) )
  }
  
  resid <- model$residuals
  for (pred in pvr_predictors){
    predictor <- data[[pred]]
    plot(predictor, resid, xlab = pred)
  }
  if (change_layout == T){dev.off()}
}

compare_qq_pvr <- function(responses, model_pred, pvr_predictors, data, name){  # compare qq and pvr between Yi transformations
  # jpeg part
  path = './plots/'
  jpegtitle <- paste("Yi_transforms_plot_", name, ".jpg", sep = "")
  jpeg(jpegtitle)
  print(jpegtitle) # sanity check
  
  par( mfrow = c( length(responses) , length(pvr_predictors) + 1 ) )
  for (response in responses) {
    predictors <- paste(model_pred, collapse = " + ")
    formula_str <- paste(response, predictors, sep = " ~ ")
    formula <- as.formula(formula_str)
    model <- lm(formula, data = data)

    plot(model, 2, main = response) # QQ plot

    plot_multi_pvr(model, data, pvr_predictors, change_layout = F)
  }
  dev.off()
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

getDFBETAS <- function(model, data){
  n <- nrow(data)
  p <- length(model$coefficients)
  X <- model.matrix(model)
  original_formula <- eval(model$call[[2]])
  
  DFBETAS <- as.data.frame(matrix(0, nrow=n, ncol=p, byrow=T))
  names(DFBETAS) <- names(model$coefficients)
  rownames(DFBETAS) <- 1:n 
  
  for_c <- solve(t(X)%*%X)
  for (i in 1:n){
    diff_fit <- lm(original_formula, data = data[-c(i),])
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
percentile_conversion <- function(vector) {           # converts a vector of numerical values into percentiles
  percentile <- ecdf(vector)
  n <- length(vector)
  
  converted_vector <- rep(0, n)
  for (i in 1:n) {
    value <- vector[i]
    converted_vector[i] <- percentile(value)
  }
  return(converted_vector)
}

###########################################################################################################
# Accessing data
df_sp <- read.csv("./data_clean/sparrow_data_Pepke_etal_EcolEvol2022clean.csv") %>%
  dplyr::select("ID", "TL", "island_name", "pop_size_meancentered", "hatch_date", "date") # IMPORTANT: this "date" represents telomere sampling date
names(df_sp)[6] <- "sample_date" # rename "date" to "sample_date"

df_wt <- read.csv("./data_clean/weather_data_Pepke_etal_EcolEvol2022clean.csv") %>%
  dplyr::select("date", "year", "yearday","temperature":"NAO")

# Merge data.frames by hatch_date
df_merged <- df_sp %>% merge(df_wt, by.x = "hatch_date", by.y = "date")   # since we are investigating early life predictors
tl_predictors <- names(df_merged)[-c(1:3, 6:8)]
#############################################################################################################
# Full model
full_pred <- c("island_name",
               "pop_size_meancentered",
               "temperature*precipitation", 
               "precipitation * atm_pressure", "NAO"
)

tl_formula <- as.formula(
  paste("TL",
        paste(full_pred, collapse = " + "),
        sep = " ~ ")
)

tl_fit <- lm(tl_formula, data = df_merged)

# Diagnostics plots 
par(mfrow = c(2,2))
plot(tl_fit)                         # 1) QQ: heavy upper tail 
dev.off()

pvr_pred <- tl_predictors[-1]
plot_multi_pvr(model = tl_fit, data = df_merged, pvr_pred = pvr_pred)      # 2) heteroskedastic: precipitation v. residuals

###############################################################################################################
# Non-Normality and Heteroskedasticity: Log v. Root v. Box-Cox on Response

# Box-Cox
bc <- boxcox(tl_fit)
lambda <- bc$x[which.max(bc$y)] # to choose value with maximum likelihood
dev.off()

# Merge transformed responses
df_merged <- df_merged %>% mutate(logTL = log(TL), bc_TL = TL^lambda)

# Compare QQ-plots and residuals
responses <- c("TL","logTL", "bc_TL")
comparepvr_pred <- c("precipitation") # for heteroskedasticity 
compare_qq_pvr(responses = responses, model_pred = full_pred, 
               pvr_predictors = comparepvr_pred, data = df_merged, name = 'raw')   
# box-cox is best
###############################################################################################################
# Refactoring Precipitation

# Binning to drop model's assumption that precipitation v. residuals are normally distr.
precip_subset <- df_merged[df_merged$precipitation!=0,]
med <- median(df_merged$precipitation)
df_merged <- df_merged %>% 
  mutate(level_of_precipitation = 
           case_when(precipitation == 0 ~ "none",
                     precipitation < med ~ "low",
                     precipitation >= med ~ "high"
           )
  )
# check if box-cox is the same?
binnedmodel_pred <- c("island_name",
                      "pop_size_meancentered",
                      "temperature*level_of_precipitation",
                      "level_of_precipitation * atm_pressure", "NAO"
)
tl_binned_formula <- as.formula(
  paste("TL", 
        paste(binnedmodel_pred, collapse = " + "),
        sep = " ~ ")
)

tl_binnedfit <- lm(tl_binned_formula, data = df_merged)
bc2 <- boxcox(tl_binnedfit)
lambda_binned <- bc2$x[which.max(bc2$y)] # to choose value with maximum likelihood
c(lambda, lambda_binned) 
# check: same bc lambda
dev.off()

pvr_pred2 <- c("pop_size_meancentered", "temperature", "atm_pressure", "NAO")
compare_qq_pvr(responses = responses, model_pred = binnedmodel_pred, 
               pvr_predictors = pvr_pred2, data = df_merged, name = 'binned')
# check: similar diagnostics as earlier
####################################################################################################################
# Possible Outliers - using bc model with binned precipitation
bin_bcformula <- as.formula(
  paste("TL^lambda",
        paste(binnedmodel_pred, collapse = " + "),
        sep = " ~ ")
)
bin_bcmodel <- lm(bin_bcformula, data = df_merged)

plot(bin_bcmodel, 4)


# Cook's
largestCooks <- sort(cooks.distance(bin_bcmodel), decreasing = T)[1:20] # top 20 values: < 1% of data
cooks_i <- as.integer(names(largestCooks))
cooks_i          
# 685 2020 1911 1632  935 1818  877  180 1581 2290 2051  227  418 1522  815  181 1693 1884  724 1536

# DFFITS
DFFITS_out <- getDFFITS(bin_bcmodel, df_merged)
sort(abs(DFFITS_out), decreasing = T)[1:100] # magnitude > 0.2 appears abnormal
pos_extreme_DFFITS <- sort(DFFITS_out, decreasing = T)[1:9]
neg_extreme_DFFITS <- sort(DFFITS_out, decreasing = F)[1:13]
largestDFFITS <- c(pos_extreme_DFFITS, neg_extreme_DFFITS) # magnitude > 0.2

DFFITS_i <- match(largestDFFITS, DFFITS_out)
DFFITS_i
# 685 2020  935 1581  227 1522 1693  724 1094 1911 1632 1818  877  180 2290 2051  418  815  181 1884 1536
# 999

DFBETAS_out <- getDFBETAS(bin_bcmodel, df_merged) 
sort(abs(unlist(DFBETAS_out, use.names=FALSE)), decreasing = T)[1:100] # magnitude > .2 appear to be abnormal
pos_extreme_DFBETAS <- sort(unlist(DFBETAS_out, use.names=FALSE), decreasing = T)[1:4]
neg_extreme_DFBETAS <- sort(unlist(DFBETAS_out, use.names=FALSE), decreasing = F)[1:5]
largestDFBETAS <- c(pos_extreme_DFBETAS, neg_extreme_DFBETAS)

DFBETAS_i_df <- data.frame(DFBETAS_out, removed_obs = 1:nrow(df_merged))
largestDFBETAS_df <- DFBETAS_i_df %>%
  filter(if_any(.cols = everything(),
                .fns = ~ .x %in% largestDFBETAS)
  )
DFBETAS_i <- largestDFBETAS_df$removed_obs
DFBETAS_i 
# DFBETAS: 685 1632 1818 2020
############################################################################################################################
# Removing outliers
possible_outliers <- unique(c(cooks_i, DFFITS_i, DFBETAS_i))

# convert variables to percentile values to investigate strange behavior
quant_variables <- c("bc_TL", "pop_size_meancentered", "temperature",
                     "precipitation", "atm_pressure", "NAO")
percentile_df = matrix(nrow=nrow(df_merged), ncol = 0)
for (variable in quant_variables) {
  converted_var <- percentile_conversion(df_merged[[variable]])
  percentile_df <- cbind(percentile_df, converted_var)
}
colnames(percentile_df) <- quant_variables
converted_df <- data.frame(percentile_df, island_name = df_merged$island_name, obs_no = 1:nrow(df_merged))
converted_df[possible_outliers,]

# Outlying behavior:  > 99th percentile and < 1st percentile
outliers_df <- converted_df[possible_outliers,] %>%
  filter(if_any(.cols = !(c(island_name, obs_no)),
                .fns = ~ .x > .99 | .x < .01)
  )                   
nrow(outliers_df)/2032 # Remove 1% of data

outliers <- outliers_df[['obs_no']]
df_merged_noout <- df_merged[-outliers,]
############################################################################################################################
# check that box-cox and diagnostics are same
tl_binnedfit_noout <- lm(tl_binned_formula, data = df_merged_noout)
bc3 <- boxcox(tl_binnedfit_noout)
lambda_noout <- bc3$x[which.max(bc3$y)] # to choose value with maximum likelihood
c(lambda, lambda_binned, lambda_noout) 
# check: different bc
dev.off()

# retransform box-cox
df_merged_noout <- df_merged_noout %>% mutate(bc_TL = TL^lambda_noout)

# check: QQ is normal and residuals are constant-looking
pvr_pred_noout <- pvr_pred2
compare_qq_pvr(responses = responses, model_pred = binnedmodel_pred, 
               pvr_predictors = pvr_pred_noout, data = df_merged_noout, name = "no_outliers_binned")
#############################################################################################################################
# Write csv file
file_out <- "./data_clean/predictors_of_tl_dataset.csv"
write.csv(x = df_merged_noout, file = file_out)

