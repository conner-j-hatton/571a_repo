library(dplyr)

df_sp <- read.csv("./data_clean/sparrow_data_Pepke_etal_EcolEvol2022clean.csv") %>%
    select(c("ID", "TL", "island_name", "pop_size_meancentered", "hatch_date"))

df_wt <- read.csv("./data_clean/weather_data_Pepke_etal_EcolEvol2022clean.csv") %>%
    select("date","temperature":"NAO")

# merge by hatch_date/date
df_merged <- df_sp %>% merge(df_wt, by.x = "hatch_date", by.y = "date")

tl_formula <- TL ~ island_name +
                   pop_size_meancentered + temperature +
                   precipitation + atm_pressure + NAO
tl_fit <- lm(tl_formula, data = df_merged)
summary(tl_fit)

# diagnostics plots 
# 1) heavy upper tail QQ-plot - transforms on Yi
# 2) heteroskedastic precipitation v. residuals - separate 0 and non-0 subsets
par(mfrow = c(2,2))
plot(tl_fit)
plot(df_merged) # only issue is QQ-plot

# predictors v. residuals
plot_pvr <- function(df, model) {
  model_resid <- get(model)$residuals
  cts_pred <- c('pop_size_meancentered', 'temperature', 'precipitation', 'atm_pressure', 'NAO')
  
  par(mfrow = c(3,2))
  for (pred in cts_pred){
    predictor <- get(df)[[pred]]
    plot(predictor, model_resid, xlab = pred)
  }
  plot(as.factor(get(df)$island_name), model_resid)
}

plot_pvr('df_merged', 'tl_fit')

# log and cox transform of the tl                  ???????????????????? need to try this
df_merged <- df_merged %>% mutate(logTL = log(TL), rootTL = TL^(1/3))
logTL_formula <- logTL ~ island_name +
  pop_size_meancentered + temperature +
  precipitation + atm_pressure + NAO
rootTL_formula <- logTL^(1/3) ~ island_name +
  pop_size_meancentered + temperature +
  precipitation + atm_pressure + NAO
logTL_fit <- lm(logTL_formula, data = df_merged)
rootTL_fit <- lm(rootTL_formula, data = df_merged)
par(mfrow = c(2,2))
plot(logTL_fit) # fixed QQ plot :)
plot(rootTL_fit)


# heteroskedasticity: look at nonzero precip subset
noprecip_subset <- df_merged[df_merged$precipitation==0,]
precip_subset <- df_merged[df_merged$precipitation!=0,]

noprecip_logTL_fit <- lm(logTL_formula, data = noprecip_subset)
precip_logTL_fit <- lm(logTL_formula, data = precip_subset)

summary(noprecip_logTL_fit)
summary(precip_logTL_fit)

par(mfrow = c(3,2))
plot(as.factor(precip_subset$island_name), precip_logTL_fit$residuals)
plot(precip_subset$pop_size_meancentered, precip_logTL_fit$residuals)
plot(precip_subset$temperature,precip_logTL_fit$residuals)
plot(precip_subset$precipitation,precip_logTL_fit$residuals) 
plot(precip_subset$atm_pressure,precip_logTL_fit$residuals)
plot(precip_subset$NAO,precip_logTL_fit$residuals)
dev.off



# heteroskedasticity: attempt 2 - inverse weighting on Yi
wt <- 1 / lm(abs(logTL_fit$residuals) ~ logTL_fit$fitted.values)$fitted.values^2
wls_fit <- lm(logTL_formula, data = df_merged, weights=wt)
summary(wls_fit)



















# investigate QQ outliers 685, 2020, 1911
df_merged[c(685,2020,1911),]
head(sort(df_merged$TL,decreasing = T))
head(sort(df_merged$TL,decreasing = F))
# 685: shortest TL .14; ID 234
# 2020: 2nd shortest TL .19; ID 230
# 1911: 2nd longest TL 3.76; ID 1891

# looking at the longest TL
df_merged[df_merged$TL == '3.78',]
logTL_fit$residuals[1818] # also had large residuals
logTL_fit$residuals[c(685,2020,1911)]
head(sort(logTL_fit$residuals, decreasing = T))
head(sort(logTL_fit$residuals, decreasing = F))

# check for other missing variables among these possible outliers
'''use id'''





# diagnostics
par(mfrow = c(3,2))
plot(as.factor(df_merged$island_name), logTL_fit$residuals)
plot(df_merged$pop_size_meancentered,logTL_fit$residuals)
plot(df_merged$temperature,logTL_fit$residuals)
plot(df_merged$precipitation,logTL_fit$residuals) # still slightly heteroskedastic
plot(df_merged$atm_pressure,logTL_fit$residuals)
plot(df_merged$NAO,logTL_fit$residuals)

# Model selection: vif and lasso
'''do this next'''
