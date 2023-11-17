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
# 1) heavy upper tail QQ-plot 
# 2) heteroskedastic precipitation v. resids plot
par(mfrow = c(2,2))
plot(tl_fit)
plot(df_merged) # only issue is QQ-plot

# predictors v. residuals
names(df_merged)
par(mfrow = c(3,2))
plot(as.factor(df_merged$island_name), tl_fit$residuals)
plot(df_merged$pop_size_meancentered,tl_fit$residuals)
plot(df_merged$temperature,tl_fit$residuals)
plot(df_merged$precipitation,tl_fit$residuals) # heteroskedastic
plot(df_merged$atm_pressure,tl_fit$residuals)
plot(df_merged$NAO,tl_fit$residuals)
dev.off()

# log of the tl
df_merged <- df_merged %>% mutate(logTL = log(TL))
logTL_formula <- logTL ~ island_name +
  pop_size_meancentered + temperature +
  precipitation + atm_pressure + NAO
logTL_fit <- lm(logTL_formula, data = df_merged)
summary(logTL_fit)
par(mfrow = c(2,2))
plot(logTL_fit) # fixed QQ plot :)

# investigate QQ outliers 685, 2020, 1911
df_merged[1911,]
df_merged[1911,]



# diagnostics
par(mfrow = c(3,2))
plot(as.factor(df_merged$island_name), logTL_fit$residuals)
plot(df_merged$pop_size_meancentered,logTL_fit$residuals)
plot(df_merged$temperature,logTL_fit$residuals)
plot(df_merged$precipitation,logTL_fit$residuals) # heteroskedastic
plot(df_merged$atm_pressure,logTL_fit$residuals)
plot(df_merged$NAO,logTL_fit$residuals)

# Model selection: vif and lasso

