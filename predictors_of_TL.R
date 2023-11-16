library(dplyr)
names(df)
df_sp <- read.csv("./data_clean/sparrow_data_Pepke_etal_EcolEvol2022clean") %>%
    select(c("ID", "TL", "island_name", "pop_size_meancentered", "hatch_date"))
# mutate full_date in weather date
df_wt <- read.csv("./data_clean/weather_data_Pepke_etal_EcolEvol2022clean") %>%
    mutate()

df_wt$
df_merged <- df %>% left_join()


