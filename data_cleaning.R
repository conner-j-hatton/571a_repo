library(dplyr)

# unzip files
zip_f <- "./doi_10.5061_dryad.612jm6463__v2.zip"
out_dir <- "./data"
unzip(zip_f, exdir=out_dir)

# clean data
remove_na <- function(filename){
    source_dir <- "./data"
    file_in <- paste(paste(source_dir, filename, sep = "/"), ".csv", sep = "")
    data <- read.csv(file_in, sep = ";")
    
    # clean data
    data_clean <- na.omit(data)
    out_dir <- "./data_clean"
    file_out <- paste(paste(out_dir, filename, sep = "/"), "clean.csv", sep = "")
    write.csv(x = data_clean, file = file_out)
}

remove_na("reproduction_data_Pepke_etal_EcolEvol2022")
remove_na("sparrow_data_Pepke_etal_EcolEvol2022")
remove_na("weather_data_Pepke_etal_EcolEvol2022")


# check proportion omitted
compare_na <- function(filename){
    file_unclean <- paste("./data/", filename, ".csv", sep = "")
    file_clean <- paste("./data_clean/", filename, "clean", sep = "")
    data <- read.csv(file_unclean)
    clean_data <- read.csv(file_clean)
    paste(c(nrow(data), nrow(clean_data)))
}

compare_na("reproduction_data_Pepke_etal_EcolEvol2022")
compare_na("sparrow_data_Pepke_etal_EcolEvol2022")
compare_na("weather_data_Pepke_etal_EcolEvol2022")