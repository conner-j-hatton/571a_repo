"
Summary:

B. Significant predictors
  5) multicollinearity?
  6) family-wise inference
  7) Ridge to compare
  8) independent residuals: 
      - checked time predictors (date or day_of_the_year or obs no.)
      - checked if need to add back any predictors

C. Visualization
"
#############################################################################################################
# Libraries
library(corrplot)

# Functions
#############################################################################################################
# Accessing data
tl_df <- read.csv("./data_clean/predictors_of_tl_dataset.csv")
head(tl_df)
#############################################################################################################
# checking for multicollinearity
X <- model.matrix(ripa_formula, data = ripa_2021)[, -1]
corrplot(cor(X), type = "lower", tl.cex = 0.65)