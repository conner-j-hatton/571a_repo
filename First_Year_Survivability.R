### Code for Regression Analyses on Effects of Early-life Factors in the Wild House Sparrow ####
### Conner Hatton



###############################################################################
### Analysis on First Year Survivability

#install.packages('tidyverse')
#library(tidyverse)

# import individual sparrow data
sparrow_data_unfiltered <- read.csv("sparrow_data_Pepke_etal_EcolEvol2022.csv")

#remove rows with missing values
sparrow_data <- sparrow_data_unfiltered[complete.cases(sparrow_data_unfiltered),]

# name variables for ease of use
TL = sparrow_data$TL
pop_size = sparrow_data$pop_size_meancentered
island = sparrow_data$island_name
clutch_size = sparrow_data$clutch_size
sex = sparrow_data$sex
tarsus = sparrow_data$tarsus
tarsus_ac = sparrow_data$tarsus_agecorrected
condition = sparrow_data$condition
hatchday_mc = sparrow_data$hatchday_meancentered
fys = sparrow_data$firstyear_survival
disperal = sparrow_data$dispersal
mass = sparrow_data$mass


##### Fit a logistic regression 
logistic_fit <- glm(fys~TL+pop_size+island+clutch_size+sex+tarsus_ac+condition+hatchday_mc, family=binomial(link='logit'))
#print summary from logistic regression fit. 
summary(logistic_fit)


#####perform lasso 

#library(glmnet)
#create a model matrix
X = model.matrix(fys~TL+pop_size+island+clutch_size+sex+tarsus_ac+condition+hatchday_mc)

model_lasso = cv.glmnet(X,fys,family = 'binomial')

best_model_min = glmnet(X,fys,family = 'binomial', lambda = model_lasso$lambda.min)
best_model_1se =glmnet(X,fys,family = 'binomial', lambda = model_lasso$lambda.1se)
best_model_3=glmnet(X,fys,family = 'binomial', lambda = 0.015)
plot(model_lasso)

#### Create Classifier using logistic regression model and fit using a subset of the data

set.seed(1234)

#create a random subset with 80% of the data
train_index = sample(1:nrow(sparrow_data),0.8*nrow(sparrow_data))

#get training and test datasets
train_data = sparrow_data[train_index,]
test_data = sparrow_data[-train_index,]

#rename variables for training
TL_train = TL[train_index]
pop_size_train = pop_size[train_index]
island_train = island[train_index]
clutch_size_train = clutch_size[train_index]
sex_train = sex[train_index]
tarsus_train = tarsus[train_index]
tarsus_ac_train = tarsus_ac[train_index]
condition_train = condition[train_index]
hatchday_mc_train = hatchday_mc[train_index]
fys_train = fys[train_index]

#rename variables for testing
TL_test = TL[-train_index]
pop_size_test = pop_size[-train_index]
island_test = island[-train_index]
clutch_size_test = clutch_size[-train_index]
sex_test= sex[-train_index]
tarsus_test = tarsus[-train_index]
tarsus_ac_test = tarsus_ac[-train_index]
condition_test = condition[-train_index]
hatchday_mc_test = hatchday_mc[-train_index]
fys_test = fys[-train_index]


### Example making predictions
#create data frame
test_dataframe = data.frame(tarsus_ac_train=tarsus_ac_test)
#fit on tarsus age corrected with training set
fit_tarsus_train = glm(fys_train~tarsus_ac_train, family=binomial)
#find the binomial probabilities
predicted_probs = predict(fit_tarsus_train,newdata=test_dataframe,type='response')
#make predictions based on maximizing expected correct 
predictions = ifelse(predicted_probs >=0.5,1,0)
#compute the error rate
error_rate = mean(predictions_1 != fys_test)

# fit on all 7 predictors
fit_7 = glm(fys_train~TL_train+pop_size_train+island_train+clutch_size_train+sex_train+tarsus_ac_train+condition_train+hatchday_mc_train, family=binomial(link='logit'))

test_data_7 = data.frame(TL_train=TL_test, pop_size_train=pop_size_test, island_train=island_test,clutch_size_train=clutch_size_test, sex_train=sex_test,tarsus_ac_train=tarsus_ac_test,condition_train=condition_test,hatchday_mc_train=hatchday_mc_test)

predicted_7_probs = predict(fit_7,newdata=test_data_7,type='response')

predictions_7 = ifelse(predicted_7_probs >=0.5,1,0)

error_rate_7 = mean(predictions_7 != fys_test)

print(1-error_rate_7)

max(predicted_7_probs)


#### Various Plots 

# add gaussian noise to first year survival for display purposes 
fys_noise = fys+rnorm(n=length(fys),mean=0, sd=0.05 )


#plot first-year survival rate
plot(fys_noise,
     main = "First-Year Survival Rate",
     xlab = "Observation",
     ylab = "Survial(1)",cex=0.5
)

#plots for prediction
prediction_vals = seq(min(tarsus_ac),max(tarsus_ac),length.out=1000)
#get bernoulli probabilities from model
p_probs = predict(fit_tarsus_train, newdata = data.frame(tarsus_ac_train=prediction_vals),type='response')
#make plot
plot(tarsus_ac,fys_noise,main = "Logistic Regression Fit on First-Year Survival Rate",
     xlab = "Age-adjusted Tarsus Length",
     ylab = "Survial(1)",cex=0.5)
lines(prediction_vals,p_probs,type='l',col='red', lwd=2)

## Effects of population size on First-year Survivability 

mean_pop_survive = mean(pop_size[survived_first_year_index])
mean_pop_death = mean(pop_size[-survived_first_year_index])
cat("Mean population for First year Survivors:",mean_pop_survive, "Mean Population for Non First Year Surivivors:", mean_pop_death)


plot(pop_size,fys_noise,main = "First-Year Survival Rate and Pop. Size",
     xlab = "Census Population Size (Mean Centered)",
     ylab = "Survial(1)",cex=0.5)
points(mean_pop_survive,1,pch=19,col='green')
points(mean_pop_death,0,pch=19,col='red')


## Effects of Mass on First-year Survivability 

mean_mass_survive = mean(mass[survived_first_year_index])
mean_mass_death = mean(mass[-survived_first_year_index])
cat("Mean mass for First year Survivors:",mean_mass_survive, "Mean Massfor Non First Year Surivivors:", mean_mass_death)


plot(mass,fys_noise,main = "First-Year Survival Rate and Sparrow Mass",
     xlab = "Mass(g)",
     ylab = "Survial(1)",cex=0.5)
points(mean_mass_survive,1,pch=19,col='green')
points(mean_mass_death,0,pch=19,col='red')


### Various Fits

##cube transform 
tarsus_cube = tarsus_ac^3
tarsus_cube_train = tarsus_ac[train_index]^3
tarsus_cube_test = tarsus_ac[-train_index]^3


fit_cube = glm(fys~tarsus_cube, family=binomial)
prediction_vals_cube = seq(min(tarsus_cube)-1000,max(tarsus_cube)+1000,length.out=1000)
p_cube_probs = predict(fit_cube, newdata = data.frame(tarsus_cube=prediction_vals_cube),type='response')

plot(tarsus_cube,fys_noise,main = "Logistic Regression Fit on First-Year Survival Rate",
     xlab = "Age-adjusted Tarsus Length",
     ylab = "Survial(1)",cex=0.5)
lines(prediction_vals,p_cube_probs,type='l',col='red', lwd=2)

##poly fit
#set degree
degree = 3
#get polynomial terms 
poly_terms = poly(tarsus_ac_train,degree,raw=TRUE)

fit_poly = glm(fys_train~poly_terms, family=binomial(link='logit'))

test_data_poly = data.frame(poly_terms=poly(tarsus_ac_test,degree,raw=TRUE))

predicted_poly_probs = predict(fit_poly,newdata=test_data_poly,type='response')

predictions_poly = ifelse(predicted_poly_probs >=0.5,1,0)
#get error rate
error_rate_poly = mean(predictions_poly != fys_test)







