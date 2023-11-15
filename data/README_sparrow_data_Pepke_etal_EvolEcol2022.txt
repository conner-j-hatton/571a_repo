README file for individual sparrow data for 
Pepke et al. (2022) Causes and consequences of variation in early-life telomere length in a bird metapopulation, Ecology and Evolution.

Below is a short explanation of variables in each coloumn. For details see the Methods section in the main text.
NA indicates missing data.

ID : individual identity number
year :	year of hatching
date :	date of telomere sampling
TL : relative telomere length
pop_size_meancentered :	spring pre-breeding census population size in the hatch year mean centered within populations
island_name : name of natal island
island_ID : identity number for natal island
brood_ID : brood identity number
sampling_age_days : age in days at telomere sampling	
clutch_size : number of chicks in the nest at telomere sampling
sex : male = 1 and female = 2
tarsus : tarsometatarsus (tarsus) length measured in mm
mass :	body mass measured in grams
tarsus_agecorrected : age-standardized tarsus length (i.e. residuals of a tarsus vs. age model, see Methods)
condition : body condition (residuals of a linear regression of log10-transformed mass against log10-transformed tarsus length)	
hatch_date : date of hatching	
hatchday_meancentered :	date of hatching (mean centered ordinal day of the year)
firstyear_survival : indicates whether an individual survived first year (0 = no, 1 = yes)
min_lifespan_days : minimum observed lifespan measured in days	
censored : indicates whether an individual is censored (i.e. still alive at the end of the study) in survival analyses (censored = 1, not censored = 2)
dispersal : indicates whether an individual dispersed to another island during its first year (= 1, natal dispersal) or not (= 0)