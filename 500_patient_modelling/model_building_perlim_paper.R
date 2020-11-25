#load data
factors_spss <- read.csv("10_month_report/data/01_2019-07-08_spss_data_with_factors.csv")

#load packages
library(caret)
library(leaps)
library(Hmisc)
library(naniar)
library(tidyverse)
library(ggplot2)
library(vtreat) 

#convert numerics to factors
factors_spss[,c(1:12, 14:20)] <- lapply(factors_spss[,c(1:12, 14:20)], as.factor)
factors_spss[,c(13, 27:30)] <- lapply(factors_spss[,c(13, 27:30)], as.numeric)
#impute missing values
#factors_spss[,30] <- lapply(factors_spss[,30], impute)
#remove bloods and intial neurology for PCA-I

#summary stats for initial neruology
perlim_inital_neurology <- describe(factors_spss[,21:25])
#remove unwanted columns - not working
perlim_inital_neurology <- as.data.frame(perlim_inital_neurology)
#perlim_inital_neurology[,c("1", "2", "6:7", "11:12")] <- NULL
save(perlim_inital_neurology, file = "10_month_report/data/perlim_initial_neruology_summary.rda")

pcai <- factors_spss
pcai[,c(1, 15, 16, 18, 20:26, 32, 38:74, 81:85)] <- NULL

#hot one encoding of variables using carets dummyVars function
test <- dummyVars(" ~ .", data = pcai, fullRank=TRUE)
pcai <- data.frame(predict(test, newdata = pcai))

#make 3-month outcome data to predict
pcai_3month <- pcai
pcai_3month[,c(27:31)] <- NULL
#make dataframes for each individual outcome measure to predict
pcai_3month_individual <- pcai_3month
pcai_3month_individual[,c(22, 26)] <- NULL 
pcai_3month_individual_motor <- pcai_3month_individual
pcai_3month_individual_motor[,c(22:23)] <- NULL
pcai_3month_individual_touch <- pcai_3month_individual
pcai_3month_individual_touch[,c(23:24)] <- NULL
pcai_3month_individual_pain <- pcai_3month_individual
pcai_3month_individual_pain[,c(22, 24)] <- NULL
#pcai 12-month dataframes 
pcai_12month <- pcai
pcai_12month[,c(22:26, 27, 31)] <- NULL
pcai_12month_motor <- pcai_12month
pcai_12month_motor[,c(22:23)] <- NULL
pcai_12month_touch <- pcai_12month
pcai_12month_touch[,c(23, 24)] <- NULL
pcai_12month_pain <- pcai_12month
pcai_12month_pain[,c(22, 24)] <- NULL
#create dataframes for pcac
pcac <- factors_spss
pcac[,c(1, 15, 16, 18, 20:26, 32, 38:80)] <- NULL
#hot one encoding of variables using carets dummyVars function
test <- dummyVars(" ~ .", data = pcac, fullRank=TRUE)
pcac <- data.frame(predict(test, newdata = pcac))

pcac_3month <- pcac
pcac_3month[,c(27:31, 23:25)] <- NULL
pcac_3month_motor <- pcac_3month
pcac_3month_motor[,c(23)] <- NULL
pcac_3month_sensory <- pcac_3month
pcac_3month_sensory[,c(22)] <- NULL
#pcac 12-month
pcac_12month <- pcac
pcac_12month[,c(22:26, 28:30)] <- NULL
pcac_12month_motor <- pcac_12month
pcac_12month_motor[,c(23)] <- NULL
pcac_12month_sensory <- pcac_12month 
pcac_12month_sensory[,c(22)] <- NULL

pcai_nas <- miss_var_summary(pcai_3month_individual_motor)



#make regression models of outcomes - note: using na.action = na.pass gives better R^2 but 
  #it only works if theres no missing outcome data, otherwise na.omit has to be used
set.seed(56)
# Fit lm model: 
model_3month_pcai_motor <- train(
  X.3mASIAMotor ~ ., 
  pcai_3month_individual_motor,
  method = "leapForward",
  trControl = trainControl(
    method = "cv", 
    number = 10,
    verboseIter = TRUE,
      ),
  preProcess = "medianImpute",
  na.action = na.omit
)

model_3month_pcai_touch <- train(
  X.3mTouch ~ ., 
  pcai_3month_individual_touch,
  method = "leapForward",
  trControl = trainControl(
    method = "cv", 
    number = 10,
    verboseIter = TRUE
  ),
  preProcess = "medianImpute",
  na.action = na.omit
)

model_3month_pcai_pain <- train(
  X.3mPain ~ ., 
  pcai_3month_individual_pain,
  method = "leapForward",
  trControl = trainControl(
    method = "cv", 
    number = 10,
    verboseIter = TRUE
  ),
  preProcess = "medianImpute",
  na.action = na.omit
)

model_12month_pcai_pain <- train(
  X.1yPain ~ ., 
  pcai_12month_pain,
  method = "leapForward",
  trControl = trainControl(
    method = "cv", 
    number = 10,
    verboseIter = TRUE
  ),
  preProcess = "medianImpute",
  na.action = na.omit
)

model_12month_pcai_touch <- train(
  X.1yTouch ~ ., 
  pcai_12month_touch,
  method = "leapForward",
  trControl = trainControl(
    method = "cv", 
    number = 10,
    verboseIter = TRUE
  ),
  preProcess = "medianImpute",
  na.action = na.omit
)

model_12month_pcai_motor <- train(
  X.1yASIAMotor ~ ., 
  pcai_12month_motor,
  method = "leapForward",
  trControl = trainControl(
    method = "cv", 
    number = 10,
    verboseIter = TRUE
  ),
  preProcess = "medianImpute",
  na.action = na.omit
)

model_3month_pcac_motor <- train(
  X.3mMotor ~ ., 
  pcac_3month_motor,
  method = "leapForward",
  trControl = trainControl(
    method = "cv", 
    number = 10,
    verboseIter = TRUE
  ),
  preProcess = "medianImpute",
  na.action = na.omit
)

model_3month_pcac_sensor <- train(
  X.3mASIASensory ~ ., 
  pcac_3month_sensory,
  method = "leapForward",
  trControl = trainControl(
    method = "cv", 
    number = 10,
    verboseIter = TRUE
  ),
  preProcess = "medianImpute",
  na.action = na.omit
)

model_12month_pcac_motor <- train(
  X.1yMotor ~ ., 
  pcac_12month_motor,
  method = "leapForward",
  trControl = trainControl(
    method = "cv", 
    number = 10,
    verboseIter = TRUE
  ),
  preProcess = "medianImpute",
  na.action = na.omit
)

model_12month_pcac_sensor <- train(
  X.1yASIASensory ~ ., 
  pcac_12month_sensory,
  method = "leapForward",
  trControl = trainControl(
    method = "cv", 
    number = 10,
    verboseIter = TRUE
  ),
  preProcess = "medianImpute",
  na.action = na.omit
)

#"lm" method test
test <- train(
  X.1yASIASensory ~ ., 
  pcac_12month_sensory,
  method = "lm",
  trControl = trainControl(
    method = "cv", 
    number = 10,
    verboseIter = TRUE
  ),
  preProcess = "medianImpute",
  na.action = na.omit
)

#below is how to get coefficients from an "lm" model, but this doesn't work with the leapForward models
coef(test$finalModel)

#use summary to get what variables were used in the model 
summary(model_12month_pcac_motor) #factor 5 (initial neurology) and CCS.1 are significant?
summary(model_12month_pcac_sensor) #factor 5 and vertabral fracture 
summary(model_3month_pcac_motor) #factor 5 and age
summary(model_3month_pcac_motor) #factor 5 and age 
summary(model_12month_pcai_motor) # factor 3, level, age, smoker
summary(model_12month_pcai_pain) # factor 3 and vertabral fracture 
summary(model_12month_pcai_touch) # factor 3 and vertebral fracture
summary(model_3month_pcai_touch) #factor 3(initial neurology) and level are significant?
summary(model_3month_pcai_pain)  #facotor 3, 4(inflam and liver) and 5(liver) + level significant 
summary(model_3month_pcai_motor) #factor 3, level and CCS significant

#to get the coefficients - note: the second part of the call (after the comma) is the number of variables
  #you want the coeffiecients of
coef(model_12month_pcac_motor$finalModel, as.numeric(model_12month_pcac_motor$bestTune))
coef(model_12month_pcac_sensor$finalModel, as.numeric(model_12month_pcac_sensor$bestTune))
coef(model_12month_pcai_motor$finalModel, as.numeric(model_12month_pcai_motor$bestTune))
coef(model_12month_pcai_pain$finalModel, as.numeric(model_12month_pcai_pain$bestTune))
coef(model_12month_pcai_touch$finalModel, as.numeric(model_12month_pcai_touch$bestTune))
coef(model_3month_pcac_motor$finalModel, as.numeric(model_3month_pcac_motor$bestTune))
coef(model_3month_pcac_sensor$finalModel, as.numeric(model_3month_pcac_sensor$bestTune))
coef(model_3month_pcai_motor$finalModel, as.numeric(model_3month_pcai_motor$bestTune))
coef(model_3month_pcai_pain$finalModel, as.numeric(model_3month_pcai_pain$bestTune))
coef(model_3month_pcai_touch$finalModel, as.numeric(model_3month_pcai_touch$bestTune))
#the coefs seem fairly similar to table 4 of the paper, but there are differences, and often there is one
  #less variable in the model here then there is in table 4 - maybe due to something about the selection
  #of "bestTune" is different to what SPSS does? 

#add id row to df 
df_individual_2 <- tibble::rowid_to_column(df_individual, "id")

