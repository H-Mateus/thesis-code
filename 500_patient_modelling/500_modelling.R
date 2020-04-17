
#load packages
library(caret)
library(leaps)
library(Hmisc)
library(naniar)
library(tidyverse)
library(ggplot2)
library(vtreat) 
library(purrr)
library(Metrics)
library(magrittr)
library(mlbench)
library(reshape2)
library(psych) # for factor analysis
library(ROCR)
library(Metrics)
library(ROSE) # for balancing data (oversampling)
library(InformationValue) # for logistic regression analysis
library(umap)

#load data - Note: use of .rda here is to maintain df column data types
load('data/df_locf_first_blooc_complete_2019-10-21.rda')
load('data/df_complete_first_blood_2019-10-21.rda')

# create version of data that contains improver col for logistic regression 
#df_locf_first_cc_logistic <- df_locf_first_blood_cc
# remove improver col from data for linear regression 
#df_locf_first_blood_cc$improver.1 <- NULL
#df_cc_first_blood$improver.1 <- NULL

# create vector of model targets
target_cols <- dplyr::select(df_locf_first_blood_cc, contains("_motor"), 
                             contains("_sensor"), contains("_scim"))
target_cols <- dplyr::select(target_cols, contains("discharge"), contains("year1"))
# add improver col
target_cols <- cbind(target_cols, df_locf_first_blood_cc$improver.1)
target_cols <- target_cols %>% 
  dplyr::rename(improver.1 = `df_locf_first_blood_cc$improver.1`)
target_cols_names <- names(target_cols)
head(target_cols)
target_cols_names


# create vector of predictor cols 
predictor_cols <- df_locf_first_blood_cc[, !(names(df_locf_first_blood_cc) %in% target_cols_names)]
predictor_cols_names <- names(predictor_cols)
head(predictor_cols)
predictor_cols_names

# select cols for logistic regression 
#df_log_reg <- (df_locf_first_cc_logistic[, c(predictor_cols_names, "improver.1")])
#names(df_log_reg)

# check the dimensions
dim(df_locf_first_blood_cc)
dim(target_cols)
dim(predictor_cols)

# how many target variables?
num_var_select <- length(target_cols_names)
print(num_var_select)

# every patient is neuro intact so column adds no value
predictor_cols$neuro_intact.Yes <- NULL
predictor_cols$discharge <- NULL

# create a character vector of abbreviated analyte names for column selection 
analyte_abbrev_list <- read.csv('data/abbreviated_names_lookup_table.csv')
analyte_abbrev_list <- as.character(analyte_abbrev_list$newValue)
# remove mononuc from vector
remove <- "mononuc"
analyte_abbrev_list <- setdiff(analyte_abbrev_list, remove)
load('data/filtered_bloods_2019-10-29.rda')
bloods_to_remove <- bloods_to_remove_df$variable
remaining_bloods <- setdiff(analyte_abbrev_list, bloods_to_remove)
print(remaining_bloods)

# create vector of initial neurology columns 
initial_cols <- names(dplyr::select(df_locf_first_blood_cc, contains("initial")))

                  # plot histogram of every column to check for errors 
    #hist.data.frame(predictor_cols) # this line halts the script as it requires a click to see second page of graphs




    ### MODELLING ###

# use modelling function - note pca with a threhold of 0.9 has been added to linreg model
source('scripts/500_functions.R')
# NOTE: this next line takes around half a minute to run - the models are stored
# in a tibble, the names of the lists within the tibble are the model outcome
modelling_function_list <- lm_glmnet_logreg_modelling(predictor_cols, target_cols)
# assign model function outputs 
df_models <- modelling_function_list[[1]]
df_log_reg <- modelling_function_list[[2]]

# save tibbles 
#write_rds(df_models, path = 'data/df_models_500_patients_2019-11-07.rds')
df_model_predictors <- df_models$predictors[[1]]
#save(df_model_predictors, file = 'data/df_models_500_patient_predictors_2019-11-07.rda')
#saveRDS(df_model_predictors, file = 'data/df_models_500_patient_predictors_2019-11-07.rds')
#write_rds(df_log_reg, path = 'data/df_log_reg_500_patients_2019-11-12.rds')

# check logistic regression models
# pull balanced dataset predictions 
# get the number of improvers in the data 
num_improvers <- sum(df_log_reg$train_data[[1]]$target)
num_observations <- nrow(df_log_reg$train_data[[1]])
num_nonimprovers <- num_observations - num_improvers
#NOTE: next line takes a few seconds to run
test <- sampleing_function(df_log_reg)
names(test)


oversample_model_predictions <- test$oversample_model_predictions
undersample_model_predictions <- test$undersample_model_predictions
both_model_predictions <- test$both_model_predictions

names(oversample_model_predictions)

# pull out knn models 
oversample_knn_model_predictions <- test$oversample_knn_model_predictions
undersample_knn_model_predictions <- test$undersample_knn_model_predictions
both_knn_model_predictions <- test$both_knn_model_predictions
oversample_knn_model_predictions[[1]]

#number of models has doubled so this must be made again
num_var_select <- length(df_models$id)
# to get all the summaries:			
for (n in 1:num_var_select){
  print(df_models$target_name[[n]])
  print(summary(df_models$model_fits[[n]]))
}


# dotplot of R^2 - note: not sure if this is adjusted R^2 or not
lattice::dotplot(Rsquare~target_name|modelName, df_models)

# plot predicted values against true values? 
pred <- map2(df_models$predictions, df_models$test_data, ~ data.frame(prediction = .x, target = .y$target))
map(pred, ~ plot(.x))

glm_interactions <- glm(target ~ .*., data = df_models$train_data[[1]])
glm_interactions
# anova(glm_interactions) # this takes a while to run - also doesn't seem to be working

lm_models <- dplyr::filter(df_models, modelName == "linreg_model")
lm_interactions <- map(lm_models$train_data, ~ lm(target ~ .*., data = .x))
map(lm_interactions, ~ anova(.x)) # no f values or p values for interactions? 
# below is anova of caret final model
map(lm_models$model_fits, ~ anova(.x$finalModel))

# plot RMSE of all models
# NOTE: the use of RMSESD as error bars may not be sensible
df_models %>% 
  ggplot(aes(x=target_name,color=modelName))+
  geom_point(aes(y=RMSE),size=2) +
  geom_errorbar(aes(ymin = RMSE-RMSESD,ymax= RMSE+RMSESD),size=.5,width=.15) 
  #ylim(c(0, 120))


# to extract coefficients of only significant variables use:
data.frame(summary(df_models$model_fits[[10]])$coef[summary(df_models$model_fits[[10]])$coef[,4] <= .05, 4])

# following function extracts coefficients of significant variables (p <= 0.05) from lm models 
lm_models_coefs_sig <- linreg_significant_coeff(df_models)
lm_models_coefs_sig


# extract variable importance 

# below function creates 2 lists, 1st a list of dataframes that filters out any
# variables that have importance of 0, and 2nd a list of plot of variable
# importance
var_imp_list <- var_importance_extraction(df_models)
var_imp_df_filtered <- var_imp_list[[1]]
var_imp_plots <- var_imp_list[[2]]

# show all dataframes of var importance
var_imp_df_filtered

# show all plots 
var_imp_plots







hist(df_models$train_data[[2]]$target)
temp <- df_models$train_data[[2]]$target
var_imp_plots[[4]]


# get model coefficeints 
coef(df_models$model_fits[[4]]$finalModel)
# are the coefficients with values all the significant variables?
coef(df_models$model_fits[[1]]$finalModel, df_models$model_fits[[1]]$bestTune$lambda)

# get all glmnet model coefficients - function filters to coef > 0
glm_coefs_filtered <- glm_coef_extractor(df_models)
glm_coefs_filtered


resamples(df_models$model_fits) %>% summary(metric = "RMSE") # why does this only have NAs?

# correlation between preditions and actual values and R^2 - Note: this is just
# a sanity check, the R^2 this code generates is the same as the Rsqaure column
# in df_models
correlations <- map2(df_models$test_data, df_models$predictions, ~ cor(.x$target, .y))
correlations <- set_names(correlations, df_models$target_name)
correlations
rsq <- map2(df_models$test_data, df_models$predictions, ~ cor(.x$target, .y)^2)
rsq <- set_names(rsq, df_models$target_name)
rsq















# get results of each model 
results <- map(df_models$model_fits, ~ print(.x$results))

glm_models <- dplyr::filter(df_models, modelName == "glmnet_model")
glm_final_model_plots <- map(glm_models$model_fits, ~ plot(.x$finalModel, xvar = "lambda", label = T))

map(glm_models$model_fits, ~ plot(.x$finalModel, xvar = "dev", label = T))
glm_crossval_rmse_plots <- map(glm_models$model_fits, ~ plot(.x))
names(glm_crossval_rmse_plots) <- glm_models$target_name
glm_crossval_rmse_plots


# note the code below is how you extract elements of a list to the global
# environment - the object names are set to the list element names
#list2env(glm_coefs_filtered ,.GlobalEnv)

# extract varaible column 
variables <- list()
for(n in 1:length(glm_coefs_filtered)) {
  variables[[n]] <- glm_coefs_filtered[[n]]$variable
}
variables
# remove intercept
for(n in 1:length(variables)) {
  variables[[n]] <- variables[[n]][-1]
}
variables

# convert models to long format
glm_models <- dplyr::filter(df_models, modelName == "glmnet_model")
glm_models_long <- map(glm_models$predictors, ~ melt(.x))
glm_models_long <- set_names(glm_models_long, glm_models$target_name)

all_vars <- as.character(unique(glm_models_long[[1]]$variable))
# create boxplots based on all variables
# function for generating boxplots
gg_boxplot <- function(df, variable_to_plot) {
  ggplot(dplyr::filter(df, variable == variable_to_plot),
         aes(x = factor(df_log_reg$predictors[[1]]$target), y = value)) +
    geom_boxplot() +
    ylab(variable_to_plot) +
    xlab("Improver")
}

# example of function use
gg_boxplot(glm_models_long[[1]], "haemocrit")

var_list <- map(variables, ~ as.list(.x))

# pmap(var_list, ~ map2(.x, glm_models_long ~ gg_boxplot(.y, .x)))

discharge_motor_glm_coef_list <- var_list[[1]]
map(discharge_motor_glm_coef_list, ~ gg_boxplot(glm_models_long[[1]], .x))



ggplot(dplyr::filter(glm_models_long[[2]], variable == variables[[2]][2]),
       aes(x = factor(df_log_reg$predictors[[1]]$target), y = value)) +
  geom_boxplot()

plots <- list()
for (i in 1:length(variables)) {
  plots[[i]] <-  ggplot(dplyr::filter(glm_models_long[[i]], variable == variables[[i]][i]),
                        aes(x = factor(df_log_reg$predictors[[1]]$target), y = value)) +
    geom_boxplot()
}


for (i in 1:length(variables)) {
  plots[[i]] <-  ggplot(dplyr::filter(glm_models_long[[i]], variable == all_vars[i]),
                        aes(x = factor(df_log_reg$predictors[[1]]$target), y = value)) +
    geom_boxplot() +
    ggtitle(all_vars[i])
}



variables[[1]][1]
df_models$predictors[[1]][,c(variables[[1]][1])]
as.list(variables[[1]])
a <- map(variables, ~ as.list(.x))

class(variables[[1]])
map2(glm_models$predictors, a[[1]], ~ ggplot(.x, 
                                             aes(x = factor(df_log_reg$predictors[[1]]$target), 
                                                 y = .x[,c(.y)])) +
       geom_boxplot())

