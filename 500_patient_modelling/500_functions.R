#### Functions ####

# function to select rows for modelling
year1_model_selecter <- function(df, variable) {
  require(tidyverse)
  # create a character vector of abbreviated analyte names for column selection
  analyte_abbrev_list <- read.csv("data/abbreviated_names_lookup_table.csv")
  analyte_abbrev_list <- as.character(analyte_abbrev_list$newValue)
  # remove mononuc from vector
  remove <- "mononuc"
  analyte_abbrev_list <- setdiff(analyte_abbrev_list, remove)

  # remove all date columns. they make the model far to accurate
  date_cols <- colnames(dplyr::select(df, contains("_date")))
  date_cols

  df_cc_no_dates <- df # change here to alter model first blood or all included
  df_cc_no_dates[, date_cols] <- NULL

  # select the columns of interest
  initial_cols <- colnames(dplyr::select(df_cc_no_dates, contains("initial")))
  neuro_cols <- colnames(dplyr::select(df_cc_no_dates, contains("neuro")))
  year1_cols <- colnames(dplyr::select(df_cc_no_dates, contains("year1")))
  asia_grade_status_columns <- colnames(dplyr::select(df_cc_no_dates, contains("by_")))
  discharge_cols <- colnames(dplyr::select(df_cc_no_dates, contains("discharge_")))
  neuro_discharge_cols <- colnames(dplyr::select(df_cc_no_dates, contains("level_discharge")))
  drug_cols <- c("amit", "aspi", "bacl", "gaba", "lact", "omep", "para", "preg", "senn", "tram", "zopi")

  discharge_cols

  day28_cols <- colnames(dplyr::select(df_cc_no_dates, contains("day28_")))
  day28_cols

  # copy from here to change model df and/or model variable
  # remove the entry to be modelled
  year1_motor_cols <- year1_cols[!year1_cols %in% variable]
  year1_motor_cols

  df_2 <- paste(df, variable, sep = "_")

  df_2 <- df_cc_no_dates
  # df_cc_year1_motor_initial <- df_locf_encoded

  # remove date cols
  df_2[, date_cols] <- NULL
  # remove other year1 cols
  df_2[, c(year1_motor_cols, discharge_cols, day28_cols, asia_grade_status_columns, neuro_discharge_cols, drug_cols)] <- NULL
  # remove improver and referal grade cols
  df_2[, c("improver.1", "from_referral", "length_of_stay_spell")] <- NULL
  colnames(df_2)
  # try build model on above, then try by filtering to only the first blood test per patient and see
  return(df_2)
}


## make another version for discharge cols
# function to select rows for modelling
discharge_model_selecter <- function(df, variable) {
  # create a character vector of abbreviated analyte names for column selection
  analyte_abbrev_list <- read.csv("data/abbreviated_names_lookup_table.csv")
  analyte_abbrev_list <- as.character(analyte_abbrev_list$newValue)
  # remove mononuc from vector
  remove <- "mononuc"
  analyte_abbrev_list <- setdiff(analyte_abbrev_list, remove)

  # remove all date columns. they make the model far to accurate
  date_cols <- colnames(dplyr::select(df, contains("_date")))
  date_cols

  df_cc_no_dates <- df # change here to alter model first blood or all included
  df_cc_no_dates[, date_cols] <- NULL

  # select the columns of interest
  initial_cols <- colnames(dplyr::select(df_cc_no_dates, contains("initial")))
  neuro_cols <- colnames(dplyr::select(df_cc_no_dates, contains("neuro")))
  year1_cols <- colnames(dplyr::select(df_cc_no_dates, contains("year1")))
  asia_grade_status_columns <- colnames(dplyr::select(df_cc_no_dates, contains("by_")))
  discharge_cols <- colnames(dplyr::select(df_cc_no_dates, contains("discharge_")))
  neuro_discharge_cols <- colnames(dplyr::select(df_cc_no_dates, contains("level_discharge")))
  drug_cols <- c("amit", "aspi", "bacl", "gaba", "lact", "omep", "para", "preg", "senn", "tram", "zopi")

  discharge_cols

  day28_cols <- colnames(dplyr::select(df_cc_no_dates, contains("day28_")))
  day28_cols

  # copy from here to change model df and/or model variable
  # remove the entry to be modelled
  discharge_cols <- discharge_cols[!discharge_cols %in% variable]


  df_2 <- paste(df, variable, sep = "_")

  df_2 <- df_cc_no_dates
  # df_cc_year1_motor_initial <- df_locf_encoded

  # remove date cols
  df_2[, date_cols] <- NULL
  # remove other year1 cols
  df_2[, c(year1_cols, discharge_cols, day28_cols, asia_grade_status_columns, neuro_discharge_cols, drug_cols)] <- NULL
  # remove improver and referal grade cols
  df_2[, c("improver.1", "from_referral", "length_of_stay_spell")] <- NULL
  colnames(df_2)
  # try build model on above, then try by filtering to only the first blood test per patient and see
  return(df_2)
}



### FUNCTION X ###

# function for perfoming linear regression, glmnet and logistic regression on
# data, with 10-fold cross validation and a 80-20 training-test data split - PCA
# has been added to linreg model with a threshold of 0.9
lm_glmnet_logreg_modelling <- function(predictor_cols, target_cols) {
  require(tidyverse)
  require(caret)

  target_cols_names <- names(target_cols)
  # how many target variables?
  num_var_select <- length(target_cols_names)
  print(num_var_select)

  # create a list of dataframes which contain all predictors
  list(predictor_cols) %>%
    rep(num_var_select) %>%
    enframe(name = "id", value = "predictors") %>%
    mutate(target_name = target_cols_names) -> starter_df
  starter_df

  # add the target variable to the dataframe of predictors (called target)
  for (n in 1:num_var_select) {
    starter_df$predictors[[n]]$target <- target_cols[, n]
  }
  starter_df
  # for some reason target was made a list, below corrects this
  for (n in 1:num_var_select) {
    starter_df$predictors[[n]]$target <- unlist(starter_df$predictors[[n]]$target)
  }



  # create 80% training-test split
  set.seed(436)
  for (n in 1:num_var_select) {
    starter_df$training_samples <- starter_df$predictors[[n]]$target %>%
      createDataPartition(p = 0.8)
  }
  # below gives error, but works correctly
  for (n in 1:num_var_select) {
    starter_df$train_data[[n]] <- starter_df$predictors[[n]][starter_df$training_samples[[n]], ]
    starter_df$test_data[[n]] <- starter_df$predictors[[n]][-starter_df$training_samples[[n]], ]
  }
  starter_df

  # row(starter_df$target_name == "improver.1")

  # extract logistic regression row
  x <- nrow(starter_df)

  df_log_reg <- starter_df[x, ]
  df_log_reg
  # remove row from linear regression sets
  starter_df <- starter_df[-x, ]
  starter_df


  # caret linear regression model with 10-fold cross validation
  linreg_model_function <- function(df) {
    ctrl <- trainControl(
      ## 10-fold CV
      method = "repeatedcv",
      number = 10
      # preProcOptions = list(thresh = 0.9)
    )
    train(
      target ~ .,
      data = df,
      method = "lm",
      # preProcess ="pca",
      trControl = ctrl,
      tuneLength = 10
    )
  }


  # glmnet model for comparison - uses combination of lasso and ridge
  # regression, also has 10-fold cross validation
  glmnet_model_function <- function(df) {
    ctrl <- trainControl(
      ## 10-fold CV
      method = "cv",
      number = 10,
      verboseIter = TRUE
    )
    train(
      target ~ .,
      data = df,
      method = "glmnet",
      trControl = ctrl,
      tuneLength = 10
    )
  }


  knn_function <- function(df) {
    ctrl <- trainControl(
      ## 10-fold CV
      method = "cv",
      number = 10,
      verboseIter = TRUE
    )

    set.seed(4837)

    train(
      target ~ .,
      data = df,
      method = "knn",
      trControl = ctrl,
      tuneLength = 10
    )
  }

  # logistic regression model
  logreg_function <- function(df) {
    ctrl <- trainControl(
      ## 10-fold CV
      method = "repeatedcv",
      number = 10
    )

    train(
      target ~ .,
      data = df,
      method = "glm",
      family = binomial(),
      trControl = ctrl
    )
  }

  # a logistic regression function from base r
  logreg_function <- function(df) {
    glm(target ~ ., family = binomial(link = "logit"), data = df)
  }


  # apply the function to every dataframe in the list - DO NOT DELETE
  df_test_model <- starter_df %>%
    mutate(model = map(.x = train_data, .f = ~ linreg_model_function(.x)))

  # create list of model function
  model_list <- list(
    linreg_model = linreg_model_function,
    glmnet_model = glmnet_model_function
  ) %>%
    enframe(name = "modelName", value = "model")
  # check
  model_list

  # combine the 2 together so train_df contains all combinations of models and targets
  train_df <-
    starter_df[rep(1:nrow(starter_df), nrow(model_list)), ]

  train_df %<>%
    bind_cols(
      model_list[rep(1:nrow(model_list), nrow(starter_df)), ] %>%
        arrange(modelName)
    ) %>%
    mutate(id = 1:nrow(.))
  # check
  train_df

  # apply the model functions to every dataframe in the list -
  # Note: the glmnet models take some time to run (less than 1 minute)
  df_models <- train_df %>%
    mutate(
      params = map(.x = train_data, .f = ~ list(df = .x)),
      model_fits = invoke_map(model, params)
    )

  df_models

  # apply logistic regression fucntion to improver data
  df_log_reg
  df_log_reg %<>%
    mutate(model = map(.x = train_data, .f = ~ logreg_function(df = .x)))
  df_log_reg



  # to get the best tuning paramater of the glmnet model
  # Note: alpha is elastic net mixing paramater - 1 = lasso, 0 = ridge, values between are a mixture
  # lambda is the shrinkage factor - determines amount of penalization
  df_models$model_fits[[1]]$bestTune
  # to see plot of the cross validation
  plot(df_models$model_fits[[1]])

  # number of models has doubled so this must be made again
  num_var_select <- length(df_models$id)


  # predict on test data
  df_models$predictions <- map2(df_models$model_fits, df_models$test_data, ~ predict(.x, .y))

  df_log_reg$predictions <- map2(df_log_reg$model, df_log_reg$test_data, ~ predict(.x, .y))

  # extract model perfomance metrics
  for (n in 1:num_var_select) {
    df_models$RMSE[[n]] <- RMSE(df_models$predictions[[n]], df_models$test_data[[n]]$target)
    df_models$Rsquare[[n]] <- R2(df_models$predictions[[n]], df_models$test_data[[n]]$target)
  }

  df_models

  # to compare models
  compare_models(df_models$model_fits[[1]], df_models$model_fits[[2]])


  best <- list()
  best_result <- list()

  # get best model stats and add RMSESD to df_models for plotting
  # NOTE: The RMSESD obtained here is taken from the internal validation of the
  # models, whereas the RMSE and R^2 cols in df_model are testing the model
  # against the test-split data
  for (n in 1:num_var_select) {
    best[[n]] <- which(rownames(df_models$model_fits[[n]]$results) == rownames(df_models$model_fits[[n]]$bestTune))
    best_result[[n]] <- df_models$model_fits[[n]]$results[best[[n]], ]
    df_models$RMSESD[[n]] <- best_result[[n]]$RMSESD
  }
  # create list of objects to be returned by function
  modelling_function_list <- list(
    "df_models" = df_models,
    "df_log_reg" = df_log_reg
  )

  return(modelling_function_list)
}


### FUNCTION x ###

# Function to oversample, undersample and do a mixture of both on data and then
# apply a logistic regrssion model to them
sampleing_function <- function(df_log_reg) {
  # oversample data
  require(ROSE)
  # get the number of improvers in the data
  num_improvers <- sum(df_log_reg$train_data[[1]]$target)
  num_observations <- nrow(df_log_reg$train_data[[1]])
  num_nonimprovers <- num_observations - num_improvers

  table(df_log_reg$train_data[[1]]$target)

  df_log_reg_oversample <- list()
  df_log_reg_undersample <- list()
  df_log_reg_both <- list()
  # df_log_reg[2:11,] <- 0
  for (n in 1:10) {
    df_log_reg_oversample[[n]] <- ovun.sample(target ~ .,
      data = df_log_reg$train_data[[1]],
      method = "over", N = (num_nonimprovers * 2)
    )$data
    df_log_reg_undersample[[n]] <- ovun.sample(target ~ .,
      data = df_log_reg$train_data[[1]],
      method = "under", N = (num_improvers * 2)
    )$data
    df_log_reg_both[[n]] <- ovun.sample(target ~ .,
      data = df_log_reg$train_data[[1]],
      method = "both", N = num_observations
    )$data
  }

  # a logistic regression function from base r
  logreg_function <- function(df) {
    glm(target ~ ., family = binomial(link = "logit"), data = df)
  }

  # make list of balanced dataframes
  df_log_reg_balanced <- list(df_log_reg_oversample, df_log_reg_undersample, df_log_reg_both)
  # apply logistic regression function to balanced data
  df_log_reg_balanced_models <- pmap(df_log_reg_balanced, .f = ~ logreg_function(df = .x)) # output should be list of lists with 3 elements - not working

  oversample_model <- map(df_log_reg_oversample, .f = ~ logreg_function(df = .x))
  oversample_model <- set_names(oversample_model, rep("oversample", 10))

  undersample_model <- map(df_log_reg_undersample, .f = ~ logreg_function(df = .x))
  undersample_model <- set_names(undersample_model, rep("undersample", 10))

  both_model <- map(df_log_reg_both, .f = ~ logreg_function(df = .x))
  both_model <- set_names(both_model, rep("both", 10))

  # predict balanced dataframes on test data
  oversample_model_predictions <- map2(oversample_model, df_log_reg$test_data, ~ predict(.x, .y, type = "response"))
  undersample_model_predictions <- map2(undersample_model, df_log_reg$test_data, ~ predict(.x, .y, type = "response"))
  both_model_predictions <- map2(both_model, df_log_reg$test_data, ~ predict(.x, .y, type = "response"))

  knn_function <- function(df) {
    ctrl <- trainControl(
      ## 10-fold CV
      method = "cv",
      number = 10,
      verboseIter = TRUE
    )

    set.seed(4837)

    train(
      target ~ .,
      data = df,
      method = "knn",
      trControl = ctrl,
      tuneLength = 10
    )
  }

  # df_knn_balanced_models <- pmap(df_log_reg_balanced, .f = ~ knn_function(df = .x)) # output should be list of lists with 3 elements - not working

  oversample_knn_model <- map(df_log_reg_oversample, .f = ~ knn_function(df = .x))
  oversample_knn_model <- set_names(oversample_knn_model, rep("oversample", 10))

  undersample_knn_model <- map(df_log_reg_undersample, .f = ~ knn_function(df = .x))
  undersample_knn_model <- set_names(undersample_knn_model, rep("undersample", 10))

  both_knn_model <- map(df_log_reg_both, .f = ~ knn_function(df = .x))
  both_knn_model <- set_names(both_knn_model, rep("both", 10))

  # predict balanced dataframes on test data
  oversample_knn_model_predictions <- map2(oversample_knn_model, df_log_reg$test_data, ~ predict(.x, .y))
  undersample_knn_model_predictions <- map2(undersample_knn_model, df_log_reg$test_data, ~ predict(.x, .y))
  both_knn_model_predictions <- map2(both_knn_model, df_log_reg$test_data, ~ predict(.x, .y))


  # put into list of lists
  df_log_reg_balanced_model_predictions <- list(
    oversample_model_predictions, undersample_model_predictions,
    both_model_predictions, oversample_knn_model_predictions,
    undersample_knn_model_predictions, both_knn_model_predictions
  )



  output_list <- list(
    "oversample_model_predictions" = oversample_model_predictions,
    "undersample_model_predictions" = undersample_model_predictions,
    "both_model_predictions" = both_model_predictions,
    "oversample_knn_model_predictions" = oversample_knn_model_predictions,
    "undersample_knn_model_predictions" = undersample_knn_model_predictions,
    "both_knn_model_predictions" = both_knn_model_predictions,
    "df_log_reg_balanced_model_predictions" = df_log_reg_balanced_models,
    "df_log_reg_oversample" = df_log_reg_oversample,
    "df_log_reg_undersample" = df_log_reg_undersample,
    "df_log_reg_both" = df_log_reg_both,
    "oversample_knn_model" = oversample_knn_model,
    "undersample_knn_model" = undersample_knn_model,
    "both_knn_model" = both_knn_model
  )

  return(output_list)
}


### FUNCTION X ###

# Function to extract variable importance
var_importance_extraction <- function(df_models) {

  # filter out knn - it breaks the following code
  df_models <- df_models %>%
    dplyr::filter(modelName != "knn_model")
  # get varaible importance for all models
  variable_importance <- lapply(df_models$model_fits, varImp)
  # give the list relevent names
  names(variable_importance) <- paste(df_models$target_name, df_models$modelName, sep = "_")

  # filter var importance dfs to only have important variables (should be the same as selecting all significant variables?)
  var_imp_df_filtered <- map(variable_importance, ~ rownames_to_column(.x$importance, "variable") %>%
    dplyr::filter(Overall > 0))

  # create list of names to apply to plots
  var_imp_names <- map2(df_models$target_name, df_models$modelName, ~ paste("Variable Importance:",
    print(.x),
    print(.y),
    sep = " "
  ))

  # show plots of variable importance for all models
  var_imp_plots <- map2(variable_importance, var_imp_names, ~ plot(.x, main = .y))
  set_names(var_imp_plots, var_imp_names)

  result_list <- list(var_imp_df_filtered, var_imp_plots)

  return(result_list)
}

### FUNCTION x ###

# Function to create list of coefficients filtered to only include significant
# variables from linear regression models

linreg_significant_coeff <- function(df_models) {
  # filter to linreg models only
  lm_models_only <- dplyr::filter(df_models, modelName == "linreg_model")
  # make list of model summaries
  lm_models_summary <- map(lm_models_only$model_fits, ~ summary(.x))
  # give list names
  lm_models_summary <- set_names(lm_models_summary, lm_models_only$target_name)

  # make list of model coefficients
  lm_models_coefs <- map(lm_models_summary, ~ coef(.x))
  # convert list of matrix to dataframes - note the optional argument keeps the rownames
  lm_models_coefs <- map(lm_models_coefs, ~ as.data.frame(.x, optional = TRUE))

  # filter coeffs to only keep significant variables
  lm_models_coefs_sig <- map(lm_models_coefs, ~ rownames_to_column(.x, "variable") %>%
    dplyr::filter(`Pr(>|t|)` <= 0.05))
  return(lm_models_coefs_sig)
}


### FUNCTION x ###

# Function to extract all coefficient above 0 from glmnet models

glm_coef_extractor <- function(df_models) {

  # filter to only glmnet models
  glm_models <- df_models %>%
    dplyr::filter(modelName == "glmnet_model")
  # get all model variables used - NOTE: if this isn't working make sure to
  # install all the packages in 500_modelling.R
  glm_coefs <- map2(glm_models$model_fits, glm_models$model_fits, ~ coef(.x$finalModel, .y$bestTune$lambda))
  # give list elements names
  glm_coefs <- set_names(glm_coefs, glm_models$target_name)

  # convert to matrix
  a <- map(glm_coefs, ~ as.matrix(.x))
  # convert to dataframe
  b <- map(a, ~ as.data.frame(.x, optional = TRUE))
  # convert rowname to column
  c <- map(b, ~ rownames_to_column(.x, "variable"))
  # rename second column form 1 to coef
  d <- map(c, ~ dplyr::rename(.x, coef = 2))
  # filter to only have coefs greater than 0
  glm_coefs_filtered <- map(d, ~ dplyr::filter(.x, coef > 0))

  return(glm_coefs_filtered)
}
