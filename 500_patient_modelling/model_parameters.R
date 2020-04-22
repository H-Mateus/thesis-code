# This is the parameters of the models used

library(caret)

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
