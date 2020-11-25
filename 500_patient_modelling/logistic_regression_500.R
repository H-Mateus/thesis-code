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

# load data 
df_log_reg <- read_rds('../data/df_log_reg_500_patients_2019-11-12.rds')


##### Below is code from first_paper_modelling.R #####

# use modelling function
source('scripts/first_paper_modelling.R')
modelling_function_list <- lm_glmnet_logreg_modelling(predictor_cols_first_paper, target_cols_first_paper)
# assign model function outputs 
df_models_first_paper <- modelling_function_list[[1]]
df_log_reg_first_paper <- modelling_function_list[[2]]

df_models
df_models_first_paper
names(df_models$predictors[[1]])
names(df_models_first_paper$predictors[[1]])

dim(df_models_first_paper$test_data[[1]])
dim(df_models_first_paper$train_data[[1]])


# pull balanced dataset predictions - function gives error, but run code
# manually and it works - need to figure out why it doesn't just work
temp <- df_log_reg # save df
df_log_reg <- df_log_reg_first_paper # rename to run function code manually

# get the number of improvers in the data - if these three lines aren't run the function below doesn't work
num_improvers <- sum(df_log_reg$train_data[[1]]$target)
num_observations <- nrow(df_log_reg$train_data[[1]])
num_nonimprovers <- num_observations - num_improvers

balanced_first_paper <- sampleing_function(df_log_reg = df_log_reg_first_paper)
df_log_reg <- temp # restore df 
names(balanced_first_paper)


oversample_model_predictions_first_paper <- balanced_first_paper$oversample_model_predictions
undersample_model_predictions_first_paper <- balanced_first_paper$undersample_model_predictions
both_model_predictions_first_paper <- balanced_first_paper$both_model_predictions

names(oversample_model_predictions_first_paper)

# pull out knn models 
oversample_knn_model_predictions_first_paper <- balanced_first_paper$oversample_knn_model_predictions
undersample_knn_model_predictions_first_paper <- balanced_first_paper$undersample_knn_model_predictions
both_knn_model_predictions_first_paper <- balanced_first_paper$both_knn_model_predictions
oversample_knn_model_predictions_first_paper[[1]]


# unsupervised learning test
knm_test <- df_log_reg$predictors[[1]]
knm_test <- knm_test[,c(remaining_bloods, "target")]

knm_df <- knm_test
knm_df$target <- NULL
kmean_result <- kmeans(knm_df, 2)

knm_test$prediction <- kmean_result$cluster

#table(knm_test$target, knm_test$prediction)
knm_df_pca <- prcomp(knm_df)
summary(knm_df_pca)


knm_df_pca_plot <- data.frame(first = knm_df_pca$x[,1], second = knm_df_pca$x[,2], improver = knm_test$target,
                              prediction = knm_test$prediction)


knm_df_pca_plot <- data.frame(first = knm_df_pca$x[,1], second = knm_df_pca$x[,2],
                              prediction = knm_test$prediction)


ggplot(knm_df_pca_plot, aes(x = first, y = second, shape = as.factor(prediction))) +
  geom_point() +
  theme_bw() +
  ggtitle('Principle Component Analysis') +
  scale_x_continuous('First Principle Component') + scale_y_continuous('Second Principle Component')

pca1_loadings <- knm_df_pca$rotation[,1]
pcal_loadings2 <- round(pca1_loadings, 3)


library(factoextra)
data("multishapes")
df <- multishapes[, 1:2]
df <- knm_df


set.seed(123)
km.res <- kmeans(df, 5, nstart = 25)
fviz_cluster(km.res, df,  geom = "point", 
             ellipse= FALSE, show.clust.cent = FALSE,
             palette = "jco", ggtheme = theme_classic())

library(fpc)
library(dbscan)
# Load the data 
data("multishapes", package = "factoextra")
#df <- multishapes[, 1:2]

# Compute DBSCAN using fpc package
library("fpc")
set.seed(123)

dbscan::kNNdistplot(df, k = 2)
abline(h = 0.15, lty = 2)

db <- fpc::dbscan(df, eps = 30, MinPts = 5)



# Plot DBSCAN results
library("factoextra")
fviz_cluster(db, data = df, stand = FALSE,
             ellipse = FALSE, show.clust.cent = FALSE,
             geom = c("point", "text"),palette = "jco", ggtheme = theme_classic())


var_list <- map(variables, ~ as.list(.x))

# pmap(var_list, ~ map2(.x, glm_models_long ~ gg_boxplot(.y, .x)))

discharge_motor_glm_coef_list <- var_list[[1]]
map(discharge_motor_glm_coef_list, ~ gg_boxplot(glm_models_long[[1]], .x))



wilcox_test <- function(df, variable_to_test, target) {
  a <- dplyr::filter(df, variable == variable_to_test) %>%
    dplyr::mutate(improver = (target)) 
  
  b <- wilcox.test(a$value ~ a$improver)
  return(b)
}

wilcox_test(glm_models_long[[1]], "mean_cell_hb", df_log_reg$predictors[[1]]$target)

map(discharge_motor_glm_coef_list, ~ wilcox_test(.x))

wilcox.test(a$value ~ a$improver)

nrow(df_log_reg$predictors[[1]])




################## logistic regression analysis ######################

# to get summaries for logreg model 
summary(df_log_reg$model[[1]])
summary(df_log_reg_first_paper$model[[1]])

anova(df_log_reg$model[[1]], test = "Chisq")
anova(df_log_reg_first_paper$model[[1]], test = "Chisq")

anova(df_log_reg$model[[1]]) # model feature interaction? 


# predict on test data - not working 
df_log_reg_first_paper_model_predictions <- predict(df_log_reg_first_paper$model[[1]], df_log_reg_first_paper$test_data[[1]], type = "response")
df_log_reg_model_predictions <- predict(df_log_reg$model[[1]], newdata = df_log_reg$test_data[[1]], type = "response")

#confint(df_log_reg_model_predictions)

#set probability threshold to 0.5
log_predict <- ifelse(df_log_reg_model_predictions >0.5, 1, 0)
log_predict_first_paper <- ifelse(df_log_reg_first_paper_model_predictions >0.5, 1, 0)

library(InformationValue) # to find optiaml cut off with log_reg model 
opt_cut_off <- optimalCutoff(df_log_reg$test_data[[1]]$target, df_log_reg_model_predictions)
opt_cut_off
# optimal cut off is found to be 0.28 - this seems very low
opt_cut_off_first_paper <- optimalCutoff(df_log_reg_first_paper$test_data[[1]]$target, df_log_reg_first_paper_model_predictions)
opt_cut_off_first_paper # cut of is 0.01...

# "Misclassification error is the percentage mismatch of predcited vs actuals,
# irrespective of 1’s or 0’s. The lower the misclassification error, the better
# is your model."
misClassError(df_log_reg$test_data[[1]]$target, df_log_reg_model_predictions, threshold = opt_cut_off)
misClassError(df_log_reg_first_paper$test_data[[1]]$target, df_log_reg_first_paper_model_predictions, threshold = opt_cut_off_first_paper)

#### ROC ####
# "Receiver Operating Characteristics Curve traces the percentage of true
# positives accurately predicted by a given logit model as the prediction
# probability cutoff is lowered from 1 to 0. For a good model, as the cutoff is
# lowered, it should mark more of actual 1’s as positives and lesser of actual
# 0’s as 1’s. So for a good model, the curve should rise steeply, indicating
# that the TPR (Y-Axis) increases faster than the FPR (X-Axis) as the cutoff
# score decreases. Greater the area under the ROC curve, better the predictive
# ability of the model."
plotROC(df_log_reg$test_data[[1]]$target, df_log_reg_model_predictions)
plotROC(df_log_reg_first_paper$test_data[[1]]$target, df_log_reg_first_paper_model_predictions)

### Concordance ###
# "Ideally, the model-calculated-probability-scores of all actual Positive’s,
# (aka Ones) should be greater than the model-calculated-probability-scores of
# ALL the Negatives (aka Zeroes). Such a model is said to be perfectly
# concordant and a highly reliable one. This phenomenon can be measured by
# Concordance and Discordance. In simpler words, of all combinations of 1-0
# pairs (actuals), Concordance is the percentage of pairs, whose scores of
# actual positive’s are greater than the scores of actual negative’s. For a
# perfect model, this will be 100%. So, the higher the concordance, the better
# is the quality of model."
Concordance(df_log_reg$test_data[[1]]$target, df_log_reg_model_predictions)
Concordance(df_log_reg_first_paper$test_data[[1]]$target, df_log_reg_first_paper_model_predictions)

# sensitivity (True positive rate) - only 57%...
sensitivity(df_log_reg$test_data[[1]]$target, df_log_reg_model_predictions, threshold = opt_cut_off)
sensitivity(df_log_reg_first_paper$test_data[[1]]$target, df_log_reg_first_paper_model_predictions, threshold = opt_cut_off_first_paper)
# specificity (False positive rate) - 84%
specificity(df_log_reg$test_data[[1]]$target, df_log_reg_model_predictions, threshold = opt_cut_off)
specificity(df_log_reg_first_paper$test_data[[1]]$target, df_log_reg_first_paper_model_predictions, threshold = opt_cut_off_first_paper)

# confusion matrix
b <- confusionMatrix(df_log_reg$test_data[[1]]$target, df_log_reg_model_predictions, threshold = opt_cut_off)
b <- as.matrix(b)
a_first_paper <- confusionMatrix(df_log_reg_first_paper$test_data[[1]]$target, df_log_reg_first_paper_model_predictions, threshold = opt_cut_off_first_paper)
a_first_paper <- as.matrix(a_first_paper)
chisq.test(a_first_paper)
chisq.test(b)

fisher.test(a_first_paper)
fisher.test(b)
#df_log_reg$test_data[[1]]$target <- as.factor(df_log_reg$test_data[[1]]$target)
#oversample_model_predictions[[1]] <- as.factor(oversample_model_predictions[[1]])
# check 'balanced' data performance 
a <- confusionMatrix(factor(oversample_model_predictions[[1]]), df_log_reg$test_data[[1]]$target)
a <- as.matrix(a)

b <- confusionMatrix(undersample_model_predictions[[1]], df_log_reg$test_data[[1]]$target)
b <- as.matrix(b)

c <- confusionMatrix(both_model_predictions[[1]], df_log_reg$test_data[[1]]$target)
c <- as.matrix(c)
#d <- pmap(test, ~ confusionMatrix(.x, df_log_reg$test_data[[1]]$target))
chisq.test(a)
chisq.test(b)
chisq.test(c)
# test below seems to crash R for some reason
#fisher.test(a)
#fisher.test(b)
#fisher.test(c)

# check the same for the first paper data 
a <- confusionMatrix(factor(oversample_model_predictions_first_paper[[1]]), df_log_reg_first_paper$test_data[[1]]$target)
a <- as.matrix(a)

b <- confusionMatrix(undersample_model_predictions_first_paper[[1]], df_log_reg_first_paper$test_data[[1]]$target)
b <- as.matrix(b)

c <- confusionMatrix(both_model_predictions_first_paper[[1]], df_log_reg_first_paper$test_data[[1]]$target)
c <- as.matrix(c)
#d <- pmap(test, ~ confusionMatrix(.x, df_log_reg$test_data[[1]]$target))
chisq.test(a)
chisq.test(b)
chisq.test(c)

## merge lists together to reduce code 
balanced_data_list <- c(oversample_model_predictions, 
                        undersample_model_predictions, both_model_predictions) 
balanced_data_list_first_paper <- c(oversample_model_predictions_first_paper, 
                                    undersample_model_predictions_first_paper, both_model_predictions_first_paper)

# optimal cut off 
opt_cut_off_balanced <- map(balanced_data_list, ~ optimalCutoff(df_log_reg$test_data[[1]]$target, .x))
opt_cut_off_balanced_first_paper <- map(balanced_data_list_first_paper, ~ optimalCutoff(df_log_reg_first_paper$test_data[[1]]$target, .x))


#d <- pmap(test, ~ optimalCutoff(df_log_reg$test_data[[1]]$target, .x)) # why doesn't this work?

# missclasification errors - seem to be fairly similar and consitant?
map2(balanced_data_list, opt_cut_off_balanced, ~ misClassError(df_log_reg$test_data[[1]]$target, .x, threshold = .y))
map2(balanced_data_list_first_paper, opt_cut_off_balanced_first_paper, ~ misClassError(df_log_reg_first_paper$test_data[[1]]$target, .x, threshold = .y))


# ROCs 
map(balanced_data_list, ~ plotROC(df_log_reg$test_data[[1]]$target, .x))
map(balanced_data_list_first_paper, ~ plotROC(df_log_reg_first_paper$test_data[[1]]$target, .x))
# one of the "both" ROCs is very good but most are bad



# concordance 
map(balanced_data_list, ~ Concordance(df_log_reg$test_data[[1]]$target, .x))
map(balanced_data_list_first_paper, ~ Concordance(df_log_reg_first_paper$test_data[[1]]$target, .x))

# check the number of improvers and non-improvers 
table(df_log_reg$train_data[[1]]$target)
table(df_log_reg$test_data[[1]]$target)

table(df_log_reg_first_paper$train_data[[1]]$target)
table(df_log_reg_first_paper$test_data[[1]]$target)

# sensitivity (True positive rate) 
map2(balanced_data_list, opt_cut_off_balanced, ~ sensitivity(df_log_reg$test_data[[1]]$target, .x, threshold = .y))
# sensitivity of all are variable but the "both" method seems to be the most variable
map2(balanced_data_list_first_paper, opt_cut_off_balanced_first_paper, ~ sensitivity(df_log_reg_first_paper$test_data[[1]]$target, .x, threshold = .y))
# first paper sensitivity 
# sensitivity (True positive rate) - they are seem to range from 0.5 - 1?


# specificity (True negative rate) - 
map2(balanced_data_list, opt_cut_off_balanced, ~ specificity(df_log_reg$test_data[[1]]$target, .x, threshold = .y))
# specificity of all are a bit higher than the base data but the "both" method seems to give the best result
map2(balanced_data_list_first_paper, opt_cut_off_balanced_first_paper, ~ specificity(df_log_reg_first_paper$test_data[[1]]$target, .x, threshold = .y))



log_predict <- ifelse(df_log_reg_model_predictions > opt_cut_off, 1, 0)

#plot ROC using ROCR package 
pr <- ROCR::prediction(log_predict, df_log_reg$test_data[[1]]$target)
perf <-ROCR::performance(pr, measure = "tpr", x.measure = "fpr")
plot(perf) 
auc(df_log_reg$test_data[[1]]$target, log_predict)



# model calibration 

mletter <- function(r,n) {
  lower <-  2 + floor(log(r/(n+1))/log(2))
  upper <- -1 - floor(log((n+1-r)/(n+1))/log(2))
  i <- 2*r > n
  lower[i] <- upper[i]
  lower
}

actual <- rbinom(length(df_log_reg_model_predictions), 1, df_log_reg_model_predictions)
classes <- mletter(rank(df_log_reg_model_predictions), length(df_log_reg_model_predictions))
pgroups <- split(df_log_reg_model_predictions, classes)
agroups <- split(actual, classes)
bincounts <- unlist(lapply(pgroups, length)) # Bin populations
x <- unlist(lapply(pgroups, mean))           # Mean predicted values by bin
y <- unlist(lapply(agroups, mean))           # Mean outcome by bin

binprop <- bincounts / max(bincounts)
colors <- -log(binprop)/log(2)
colors <- colors - min(colors)
colors <- hsv(colors / (max(colors)+1))

abline(0,1, lty=1, col="Gray")                           # Reference curve
points(x,y, pch=19, cex = 3 * sqrt(binprop), col=colors) # Solid colored circles
points(x,y, pch=1, cex = 3 * sqrt(binprop))              # Circle outlines


#### below is based on code practicles from Keele modelling course ####
# Obtain the c statistic / AUC
library(pROC)
c1 <- roc(df_log_reg$test_data[[1]]$target~df_log_reg$predictions[[1]],ci=TRUE)
c1

## Another measure of separation
# Obtaining Somers' D statistic
Dstat <- rcorr.cens(df_log_reg$predictions[[1]],df_log_reg$test_data[[1]]$target)
Dstat["Dxy"]

# Confidence interval calculation is based on c-statistic and then converted to D as D=2c-1
# NOTE: D statistic ranges from -1 to 1 where -1 means all pairs disagree and 1
# means all pairs agree - so larger values are better.
(2*(Dstat['C Index']-(1.96*(Dstat['S.D.']/2))))-1
(2*(Dstat['C Index']+(1.96*(Dstat['S.D.']/2))))-1


# Calculate apparent calibration performance of the model in development data
library(DescTools)
# Estimate R squared
PseudoR2(df_log_reg$model[[1]],"Nagelkerke")

# And the Brier score - between 0 and 1 where 0 is perfect accuracy and 1 is worst perfomance
BrierScore(df_log_reg$model[[1]])

# CITL - note that the test data had to be used for this as they don't split data on the keele course
mod_log_1 <- glm(df_log_reg$test_data[[1]]$target~offset(df_log_reg$predictions[[1]]),family="binomial")
mod_log_1$coef
confint(mod_log_1)

# E/O - expected/observed ratio - not sure how to interpret
prop.table(table(df_log_reg$train_data[[1]]$target))
O <- prop.table(table(df_log_reg$train_data[[1]]$target))[2]
E <- mean(df_log_reg$predictions[[1]])
E/O


# Visual assessment of calibration by risk groups
# create objects they use to check code 
pred_prob <- df_log_reg$predictions[[1]]
DAY30 <- df_log_reg$test_data[[1]]$target
# create 10 risk groups
groups <- cut(pred_prob,breaks=quantile(pred_prob, prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),labels=c(1:10),include.lowest=TRUE)
# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(df_log_reg$test_data[[1]],groups,pred_prob)
obs <- (ddply(gpdata,~groups,summarise,mean=mean(as.numeric(DAY30)))[,2])-1
exp <- ddply(gpdata,~groups,summarise,mean=mean(pred_prob))
obsn <- table(DAY30,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# graph a simple calibration plot over 10 risk groups
par(pty="s")
plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),col="red",ylab="Observed",xlab="Expected")
plot(obs~exp[,2],col="red",ylab="Observed",xlab="Expected")

lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="green")
}
h <- hist(pred_prob, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(DAY30)-1)~pred_prob,span=1))
lines_data <- data.frame(pred_prob,obs_all)
lines_data2 <- lines_data[order(pred_prob),] 
lines(lines_data2[,1],lines_data2[,2],col="blue")
legend(0.0,0.9,c("Risk groups","Reference line","95% CI","Loess"),col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA),bty="n")





# You can obtain c-statistic and c-slope directly from Harrell's "validate" programme within 'rms'.
# However, CITL is not the intercept value - it is from the calibration model without the LP set as an offset
# So only use this code for c-statistics & c-slopes
library(rms)
mod1 <- lrm(df_log_reg$train_data[[1]]$target~ ., x=TRUE,y=TRUE)
set.seed(231398)
boot_1 <- validate(mod1,method="boot",B=500)
boot_1


##### ASSESSING KNN MODELS ######
d <- confusionMatrix(df_models$predictions[[10]], df_models$test_data[[10]]$target)
chisq.test(d)

print(df_models$model_fits[[10]])
df_models$target_name[[10]]
# filter to only have knn models 
knn_models_continuous <- dplyr::filter(df_models, modelName == "knn_model")
knn_models_continuous$model_fits <- set_names(knn_models_continuous$model_fits, knn_models_continuous$target_name)
map(knn_models_continuous$model_fits, ~ print(.x))
map(knn_models_continuous$model_fits, ~ plot(.x))


# get model summaries from "balanced" data
map(test$oversample_knn_model, ~ print(.x))
map(test$undersample_knn_model, ~ print(.x))
map(test$both_knn_model, ~ print(.x))

map(test$oversample_knn_model, ~ plot(.x))
map(test$undersample_knn_model, ~ plot(.x))
map(test$both_knn_model, ~ plot(.x))

oversample_confusion <- map(test$oversample_knn_model_predictions, ~ confusionMatrix(.x, df_log_reg$test_data[[1]]$target))

test$oversample_knn_model_predictions[[1]]
table(test$oversample_knn_model_predictions[[1]], df_log_reg$test_data[[1]]$target)




