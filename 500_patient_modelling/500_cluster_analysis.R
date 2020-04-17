
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

#load data - Note: if you load the .csv files instead it changes the format of several columns
load('data/df_locf_first_blooc_complete_2019-10-21.rda')
load('data/df_complete_first_blood_2019-10-21.rda')

# create version of data that contains improver col for logistic regression 
#df_locf_first_cc_logistic <- df_locf_first_blood_cc
# remove improver col from data for linear regression 
#df_locf_first_blood_cc$improver.1 <- NULL
#df_cc_first_blood$improver.1 <- NULL

# create vector of model targets
target_cols <- dplyr::select(df_locf_first_blood_cc, contains("_motor"), contains("_sensor"), contains("_scim"))
target_cols <- dplyr::select(target_cols, contains("discharge"), contains("year1"))
# add improver col
target_cols <- cbind(target_cols, df_locf_first_blood_cc$improver.1)
target_cols <- target_cols %>% dplyr::rename(improver.1 = `df_locf_first_blood_cc$improver.1`)
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

# create vector of initial neurology columns 
initial_cols <- names(dplyr::select(df_locf_first_blood_cc, contains("initial")))


# create df of bloods and intital neurology for replicating factor analysis
df_blood_initial_neurology <- df_locf_first_blood_cc[, c(remaining_bloods, initial_cols)]
names(df_blood_initial_neurology)

# create kendal cor matrix 
bloods_cor <- cor(df_blood_initial_neurology, method = "kendall") 
bloods_cor

# filter to correlations between 0.3 and 0.8 
# Note: check if it was 0.7 in first paper

#convert negative to absolute values for filtering 
bloods_cor <- abs(bloods_cor)
#convert matrix to long format 
bloods_cor_long <- reshape2::melt(bloods_cor)
#filter to correlations between 0.3 and 0.8
bloods_cor_filtered <- bloods_cor_long %>% 
  dplyr::filter(between(value, 0.3, 0.8))

# see what is left with these correlations
remaining_cors <- distinct(bloods_cor_filtered, Var1)
#NOTE: the variables with these correlations don't match those from the first
#paper

# convert back to wide matrix 
bloods_cor_filtered_wide <- acast(bloods_cor_filtered, Var1 ~ Var2)
# save correlated cols as character vector for column selection
kendall_cor_cols <- as.character(remaining_cors$Var1)

#Factor analysis of the data
factors_data <- fa(r = bloods_cor, nfactors = 6, rotate = "oblimin", fm = "pa")
#Getting the factor loadings and model analysis 
factors_data

# plot histogram of every column to check for errors 
#hist.data.frame(predictor_cols) # this line halts the script as it requires a click to see second page of graphs




pca_test <- predictor_cols
pca_test$improver <- target_cols$improver.1



# PCA on predictor_cols 
pca_all_features <- prcomp(predictor_cols)
summary(pca_all_features)
screeplot(pca_all_features, type = 'lines')
# NOte: the 1st pca describes 86% of the varience in the data. PC1 and PC2 describe 95%...

# try pca on bloods only 
pca_bloods_only <- prcomp(pca_test[,c(remaining_bloods, "improver")])
summary(pca_bloods_only)
screeplot(pca_bloods_only, type = 'lines')

remaining_bloods

iris_pca_df <- data.frame(first = pca_all_features$x[,1], second = pca_all_features$x[,2], improver = target_cols$improver.1)

ggplot(iris_pca_df, aes(x = first, y = second, colour = factor(improver))) +
  geom_point() +
  theme_bw() +
  ggtitle('Principle Component Analysis') +
  scale_x_continuous('First Principle Component') + scale_y_continuous('Second Principle Component')

# umap 
umap_test <- umap(predictor_cols)
umap_df <- data.frame(umap_test$layout)

ggplot(umap_df, aes(x = X1, y = X2, color = factor(predictor_cols$admission_asia), size = predictor_cols$initial_scim, alpha = factor(predictor_cols$neuro_level.L))) +
  geom_point()

# create character vector for umap plotting 
all_cols <- names(predictor_cols)
cols_no_blood <- setdiff(all_cols, remaining_bloods)
# remove admission asia
remove <- "admission_asia"
cols_no_blood <- setdiff(cols_no_blood, remove)
# add neuro c level injury 
cols_no_blood[(length(cols_no_blood)+1)] <- "neuro_level.C"
cols_no_blood <- set_names(cols_no_blood, cols_no_blood)
#cols_no_blood <- as.list(cols_no_blood)


umap_plot_df <- predictor_cols
umap_plot_df$x1 <- umap_df$X1
umap_plot_df$x2 <- umap_df$X2
# add cerviacal injury cols back
umap_plot_df$neuro_level.C <- ifelse(umap_plot_df$neuro_level.T == 0 & 
                                       umap_plot_df$neuro_level.S == 0 & 
                                       umap_plot_df$neuro_level.L == 0, 1, 0)

# umap plot function 
umap_plot <- function(df, x) {
  ggplot(df, aes(x = x1, y = x2, color = factor(admission_asia), size = df[[x]])) +
    geom_point() +
    labs(color = "Admin ASIA", size = x) +
    ggtitle("Umap plots - all variables") +
    theme_bw()
}

# create plots of all variables
umap_all_var_plots <- map(.x = cols_no_blood, .f = ~ umap_plot(umap_plot_df, .x))
# show plots 
umap_all_var_plots


set.seed(333)

# umap with just bloods
umap_test <- umap(pca_test[,c(remaining_bloods, "improver")])
umap_test <- umap(pca_test[,c(remaining_bloods)])
umap_df_bloods <- data.frame(umap_test$layout)

umap_plot_df$x1 <- umap_df_bloods$X1
umap_plot_df$x2 <- umap_df_bloods$X2

umap_bloods_plot <- map(.x = cols_no_blood, .f = ~ umap_plot(umap_plot_df, .x))






ggplot(umap_df, aes(x = X1, y = X2, size = factor(target_cols$improver), color = factor(predictor_cols$admission_asia))) +
  geom_point()

umap_df$admin_asia <- as.factor(predictor_cols$admission_asia)
umap_df$improver <- as.factor(target_cols$improver)




summarised_umap_df <- umap_df %>%
  dplyr::group_by(admin_asia, improver) %>%
  dplyr::summarise(mean_x1 = mean(X1), mean_x2 = mean(X2))

t.test(data = summarised_umap_df, mean_x2 ~ improver)

anova(data = summarised_umap_df, mean_x2 ~ admin_asia)

library(plotly)

plot_ly(umap_df, x = ~X1, y = ~X2, color = factor(predictor_cols$admission_asia), fillcolor = factor(target_cols$improver), type = "scatter")

# linear discriminant analysis 
library(MASS)
lda_all_features <- lda(target_cols$improver.1 ~., predictor_cols) # gives error about collinearity
lda_all_features

lda_all_features_predict <- predict(lda_all_features, predictor_cols)

a <- df_locf_first_blood_cc
a$predict <- lda_all_features_predict$class
head(a)
a$improver.1 <- as.factor(a$improver.1)
a$predict <- as.factor(a$predict)
# has accuracy of 84% - but was build and tested with same data, how would one
# know it isn't overfited?
caret::confusionMatrix(a$improver.1, a$predict)

