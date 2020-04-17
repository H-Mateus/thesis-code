### Normality checks

library(tidyverse)
load("../../mateus/500_patient/data/df_long.rda")

# select adj calcium
calcium <- df_long %>%
  dplyr::filter(analyte == "adj_calc_est")
# summarise calcium: are mean and median similar??
calcium %>%
  summarise(count = n(), mean = mean(value, na.rm = TRUE),
            median = median(value), IQR = IQR(value))
# so median and mean similar, good start for normality!
# create a histogram:
calcium %>%
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = 'red', colour = 'black')
# so not too bad, but somewhat right skew
# next a qqnorm plot:
# dev.new()
qqnorm(calcium$value)
qqline(calcium$value)
# some are above the line at higher quantiles, not so good!
# shapiro.test
shapiro.test(calcium$value)
# p < 2.2 e-16, so significantly different to normal distribution; normality not agreed
# but shapiro test best for small numbers. Try
# ks test:
ks.test(calcium$value, "pnorm", mean = mean(calcium$value), sd = sd(calcium$value))
# shouldn't have ties:
# remove ties:
ks.test(unique(calcium$value), "pnorm", mean = mean(unique(calcium$value)), sd = sd(unique(calcium$value)))
# but ks is not particularly powerful
# standardised range??
q <- (max(calcium$value) - min(calcium$value)) / sd(calcium$value)
q
# Statistical tables show the critical value for the test statistic (q) to be between 3.34 and 4.71 when p = 0.05
# if q is more extreme (not in between those values), maybe normality assumption not reasonable, but this test is not very powerful

# now do the same with the wide data frame:
load("../../mateus/500_patient/data/df_wide.rda")
str(df_wide)
# shapiro
shapiro.test(df_wide$"adj_calc_est")
# so the same.

# how to do this automatically for all analytes??
# although there are less than 5000 values and shapiro wilk should work,
# the data frame contains more values for others, so when you do lapply, it doesn't
# work; better to use other tests:

lapply(df_wide[,3], shapiro.test)

# but try the anderson darling test instead:
library(nortest)
ad_tests <- lapply(df_wide[,3:49], ad.test)
# get the p.values
ad_tests_pvalues <- lapply(ad_tests, '[[', 'p.value')
# unlist
ad_tests_pvalues <- unlist(ad_tests_pvalues)
str(ad_tests_pvalues)
head(ad_tests_pvalues)
ad_tests_pvalues[ad_tests_pvalues >= 0.0001]
# so it doesn't look like normality assumption is reasonable in any!!


# how to do all methods in one go (except shapiro.test)??
str(df_wide)
library(dplyr)
library(nortest)
# define a functions that does everything (but not the saphiro-wilk test, which is for small numbers)
check_normality <- function(input_variable){
  # remove nas
  input_variable <- input_variable[!is.na(input_variable)]
  mean_input <-  mean(input_variable, na.rm = TRUE)
  median_input <- median(input_variable)
  max_input <- max(input_variable)
  min_input <- min(input_variable)
  range_input <- max_input - min_input
  sd_input <-  sd(input_variable)
  q_input <-  (max(input_variable) - min(input_variable)) / sd(input_variable)
  summary <- data.frame(mean_input, median_input, max_input, min_input, range_input, sd_input, q_input)
  norm_plot <- qqnorm(input_variable)
  histogram <- hist(input_variable)
  # need to use unique values for ks tests as ties give a warning
  kolmagorov_smirnov_test <- ks.test(unique(input_variable),
                                     'pnorm', mean = mean(unique(input_variable)), sd = sd(unique(input_variable)))
  anderson_darling_test <- ad.test(input_variable)
  answer <- list(summary, kolmagorov_smirnov_test, anderson_darling_test,
                 histogram, norm_plot)
  return(answer)
}

# see if it works:
result <- check_normality(df_wide[,4])
# to look at the contents:
result[[1]]
result[[2]]
result[[3]]
plot(result[[4]])
plot(result[[5]])

# check more / all:
results_all <- lapply(df_wide[, 3:49], check_normality)
str(results_all)
# get the first summary
results_all[[1]][[1]]
# and the third summary
results_all[[3]][[1]]
# and the last summary:
results_all[[length(results_all)]][[1]]
# to get the first histogram
plot(results_all[[1]][[4]])
# and the fourth normplot
plot(results_all[[4]][[5]])
# and the last normplot
plot(results_all[[length(results_all)]][[5]])

# get all the summaries in a list and unlist for example the q values
summaries <- lapply(results_all, '[[', 1)
# get the q values:
q_values_listed <- lapply(summaries, '[[', 'q_input')
q_values <- unlist(q_values_listed)
q_values

# get all the p-values of the ks tests:
ks_results <- lapply(results_all, '[[', 2)
ks_results_pvalues_listed <- lapply(ks_results, '[[', 'p.value')
ks_pvalues <- unlist(ks_results_pvalues_listed)
ks_pvalues <- ks_pvalues[order(-ks_pvalues)] # show in reverse order
normality_reasonable <- ks_pvalues[ks_pvalues >= 0.05]
normality_not_reasonable <- ks_pvalues[ks_pvalues <= 0.05]
normality_reasonable
normality_not_reasonable # so you need to have a closer look at these variables
# but ks is not so powerful!!!

# get all the p-values of the ad tests:
ad_results <- lapply(results_all, '[[', 3)
ad_results_pvalues_listed <- lapply(ad_results, '[[', 'p.value')
ad_pvalues <- unlist(ad_results_pvalues_listed)
ad_pvalues
# so none are normally distributed as per ad test!!

# can you transform data to normality?????
# use Box Cox transformation to find the best value for lambda:
library(MASS)
box <- boxcox(calcium$value ~ 1, lambda = seq(-10, 10, 0.1)) # try values from -10 to 10
cox <- data.frame(box$x, box$y) # data frame with results
cox2 = cox[with(cox, order(-cox$box.y)),] # order data frame by decreasing y
cox2[1,] # what is the best value?
lambda <- cox2[1, 'box.x'] # extract that value
transformed_calcium_values <- (calcium$value ^ lambda - 1) / lambda # transform the data
hist(transformed_calcium_values)
shapiro.test(transformed_calcium_values)
# so, transformation doesn't really help. You could also try log, sqrt and square: same result

# the reason the plot tests are probably not normally distributed is because patients are abnormal / ill
# the normal range and reference range described is for people who are not ill.
# bloods probably improve over time???
# perhaps check blood tests on discharge (when they are better) and see how they are distributed

# also bear in mind that some of the blood test are related:
# mcv (mean cell volume) = 10 * haematocrit(%) / red blood cell count (10^12 / L)
# mch (mean cell hb)= 10 * haemoglobin (g/dL) / red blood cell count (10^12 / L)
# mchc (mean cell hb concentration) = 100 * haemoglobin (g/dL) / haematocrit (%)
# have a look at:
# https://www.labce.com/spg579119_red_blood_cell_rbc_indices_definitions_and_calcula.aspx
# (so, if you do model forming, probably best to just take haemaglobin, haematocrit and red blood cell count)

# you can use this as a check:
mcv_df <- df_wide$'mean_cell_vol'
mcv_calc <- 10 * (df_wide$'haemocrit') / df_wide$'rbc_count'
mcv <- mcv_df - mcv_calc # this should be close to zero for all values (subject to rounding error)!!
summary(mcv) # but it isn't zero!!!
# check with a histogram:
hist(mcv) # so bimodal
# I suspect the dimension are different in some of the variables.
# you need to carefully check all the values to make sure they make sense.
# did you change dimensions???
# you can do that by looking at the histograms created above and check with the reference range
# then, redo the normality checks, as this may explain some of the problems!!
# for example, to look at the haematocrit (variable 22) histogram:
plot(results_all[[21]][[4]])
# which shows a problem.
# haematocrit is percentage of blood cells in the blood and is normally about 40%
# the histogram shows values as low a 20%, you could get this by taking blood out of a drip arm!
# More than 80% cells in your blood will give you serious problems and blood test is unlikely correct!
# it is impossible for the haematocrit to be over 100%!!
# as you can see, quite a few are!
# check out the variables like this

# have a look at the mch:
mch_df <- df_wide$'mean_cell_hb)'
mch_calc <- (df_wide$'haemoglob') / df_wide$'rbc_count'
mch <- mch_df - mch_calc
summary(mch)
hist(mch)
# this looks better, but there is an extreme value at 1!!
# so, some checking to do.

# just to show the anderson darling test works with random numbers from the normal distribution:
example <- rnorm(5000, mean = mean(df_wide$'mean_cell_hb', na.rm = TRUE), sd = sd(df_wide$'mean_cell_hb', na.rm = TRUE))
ad.test(example)
