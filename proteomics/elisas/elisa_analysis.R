# load packages
library(tidyverse)
library(drc)
library(readxl)

# load data
a2m <- read.csv("../data/elisa_a2m_long_2020-10-27.csv")
rbp4 <- read.csv("../data/elisa_rbp4_long_2020-10-27.csv")
saa1 <- read.csv("../data/elisa_saa1_long_2020-10-27.csv")
head(a2m)
str(a2m)
# add standard concentration
standard_conc_a2m <- c(40, 20, 10, 5, 2.5, 1.25, 0.625, 0)
standard_conc_rbp4 <- c(1500, 750, 375, 188, 93.8, 46.9, 23.4, 0)
standard_conc_saa1 <- c(100, 50, 25, 12.5, 6.25, 3.13, 1.56, 0)

# define well content by row
well_content <- c(
  "standard", "standard_rep", "sample1_spike_sample", "sample1_spike_control",
  "sample1_unspiked_sample", "sample2_spike_sample", "sample2_spike_control",
  "sample2_unspiked_sample", "empty", "empty", "empty", "empty"
)
well_content2 <- c(
  "standard", "standard_rep", "sample1_spike_sample_rep", "sample1_spike_control_rep",
  "sample1_unspiked_sample_rep", "sample2_spike_sample_rep", "sample2_spike_control_rep",
  "sample2_unspiked_sample_rep", "empty", "empty", "empty", "empty"
)
# well_content2 <- c("standard", "standard_rep", rep("empty", 10))
# define dilutions
dilutions <- rep(c(1, 2, 4, 8), 2)
dilutions_saa1 <- c(10, 100, 500, 1000, rep(NA, 4))
# stop scientific notation
options(scipen = 999)
# function for adding standart concentrations and generating standard curve model
update_data <- function(data, standard_concentration, standard_units, dilutions) {
  # rename absorbance col
  names(data)[4] <- "absorbance"
  # update content col with relevent well content
  data$Content <- c(rep(well_content, 4), rep(well_content2, 4))
  # add units
  data$units <- standard_units
  # initiate standard concentration and dilution columns
  data$standard_con <- NA
  data$dilution <- NA
  # add standard concentration to relevent rows
  for (i in 1:8) {
    data$standard_con <- ifelse(grepl("standard", data$Content) &
      data$Well.Row == LETTERS[i], standard_concentration[i], data$standard_con)
    data$dilution <- ifelse(grepl("sample", data$Content) & data$Well.Row == LETTERS[i],
      dilutions[i], data$dilution
    )
  }

  # get average blank absorbance to subtract
  mean_blank <- data %>%
    dplyr::filter(standard_con == 0) %>%
    summarise(mean = mean(absorbance))

  data$absorbance_zero_adjust <- data$absorbance - mean_blank$mean

  # add indicator for when absorbance_zero_adjust is in standard curve
  standard_rows <- data %>%
    dplyr::filter(grepl("standard", data$Content))
  data$standard_range_indicator <- ifelse(data$absorbance_zero_adjust <=
    max(standard_rows$absorbance_zero_adjust)
  & data$absorbance_zero_adjust >=
      min(standard_rows$absorbance_zero_adjust), 1, 0)

  return(data)
}

rbp4 <- update_data(rbp4, standard_conc_rbp4, "pg/ml", dilutions)
a2m <- update_data(a2m, standard_conc_a2m, "ng/ml", dilutions)
saa1 <- update_data(saa1, standard_conc_saa1, "ng/ml", dilutions)
head(rbp4)

# save tables
write.csv(rbp4, "../data/rbp4_elisa_long_updated_2020-11-02.csv")
write.csv(a2m, "../data/a2m_elisa_long_updated_2020-11-02.csv")
write.csv(saa1, "../data/saa1_elisa_long_updated_2020-11-02.csv")
test <- pivot_wider(rbp4, names_from = Well.Col, values_from = c(absorbance_zero_adjust, Well.Row))

# get standard curve and plot samples on it
plot_standard_curve <- function(data, plot_file_name, plot_title) {
  # filter out blanks
  data <- data %>%
    dplyr::filter(standard_con != 0 | is.na(standard_con))
  # generate standard curve
  model1 <- drm(absorbance_zero_adjust ~ standard_con,
    fct = LL.4(names = c("Slope", "Lower", "Upper", "ED50")),
    data = data
  )
  # save plot as pdf
  pdf(plot_file_name)
  # paste contents of content and dilution
  cols <- c("Content", "dilution")
  data$pasted_col <- do.call(paste, c(data[cols], sep = "_"))

  plot(model1,
    main = plot_title,
    xlab = paste("Concentration", data$units[1], sep = " "),
    ylab = "Absorbance (450 nm)"
  )

  # filter to absorbances of samples
  response <- data %>%
    dplyr::filter(!grepl("standard", data$Content) & !grepl("empty", data$Content))
  response
  # get mean absorbance of replicates
  mean_abs <- tapply(response$absorbance_zero_adjust, response$pasted_col, mean)

  # Estimate the concentration
  DOSEx <- ED(model1, mean_abs, type = "absolute", display = F)
  # And here are the estimates including standard errors
  DOSEx

  # We can add those to our plot
  points(y = mean_abs, x = DOSEx[, 1], col = "blue", pch = 19, cex = 2)
  # With error bars
  arrows(DOSEx[, 1], mean_abs, DOSEx[, 1] + DOSEx[, 2] * 1.96, mean_abs, length = 0.1, angle = 90, lwd = 3, col = "blue")
  arrows(DOSEx[, 1], mean_abs, DOSEx[, 1] - DOSEx[, 2] * 1.96, mean_abs, length = 0.1, angle = 90, lwd = 3, col = "blue")

  # add sample labels to plot
  # with(data, text(absorbance_zero_adjust~standard_con, labels = Content, pos = 4))
  dev.off()
  return(mean_abs)
}
# use function
plot_standard_curve(rbp4, "../figures/rbp4_standard_curve_2020-11-02.pdf", "RBP4 Standard curve")
plot_standard_curve(a2m, "../figures/a2m_standard_curve_2020-11-02.pdf", "A2M Standard curve")
plot_standard_curve(saa1, "../figures/saa1_standard_curve_2020-11-02.pdf", "SAA1 Standard curve")

# ggplot version of above function
ggplot_version <- function(data, plot_file_name, plot_title, sample_dilution_range, sample_exclude) {
  require(ggrepel)
  # filter out blanks
  data <- data %>%
    dplyr::filter(standard_con != 0 | is.na(standard_con))
  # generate standard curve
  model1 <- drm(absorbance_zero_adjust ~ standard_con,
    fct = LL.4(names = c("Slope", "Lower", "Upper", "ED50")),
    data = data
  )
  # paste contents of content and dilution
  cols <- c("Content", "dilution")
  data$pasted_col <- do.call(paste, c(data[cols], sep = "_"))

  # filter to absorbances of samples
  response <- data %>%
    dplyr::filter(!grepl("standard", data$Content) & !grepl("empty", data$Content))
  response
  # get mean absorbance of replicates
  mean_abs <- tapply(response$absorbance_zero_adjust, response$pasted_col, mean)

  # Estimate the concentration
  DOSEx <- ED(model1, mean_abs, type = "absolute", display = F)
  # And here are the estimates including standard errors
  DOSEx
  # filter to where standard_con is not 0
  data_filter <- data %>%
    dplyr::filter(standard_con != 0)
  # add model predictions to data
  data.predict <- cbind(data_filter, predict(model1, interval = "confidence"))

  # pivot response wide
  response_wide <- pivot_wider(response, names_from = pasted_col, values_from = absorbance_zero_adjust)
  # get sample means
  sample_means <- response_wide %>%
    summarise_at(vars({{ sample_dilution_range }}), mean, na.rm = TRUE)
  # convert sample means to long format
  sample_means_df <- as.data.frame(t(sample_means)) %>%
    rownames_to_column(var = "Content")
  names(sample_means_df)[2] <- "absorbance_zero_adjust"
  # add predicted sample conc and conf intervals
  sample_means_df$standard_con <- DOSEx[, 1]
  sample_means_df$Lower <- DOSEx[, 2]
  sample_means_df$Upper <- DOSEx[, 2]
  sample_means_df$well_type <- "sample"
  # add well type to data.predict
  data.predict$well_type <- "standard"
  # filter to predicted sample concentration less than the highest standard
  sample_filter <- sample_means_df %>%
    dplyr::filter(standard_con < max(data.predict$standard_con) & absorbance_zero_adjust < max(data.predict$absorbance_zero_adjust))
  # bind rows
  test <- dplyr::bind_rows(sample_filter, data.predict)

  # filter to a specific sample
  test$sample_num <- ifelse(grepl(sample_exclude, test$Content), 1, 0)
  test <- test %>%
    dplyr::filter(sample_num == 0)
  # add col to denote well type
  test <- test %>%
    dplyr::mutate(Sample = case_when(
      grepl("spike_sample", test$Content) ~ "Spiked Sample",
      grepl("spike_control", test$Content) ~ "Spiked Control",
      grepl("unspiked_sample", test$Content) ~ "Unspiked Sample",
      TRUE ~ "standard"
    ))
  plot <- test %>%
    ggplot(aes(x = standard_con, y = absorbance_zero_adjust, colour = well_type)) +
    geom_point() +
    # add line for standard
    geom_line(aes(x = standard_con, y = Prediction)) +
    # add confidence intervals to line
    geom_ribbon(
      data = dplyr::filter(test, well_type == "standard"),
      aes(ymin = Lower, ymax = Upper), alpha = 0.3
    ) +
    # add error bars to sample concentrations
    # geom_errorbar(data = dplyr::filter(test, well_type == "sample"), aes(ymin = Lower, ymax = Upper)) +
    # add labels for samples
    geom_label_repel(
      data = dplyr::filter(test, well_type == "sample"),
      aes(label = Content, fill = factor(Sample)),
      color = "black",
      size = 2.1
    ) +
    ggtitle(plot_title) +
    xlab(paste("Concentration", data.predict$units[1], sep = " ")) +
    ylab("Absorbance (450 nm)") +
    theme_bw()

  plot
  ggsave(plot_file_name)
  return(test)
}

# note that rbp4 has 3 observations that get cut due to the high absorbance, but their estimated concentration is in line
ggplot_df <- ggplot_version(
  rbp4, "../figures/rbp4_standard_curve_ggplot_sample_1_2020-11-02.pdf",
  "RBP4 Standard curve",
  sample1_spike_sample_1:sample2_unspiked_sample_rep_8,
  sample_exclude = "sample2"
)
ggplot_version(
  rbp4, "../figures/rbp4_standard_curve_ggplot_sample_2_2020-11-02.pdf",
  "RBP4 Standard curve",
  sample1_spike_sample_1:sample2_unspiked_sample_rep_8,
  sample_exclude = "sample1"
)
ggplot_version(
  a2m, "../figures/a2m_standard_curve_ggplot_sample1_2020-11-02.pdf",
  "A2M Standard curve",
  sample1_spike_sample_1:sample2_unspiked_sample_rep_8,
  sample_exclude = "sample2"
)
ggplot_version(
  a2m, "../figures/a2m_standard_curve_ggplot_sample2_2020-11-02.pdf",
  "A2M Standard curve",
  sample1_spike_sample_1:sample2_unspiked_sample_rep_8,
  sample_exclude = "sample1"
)
ggplot_version(
  saa1, "../figures/saa1_standard_curve_ggplot_sample1_2020-11-02.pdf",
  "SAA1 Standard curve",
  sample1_spike_sample_1:sample2_unspiked_sample_rep_8,
  sample_exclude = "sample2"
)
ggplot_version(
  saa1, "../figures/saa1_standard_curve_ggplot_sample2_2020-11-02.pdf",
  "SAA1 Standard curve",
  sample1_spike_sample_1:sample2_unspiked_sample_rep_8,
  sample_exclude = "sample1"
)

# testing:
# pivot response wide
response_wide <- pivot_wider(response, names_from = pasted_col, values_from = absorbance_zero_adjust)
# get sample means
sample_means <- response_wide %>%
  summarise_at(vars(sample1_500:sample4_10000), mean, na.rm = TRUE)
# convert sample means to long format
sample_means_df <- as.data.frame(t(sample_means)) %>%
  rownames_to_column(var = "Content")
names(sample_means_df)[2] <- "absorbance_zero_adjust"
# add predicted sample conc and conf intervals
sample_means_df$standard_con <- DOSEx[, 1]
sample_means_df$Lower <- DOSEx[, 2]
sample_means_df$Upper <- DOSEx[, 2]
sample_means_df$well_type <- "sample"
# add well type to data.predict
data.predict$well_type <- "standard"
# filter to predicted sample concentration less than the highest standard
sample_filter <- sample_means_df %>%
  dplyr::filter(standard_con < max(data.predict$standard_con))
# bind rows
test <- dplyr::bind_rows(sample_filter, data.predict)

plot <- ggplot_df %>%
  ggplot(aes(x = standard_con, y = absorbance_zero_adjust, colour = well_type)) +
  geom_point() +
  geom_line(aes(x = standard_con, y = Prediction)) +
  geom_ribbon(data = dplyr::filter(ggplot_df, well_type == "standard"), aes(ymin = Lower, ymax = Upper), alpha = 0.3) +
  # geom_errorbar(data = dplyr::filter(test, well_type == "sample"), aes(ymin = Lower, ymax = Upper)) +
  geom_text(data = dplyr::filter(ggplot_df, well_type == "sample"), aes(label = Content), nudge_y = 0.02) +
  ggtitle("rbp4") +
  xlab(paste("Concentration", "units", sep = " ")) +
  ylab("Absorbance (450 nm)") +
  theme_bw()
plot

?geom_text
