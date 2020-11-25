## On macOS there might be problems installing from Bioconductor
## via KNIME. Make sure you installed it manually in the referenced
## R library used by KNIME (see Settings->R).
options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = FALSE)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# BiocManager::install(version = "3.9")
# BiocManager::install("MSstats")
# BiocManager::install("MSstatsTMT")

library(MSstatsTMT) # for analysis
library(MSstats) # for plotting
library(tidyverse) # for data wrangling
library(ggrepel) # for arranging plot labels

## file <- "../../../../proteomic_data/knime_outputs/all_fraction_both_runs_itraq_workflow_2020-10-23/both_runs_msstats.csv"
## file <- "../../../../proteomic_data/knime_outputs/mstats_run2.csv"
## file <- "../../../../proteomic_data/knime_outputs/fractions_merged_itraq_workflow/msstatstmt_both_runs_2020-10-27.csv"
file <- "msstatstmt_both_runs_2020-10-27.csv"

output_directory <- "plots/"
print(output_directory)

#################################################################################
########################## MSstatsTMT converter for OpenMS output ###############
#################################################################################
msstatstmt_data <- read.csv(file)

data <- msstatstmt_data
## check data
colnames(data)
head(data)
length(unique(data$ProteinName)) ## 494 proteins
sum(is.na(data$Intensity)) ## no NAs
sum(data$Intensity == 0) ## ~14K with intensity of 0
unique(data$Condition)
## compare RetentionTime and Reference
length(unique(data$RetentionTime)) # 19855
length(unique(data$Reference)) # 19786

## combine the fractions for each mixture
processed.data <- OpenMStoMSstatsTMTFormat(data)
head(processed.data)
dim(processed.data)
length(unique(processed.data$ProteinName)) ## 270 proteins
sum(processed.data$Intensity == 0) ## none with 0 anymore

#################################################################################
########################## protein summarization ################################
#################################################################################
## protein summarization: msstats + normalization
quant.data <- proteinSummarization(processed.data,
  method = "msstats",
  global_norm = FALSE, # don't want to use this where there may be biological
  # variation in protein level or in cases of a small number of proteins
  reference_norm = FALSE, # we don't have a reference channel in our experiment
  MBimpute = TRUE,
  maxQuantileforCensored = NULL,
  remove_norm_channel = TRUE,
  remove_empty_channel = TRUE
)

head(quant.data)
# save(quant.data, file='quant.data.rda')

dataProcessPlotsTMT(
  data.peptide = processed.data,
  data.summarization = quant.data,
  type = "ProfilePlot",
  width = 21,
  height = 7,
  address = "../data/both_runs_"
)
#################################################################################
########################## statistical testing ##################################
#################################################################################
## prepare contrast matrix
unique(quant.data$Condition)

comparison <- matrix(c(
  -1, 0, 0, 1,
  0, -1, 0, 1,
  0, 0, -1, 1,
  -1, 1, 0, 0,
  0, 1, -1, 0
), nrow = 5, byrow = T)

# Set the names of each row
row.names(comparison) <- contrasts <- c(
  "c_acute_improve-c_subacute_improve",
  "c_acute_improve-c_subacute_nonimprove",
  "c_acute_improve-c_acute_nonimprove",
  "c_subacute_nonimprove-c_subacute_improve",
  "c_subacute_nonimprove-c_acute_nonimprove"
)
# Set the column names
colnames(comparison) <- c(
  "c_subacute_improve",
  "c_subacute_nonimprove",
  "c_acute_nonimprove",
  "c_acute_improve"
)

comparison

## do pairwise comparison
data.res <- groupComparisonTMT(
  data = quant.data,
  contrast.matrix = "pairwise",
  moderated = TRUE, # do moderated t test
  adj.method = "BH"
) # multiple comparison adjustment

## data.res <- groupComparisonTMT(
##   data = quant.data,
##   contrast.matrix = comparison,
##   moderated = TRUE, # do moderated t test
##   adj.method = "BH"
## ) # multiple comparison adjustment

head(data.res) # all the stats are NAs....
# perhaps there's something wrong with the experimental design file?


data.res <- data.res %>% filter(!is.na(Protein))
# save(data.res, file = "testing.data.rda")
prots <- as.character(data.res$Protein)
isups <- sapply(prots, function(x) {
  grepl(x, pattern = "ups", fixed = TRUE)
})
prots[isups] <- "ups"
prots[!isups] <- NA

data.res.mod <- data.res
data.res.mod$Protein <- prots
head(data.res.mod)
library(MSstats)
# use volcano plot from msstats package
groupComparisonPlots(data = data.res, type = "VolcanoPlot", address = paste(output_directory, "pairwise_", sep = ""), sig = 0.05)
## groupComparisonPlots(data = data.res, type = "VolcanoPlot", address = paste(output_directory, "c_acute_improve-c_subacute_improve_", sep = ""), which.Comparison = "c_acute_improve-c_subacute_improve", sig = 1e-5) # 1e-5
## groupComparisonPlots(data = data.res, type = "VolcanoPlot", address = paste(output_directory, "c_acute_improve-c_subacute_nonimprove_", sep = ""), which.Comparison = "c_acute_improve-c_subacute_nonimprove", sig = 1e-5) # 1e-5
## groupComparisonPlots(data = data.res, type = "VolcanoPlot", address = paste(output_directory, "c_acute_improve-c_acute_nonimprove_", sep = ""), which.Comparison = "c_acute_improve-c_acute_nonimprove", sig = 0.05) # 1e-5
## groupComparisonPlots(data = data.res, type = "VolcanoPlot", address = paste(output_directory, "c_subacute_nonimprove-c_subacute_improve_", sep = ""), which.Comparison = "c_subacute_nonimprove-c_subacute_improve", sig = 0.05) # 1e-5
## groupComparisonPlots(data = data.res, type = "VolcanoPlot", address = paste(output_directory, "c_subacute_nonimprove-c_acute_nonimprove_", sep = ""), which.Comparison = "c_subacute_nonimprove-c_acute_nonimprove", sig = 1e-5) # 1e-5
## groupComparisonPlots(data = data.res.mod, type = "VolcanoPlot", address = paste(output_directory, "0125-05_", sep = ""), which.Comparison = "0125-05", sig = 0.05)

## # to show in view output
## groupComparisonPlots(data = data.res.mod, type = "VolcanoPlot", address = F, which.Comparison = "0125-05", sig = 0.05)

# only 7 are significant
sum(data.res$adj.pvalue <= 0.05, na.rm = TRUE)
data.res %>% filter(adj.pvalue <= 0.05)



### Manual attempt at analysing the data
head(msstatstmt_data)
## add log2 transformed intensity
msstatstmt_data$log2_intensity <- log2(msstatstmt_data$Intensity)
## rename run col
msstatstmt_data$Run <- ifelse(msstatstmt_data$Run == "1_1_1", 1, 2)
## get summary stats
msstats_summary <- msstatstmt_data %>%
  dplyr::group_by(ProteinName, Run, Condition) %>%
  dplyr::summarise(
    mean = mean(log2_intensity),
    sd = sd(log2_intensity),
    len = n()
  )
## add standard error
msstats_summary <- mutate(msstats_summary,
  se = sqrt(sd^2 / len),
  lower = mean - qt(0.975, len - 1) * se,
  upper = mean + qt(0.975, len - 1) * se
)

msstats_summary


## save a plot of each protein to a single pdf
pdf("plots/all_proteins_2020-11-16.pdf", width = 6, height = 5)
## filter out where mean is infinate
inf_filter <- msstats_summary %>%
  dplyr::filter(!is.infinite(mean))

## get list of all proteins
prolist_filtered <- unique(inf_filter$ProteinName)

for (i in 1:length(prolist_filtered)) {
  ## filter to a single protein
  sub <- inf_filter %>%
    filter(ProteinName %in% prolist_filtered[i])

  ## plot the mean log2_intensity with error bars for each condition
  p <- ggplot(
    data = sub,
    aes(x = Condition, y = mean, shape = factor(Run), colour = Condition)
  ) +
    geom_point(position = position_dodge(width = 0.3)) +
    labs(
      title = unique(sub$ProteinName),
      y = "mean of Log2(Intensity)",
      colour = "",
      shape = "Run"
    ) +
    scale_x_discrete(guide = guide_axis(n.dodge = 3)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(colour = "black"),
      legend.key = element_rect(colour = "white")
    )

  p <- p + geom_errorbar(aes(ymax = upper, ymin = lower),
    width = 0.1,
    position = position_dodge(width = 0.3)
  )

  print(p)
}

dev.off()

## function to perfom a ttest on each protein and then get the log2 foldchange
fold_change_ttest <- function(df, conditions, run) {
  require(tidyverse)
  require(ggrepel)
  ## t.test for each protein
  ## gives error about not enough x observations

  # first, set up the table to record the results
  allresult <- data.frame(
    Protein = NULL,
    tstat = NULL,
    pvalue = NULL
  )
  ## filter to relevant run
  df <- df %>%
    dplyr::filter(Run == run)
  ## get list of all proteins
  prolist <- unique(df$ProteinName)

  for (i in 1:length(prolist)) {
    sub <- df %>%
      dplyr::filter(ProteinName %in% prolist[i]) %>%
      dplyr::filter(Condition %in% conditions) %>%
      dplyr::filter(!is.infinite(log2_intensity))

    ## use tryCatch to prevent loop from stopping on error
    ## the error seems to occur due to difference in
    ## number of values in the comparision
    tryCatch(
      {
        ## perfom t.test
        result <- t.test(log2_intensity ~ Condition,
          data = sub
        )
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )

    ## skip data.frame grouping error
    tryCatch(
      {
        tmp <- data.frame(
          Protein = unique(sub$Protein),
          log2fc = result$estimate[1] - result$estimate[2],
          tstat = result$statistic,
          pvalue = result$p.value
        )
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )

    ## build a dataframe of results
    allresult <- rbind(
      allresult,
      tmp
    )
    print(i)
  }

  head(allresult)
  ## adjust p value
  allresult$adj.pvalue <- p.adjust(allresult$pvalue, method = "BH")

  return(allresult)
}

## make list of comparisons of interest
conditions <- unique(msstatstmt_data$Condition)
conditions
comparisons_of_interest <- list(
  "acute_c_improvers_vs_nonimprovers_run1" = c("c_acute_improve", "c_acute_nonimprove"),
  "subacute_c_improvers_vs_nonimprovers_run1" = c("c_subacute_improve", "c_subacute_nonimprove"),
  "acute_c_improvers_vs_nonimprovers_run2" = c("c_acute_improve", "c_acute_nonimprove"),
  "a_vs_d_run2" = c("a", "d"),
  "acute_c_improvers_vs_a_run2" = c("c_acute_improve", "a"),
  "acute_c_improvers_vs_d_run2" = c("c_acute_improve", "d"),
  "acute_c_nonimprovers_vs_a_run2" = c("c_acute_improve", "a"),
  "acute_c_nonimprovers_vs_d_run2" = c("c_acute_improve", "d")
)
## make lists of relevant runs
runs <- c(rep(1, 2), rep(2, 6))
## use function
ttests_fold_changes <- map2(
  comparisons_of_interest, runs,
  ~ fold_change_ttest(msstatstmt_data, .x, .y)
)

## make plots and save them
pdf("plots/volcano_plots_manual_2020-11-16.pdf", width = 6, height = 5)

map2(ttests_fold_changes, names(ttests_fold_changes), ~
.x %>%
  mutate(Significance = ifelse(adj.pvalue < 0.05, "P < .05", "P ≥ .05")) %>%
  ## use log10 adjusted pval
  ggplot(aes(x = log2fc, y = -log10(adj.pvalue))) +
  ## add labels for significant proteins
  geom_text_repel(
    data = filter(.x, adj.pvalue < .05),
    mapping = aes(label = Protein)
  ) +
  geom_point(aes(color = Significance)) +
  ## add lines for significance cut off and +- 1 fold change
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted") +
  labs(
    x = expression(log[2] ~ fold - change),
    y = expression(-log[10] ~ adjusted ~ p - value),
    title = .y
  ) +
  theme_minimal())

dev.off()


## heatmaps
## need wide format for heatmaps
msstatstmt_data_wide <- msstatstmt_data %>%
  dplyr::filter(!is.infinite(log2_intensity)) %>%
  dplyr::select("ProteinName", "Run", "Condition", "log2_intensity") %>%
  pivot_wider(
    names_from = c(Condition, Run),
    values_from = log2_intensity, values_fn = mean # take mean of duplicates
  )

## get matrix without
crc <- t(as.matrix(msstatstmt_data_wide[, -1]))
crc[1:8, 1:10]
colnames(crc) <- msstatstmt_data_wide$ProteinName
crc[1:8, 1:10]
sum(is.na(crc))


## fill in NA values
crc <- apply(crc, 2, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))
## so many proteins it's hard to analyse
heatmap(crc)
## filter to proteins of interest
example_proteins <- ttests_fold_changes$acute_c_nonimprovers_vs_a_run2 %>%
  dplyr::filter(adj.pvalue < 0.05)


wide_filter <- msstatstmt_data_wide %>% filter(ProteinName %in% example_proteins$Protein)

crc.sub <- t(as.matrix(wide_filter[, -1]))
crc.sub
colnames(crc.sub) <- wide_filter$ProteinName
crc.sub
sum(is.na(crc.sub))
crc.sub <- apply(crc.sub, 2, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

heatmap(crc.sub)

condition <- substr(rownames(crc.sub), start = 13, stop = 19)
condition.color <- as.character(factor(condition, labels = c("red", "blue", "yellow", "green")))

heatmap(crc.sub, RowSideColors = condition.color)

library("RColorBrewer")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap(crc.sub,
  col = col,
  # scale = 'col',
  RowSideColors = condition.color
)
