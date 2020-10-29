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

library(dplyr) # for data wrangling
library(MSstatsTMT) # for analysis
library(MSstats) # for plotting

file <- "../../../../proteomic_data/knime_outputs/all_fraction_both_runs_itraq_workflow_2020-10-23/both_runs_msstats.csv"
file <- "../../../../proteomic_data/knime_outputs/mstats_run2.csv"
file <- "../../../../proteomic_data/knime_outputs/fractions_merged_itraq_workflow/msstatstmt_both_runs_2020-10-27.csv"

output_directory <- "../figures/msstatstmt_plots/"
print(output_directory)

#################################################################################
########################## MSstatsTMT converter for OpenMS output ###############
#################################################################################
MSstatsConverter_OpenMS_out <- read.csv(file)

data <- MSstatsConverter_OpenMS_out
colnames(data)

## compare RetentionTime and Reference
length(unique(data$RetentionTime)) # 19855
length(unique(data$Reference)) # 19786

## combine the fractions for each mixture
processed.data <- OpenMStoMSstatsTMTFormat(data)
head(processed.data)

#################################################################################
########################## protein summarization ################################
#################################################################################
## protein summarization: msstats + normalization
quant.data <- proteinSummarization(processed.data,
  method = "msstats",
  global_norm = FALSE, # don't want to use this where there may be biological variation in protein level or in cases of a small number of proteins
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

data.res <- groupComparisonTMT(
  data = quant.data,
  contrast.matrix = "pairwise",
  moderated = TRUE, # do moderated t test - if TRUE code gives an error
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

# nothing is significant...
sum(data.res$adj.pvalue <= 0.05, na.rm = TRUE) # 69
data.res %>% filter(adj.pvalue <= 0.05)
