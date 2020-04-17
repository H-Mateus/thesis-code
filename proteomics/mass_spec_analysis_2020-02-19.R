# install packages
# devtools::install_github("sneumann/mzR") # this installs the newest version of
      # the package which solves the cvParams error - note you need to restart the R
      # session after installing this if you already loaded the normal version of mzR

# BiocManager::install("mzR") # This requiers package: libnetcdf-dev in the
      # terminal on Debian 10 - needs ncdf package in r - install netcdf in arch
# BiocManager::install("MSnb  ase") # note this need gcc-fortran package on arch (use pacman)
# BiocManager::install("RforProteomics")
# BiocManager::install("MSGFplus") # some sort of java error
# devtools::install_github('thomasp85/MSGFplus') # this works on Debain 10 and
      # arch - note it needs java to be installed to work

# BiocManager::install("BiocParallel")
# BiocManager::install("IPPD")
# BiocManager::install("biomaRt")

# load packages 
library("mzR")
library("mzID")
library("MSnID")
library("MSnbase")
library("rpx")
library("MLInterfaces")
library("pRoloc")
library("pRolocdata")
library("MSGFplus")
library("rols")
library("hpar")
library("RforProteomics")
library("purrr")

# handling raw ms data

# inspect file
# save file path
fl = "/run/media/mateus/4TB-Seagate/proteomic_mass_spec_data/wiff_files/mzXML_files/Mateus_iTRAQ_run1_Fr01_4ul.mzXML"

ms <- openMSfile(fl)
ms
instrumentInfo(ms)
runInfo(ms) # this takes a whlie to run 

# extract metadata for scan 1000
hd <- header(ms)
dim(hd)
names(hd)
head(hd)

head(peaks(ms, 1000), 20)
plot(peaks(ms, 1000), type = "h")
# save plot - note this code will save the currently displayed plot
dev.print(pdf, 'ms_spectra_plot_2020-03-31.pdf')

## a set of spectra of interest: MS1 spectra eluted
## between 30 and 35 minutes retention time
ms1 <- which(hd$msLevel == 1)
rtsel <- hd$retentionTime[ms1] / 60 > 30 &
  hd$retentionTime[ms1] / 60 < 35

## the map
M <- MSmap(ms, ms1[rtsel], 521, 523, .005, hd)
plot(M, aspect = 1, allTicks = FALSE)
#dev.print(pdf, 'ms1_spectra_eluted_between_30-35mins_2020-03-31.pdf')

plot3D(M)
#dev.print(pdf, 'ms1_spectra_eluted_between_30-35mins_3d_2020-03-31.pdf')

## With some MS2 spectra
i <- ms1[which(rtsel)][1]
j <- ms1[which(rtsel)][2]
M2 <- MSmap(ms, i:j, 100, 1000, 1, hd)
plot3D(M2)
#dev.print(pdf, 'ms2_spectra_eluted_between_30-35mins_3d_2020-03-31.pdf')



### Handling indentification
# MSGF+ engine database search to generate mzid file - needs a fasta file of
# proteins - I've used uniprots fasta of human proteins
# set up search paramaters
msgfpar <- msgfPar(database = "../proteomics/data/human_uniprot_UP000005640_9606.fasta",
                   instrument = 'HighRes',
                   tda = TRUE,
                   enzyme = 'Trypsin',
                   protocol = 'iTRAQ')

# generate the mzid file - NOTE: this will create the mzid file and save it in
# the location of the mzML file. You can set the file name with savenames =
# "file_name", but it defaults to the same file name but with the .mzid file
# extention

# 2 lines below is running a single file 
#idres <- runMSGF(msgfpar, fl, memory=15000) 
#mzid_file <- mzID::files(idres)$id
#mzid_file


# filter mzid files 
mzid_files <- list.files("/run/media/mateus/4TB-Seagate/proteomic_mass_spec_data/wiff_files/mzXML_files", pattern = "mzid", full.names = TRUE)
mzid_files
# get run1 files
mzid_files_run1 <- mzid_files[1:12]


msnid <- MSnID(".")
# note: the following function can take multiple files at once
msnid <- read_mzIDs(msnid, mzid_files_run1)
# note: the FDR seems really high here
show(msnid)

msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")

## Trimming the data
# Now, we can use the apply_filter function to effectively apply filters. The
# strings passed to the function represent expressions that will be evaludated,
# this keeping only PSMs that have 0 irregular cleavages and 2 or less missed
# cleavages.
msnid <- apply_filter(msnid, "numIrregCleavages == 0")
msnid <- apply_filter(msnid, "numMissCleavages <= 2")
show(msnid)

# parent ion mass errors
# Using "calculatedMassToCharge" and "experimentalMassToCharge", the
# mass_measurement_error function calculates the parent ion mass measurement
# error in parts per million.
summary(mass_measurement_error(msnid))
# We then filter any matches that do not fit the +/- 20 ppm tolerance
msnid <- apply_filter(msnid, "abs(mass_measurement_error(msnid)) < 20")
summary(mass_measurement_error(msnid))
show(msnid)

## filtering criteria 
#  Filtering of the identification data will rely on: 
# - log10 transformed MS-GF+ Spectrum E-value, reflecting the goodness of match
# experimental and theoretical fragmentation patterns
msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
# - the absolute mass measurement error (in ppm units) of the parent ion
msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))
# MS2 fiters are handled by a special MSnIDFilter class objects, where
# individual filters are set by name (that is present in names(msnid)) and
# comparison operator (>, <, = , …) defining if we should retain hits with
# higher or lower given the threshold and finally the threshold value itself.
filtObj <- MSnIDFilter(msnid)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
filtObj$msmsScore <- list(comparison=">", threshold=10.0)
show(filtObj)
# We can then evaluate the filter on the identification data object, which
# return the false discovery rate and number of retained identifications for the
# filtering criteria at hand.
evaluate_filter(msnid, filtObj)

filtObj <- MSnIDFilter(msnid)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
filtObj$msmsScore <- list(comparison=">", threshold=10.0)
show(filtObj)

filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01,
                                method="Grid", level="peptide",
                                n.iter=500)

show(filtObj.grid)
evaluate_filter(msnid, filtObj.grid) # NOTE: this results in nothing left... Maybe one of the files is parcitularly bad?

# Filters can eventually be applied (rather than just evaluated) using the apply_filter function.
msnid_filtered <- apply_filter(msnid, filtObj)
show(msnid_filtered)
msnid_filtered <- apply_filter(msnid_filtered, "isDecoy == FALSE")
msnid_filtered <- apply_filter(msnid_filtered, "!grepl('Contaminant',accession)")
show(msnid_filtered) # 446 peptides left for run1 but there are duplicates

# peptide sequence - note: there are duplicates for some reason. Maybe multiple
# scans led to the same protein match?
unique(msnid_filtered$pepSeq)
# identified proteins
unique(msnid_filtered$accession)

### Get protein names and function from msnid
# code to extract string between two characters 
try <- unique(msnid_filtered$accession)
try
library(qdapRegex)
try <- rm_between(try, "|", "|", extract = TRUE)
protein_codes <- unlist(try)
protein_codes


# get protein name and description from uniprot
library(purrr)
# source uniprot search functions
source("../current-thesis-template/thesis/scripts/functions.R")
protein_names <- map(protein_codes, ~find_protein_name(.x))
protein_names
protein_descriptions <- map(protein_codes, ~find_protein_description(.x))
protein_descriptions
protein_df <- data.frame(name = unlist(protein_names), description = unlist(protein_descriptions))

# filter out keratin
protein_df_filtered <- dplyr::filter(protein_df, !grepl("Keratin", name))
# save dataframe 
save(protein_df_filtered, file = "../proteomics/data/mass_spec_protein.Rda")

