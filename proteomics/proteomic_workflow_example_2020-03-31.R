# example from: https://master.bioconductor.org/packages/release/workflows/vignettes/proteomics/inst/doc/proteomics.html
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

# explore aviable proteomic packages 
pp <- proteomicsPackages()
pp


# handling raw ms data
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
dev.print(pdf, 'ms1_spectra_eluted_between_30-35mins_2020-03-31.pdf')

plot3D(M)
dev.print(pdf, 'ms1_spectra_eluted_between_30-35mins_3d_2020-03-31.pdf')

## With some MS2 spectra
i <- ms1[which(rtsel)][1]
j <- ms1[which(rtsel)][2]
M2 <- MSmap(ms, i:j, 100, 1000, 1, hd)
plot3D(M2)
dev.print(pdf, 'ms2_spectra_eluted_between_30-35mins_3d_2020-03-31.pdf')


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
idres <- runMSGF(msgfpar, fl, memory=15000) 
mzid_file <- mzID::files(idres)$id
mzid_file

# file path to mzid file
mzid_file = "/run/media/mateus/4TB-Seagate/proteomic_mass_spec_data/wiff_files/mzXML_files/Mateus_iTRAQ_run1_Fr01_4ul.mzid"


  
# Various data can be extracted from the mzID object, using one the accessor
# functions such as database, scans, peptides, … The object can also be
# converted into a data.frame using the flatten function.
id <- mzID(mzid_file)
id
fid <- flatten(id)
names(fid)
dim(fid)
head(fid)
# The mzR package also provides support fasta parsing mzIdentML files with the
# openIDfile function. As for raw data, the underlying C/C++ code comes from the
# proteowizard.

id1 <- openIDfile(mzid_file)
fid1 <- mzR::psms(id1)

head(fid1)
# meta-data in mzid file
softwareInfo(id1)
enzymes(id1)
names(psms(id1))
head(psms(id1))[,1:13]
# below is some stuff about manual filtering of data - Don't currently
# understand the data well enough to use do things this way - see below for a
# different high-level interface

# The MSnID package can be used for post-search filtering of MS/MS
# identifications. One starts with the construction of an MSnID object that is
# populated with identification results that can be imported from a data.frame
# or from mzIdenML files. Here, we will use the example identification data
# provided with the package.

# We start by loading the package, initialising the MSnID object, and add the
# identification result from our mzid file (there could of course be more that
# one).
library("MSnID")
msnid <- MSnID(".")

msnid <- read_mzIDs(msnid, mzid_file)
# note: the FDR seems really high here - 26499 peptides here
show(msnid)

# Printing the MSnID object returns some basic information such as:
# Working directory.
# Number of spectrum files used to generate data.
# Number of peptide-to-spectrum matches and corresponding FDR.
# Number of unique peptide sequences and corresponding FDR.
# Number of unique proteins or amino acid sequence accessions and corresponding FDR.

names(msnid)

## analysis of peptide sequences
# Cleaning irregular cleavages at the termini of the peptides and missing
# cleavage site within the peptide sequences. The following two function call
# create the new numMisCleavages and numIrrCleabages columns in the MSnID object
msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")

## Trimming the data
# Now, we can use the apply_filter function to effectively apply filters. The
# strings passed to the function represent expressions that will be evaludated,
# this keeping only PSMs that have 0 irregular cleavages and 2 or less missed
# cleavages.
msnid <- apply_filter(msnid, "numIrregCleavages == 0")
msnid <- apply_filter(msnid, "numMissCleavages <= 2")
show(msnid) # 20601 peptides now

# parent ion mass errors
# Using "calculatedMassToCharge" and "experimentalMassToCharge", the
# mass_measurement_error function calculates the parent ion mass measurement
# error in parts per million.
summary(mass_measurement_error(msnid))
# We then filter any matches that do not fit the +/- 20 ppm tolerance
msnid <- apply_filter(msnid, "abs(mass_measurement_error(msnid)) < 20")
summary(mass_measurement_error(msnid))
show(msnid) # 10822 peptides now

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

## Filter optimisation
# Rather than setting filtering values by hand, as shown above, these can be set
# automativally to meet a specific false discovery rate.
filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01,
                                method="Grid", level="peptide",
                                n.iter=500)
show(filtObj.grid)
evaluate_filter(msnid, filtObj.grid)
# NOTE: The number after filtering is 91. This seems like a massive reduction
# from the thousands before filtering. Considering there were only around 65
# proteins left after cleaning up the ProteinPilot data I suppose this is a
# similar output. But why is the FDR so high as to require over 90% of the data
# to be filtered out?

# Filters can eventually be applied (rather than just evaluated) using the apply_filter function.
msnid <- apply_filter(msnid, filtObj.grid)
show(msnid)
msnid <- apply_filter(msnid, "isDecoy == FALSE")
msnid <- apply_filter(msnid, "!grepl('Contaminant',accession)")
show(msnid) # 52 peptides left

# peptide sequence - note: there are duplicates for some reason. Maybe multiple
# scans led to the same protein match?
uniques(msnid$pepSeq)
# identified proteins
unique(msnid$accession)


### Get protein names and function from msnid
# code to extract string between two characters 
try <- unique(msnid$accession)
try
library(qdapRegex)
try <- rm_between(try, "|", "|", extract = TRUE)
protein_codes <- unlist(try)
protein_codes

# code to extract gene
library(stringr)
gene_try <- as.character(unique(msnid$accession))
# select all characters after 8th
gene_names <- str_sub(unlist(gene_try), 8)
gene_names

# get protein name and description from uniprot
library(purrr)
# source uniprot search functions
source("../current-thesis-template/thesis/scripts/functions.R")
protein_names <- map(protein_codes, ~find_protein_name(.x))
protein_names
protein_descriptions <- map(protein_codes, ~find_protein_description(.x))
protein_descriptions
protein_df <- data.frame(name = unlist(protein_names), description = unlist(protein_descriptions))

#### High-level data interface
# The above sections introduced low-level interfaces to raw and identification
# results. The MSnbase package provides abstractions for raw data through the
# MSnExp class and containers for quantification data via the MSnSet class. Both
# store:
# - the actual assay data (spectra or quantitation matrix, see below), accessed with spectra (or the [, [[ operators) or exprs;
# - sample metadata, accessed as a data.frame with pData;
# - feature metadata, accessed as a data.frame with fData.

# Another useful slot is processingData, accessed with processingData(.), that
# records all the processing that objects have undergone since their creation
# (see examples below).

# The readMSData will parse the raw data, extract the MS2 spectra (by default)
# and construct an MS experiment object of class MSnExp.
# (Note that while readMSData supports MS1 data, this is currently not
# convenient as all the data is read into memory.)
#NOTE: if you don't specify the mode it extracts the MS2 spectra by default - it also takes a lot longer to run (~a few minutes)

msexp <- readMSData(fl, verbose = FALSE, msLevel = 2)
msexp
# MS2 spectra can be extracted as a list of Spectrum2 objects with the spectra
# accessor or as a subset of the original MSnExp data with the [ operator.
# Individual spectra can be accessed with [[.
length(msexp)
msexp[1:2]
msexp[[2]]
# The identification results stemming from the same raw data file can then be
# used to add PSM matches.
head(fData(msexp))
head(pData(msexp))

# may need mslevel 1 to get chormatogram objects
mtof_tic <- chromatogram(msexp)
mtof_tic


iTRAQ4
mz(iTRAQ4)
width(iTRAQ4)
# Typically, identification data is produced by a search engine and serialised
# to disk in the mzIdentML (or mzid) file format. This format can be parsed by
# openIDfile from the mzR package or mzID from the mzID package. The MSnbase
# package relies on the former (which is faster) and offers a simplified
# interface by converting the dedicated identification data objects into
# data.frames.
iddf <- readMzIdData(mzid_file)
str(iddf)

table(iddf$isDecoy)
table(iddf$chargeState)
library("ggplot2")
ggplot(data = iddf, aes(x = MS.GF.RawScore, colour = isDecoy)) +
  geom_density() +
  facet_wrap(~chargeState)

# filtering mzid with MSbase package

# The filterIdentificationDataFrame function can be used to remove - PSMs that
# match decoy entries - PSMs of rank > 1 - PSMs that match non-proteotypic
# proteins
iddf <- filterIdentificationDataFrame(iddf)
nrow(iddf)
# This data.frame can be now be further reduced so that individual rows
# represent unique spectra, which can be done with the reduce method.
iddf2 <- IRanges::reduce(iddf, key = "spectrumID")
nrow(iddf2)
# This reduces the number of rows from 36263 to 12312.

# The first duplicated spectrum mentioned above is now unique as is matched a
# decoy protein that was filtered out with filterIdentificationDataFrame.




# save unmodififed data
msexp_save <- msexp

# the addIdentificationData function gives an error about argument "decoy" being
# missing when using the iddf dataframe
#msexp <- addIdentificationData(msexp, iddf2)

msexp <- addIdentificationData(msexp, mzid_file)

# try adding the flattened mzid file that contains the accession data
# msexp_test <- addIdentificationData(msexp, fid)

head(fData(msexp))
idSummary(msexp)

itraqdata2 <- pickPeaks(msexp, verbose=FALSE)
i <- 16 # note number needs to be rows where sequence col isn't NA
s <- as.character(fData(itraqdata2)[i, "sequence"])
plot(itraqdata2[[i]], s, main = s)

# The readMSData and addIdentificationData make use of mzR and mzID packages to
# access the raw and identification data.

# Spectra and (parts of) experiments can be extracted and plotted.
msexp[[1]]
plot(msexp[[1]], full=TRUE)
msexp[1:3]
plot(msexp[1:3], full=TRUE)
dev.print(pdf, 'mzxml_spectra_2020-04-01.pdf')

### Filtering identification data
# One can remove the features that have not been identified using removeNoId.
# This function uses by default the pepseq feature variable to search the
# presence of missing data (NA values) and then filter these non-identified
# spectra.
idSummary(msexp)
head(fData(msexp)$sequence)
msexp <- removeNoId(msexp)
head(fData(msexp)$sequence)
idSummary(msexp) # the mzid file only appears in this output when the NAs have been removed

# Similarly, the removeMultipleAssignment method can be used to filter out
# non-unique features, i.e. that have been assigned to protein groups with more
# than one member. This function uses by default the nprot feature variable.
# Note: removeNoId and removeMultipleAssignment methods can also be called
# on MSnExp instances.
head(fData(msexp)$sequence, 10)
msexp <- removeMultipleAssignment(msexp)
head(fData(msexp)$sequence, 10)
msexp

## calculate fragments 
pepseq <- fData(msexp)$sequence[1]
calculateFragments(pepseq, msexp[[1]], type=c("b", "y"))


### Quality control 
# MSnbase allows easy and flexible access to the data, which
# allows to visualise data features to assess it’s quality. Some methods are
# readily available, although many QC approaches will be experiment specific and
# users are encourage to explore their data.

# The plot2d method takes one MSnExp instance as first argument to produce
# retention time vs. precursor MZ scatter plots. Points represent individual MS2
# spectra and can be coloured based on precursor charge (with second argument
# z="charge"), total ion count (z="ionCount"), number of peaks in the MS2
# spectra z="peaks.count") or, when multiple data files were loaded, file
# z="file"), as illustrated on the next figure. The lower right panel is
# produced for only a subset of proteins. See the method documentation for more
# details.
MSnbase::plot2d(msexp, z = "ionCount")
MSnbase::plot2d(msexp, z = "peaks.count")
MSnbase::plot2d(msexp, z = "charge") # only one charge? Maybe need msLevel 1 for charge?
MSnbase::plotDensity(msexp, z = "peaks.count")
MSnbase::plotDensity(msexp, z = "ionCount")
MSnbase::plotDensity(msexp, z = "precursor.mz")
#NOTE: this creates a 1.7Gb plot. Attempting to display the plot resulted in
#crashing R even though the system should have had enough memory
# mzdplot <- MSnbase::plotMzDelta(msexp,
#                      subset = 0.5,
#                      reporters = iTRAQ4,
#                      verbose = FALSE,
#                      plot = FALSE)




## Quantitative proteomics

# There are a wide range of proteomics quantitation techniques that can broadly
# be classified as labelled vs. label-free, depending whether the features are
# labelled prior the MS acquisition and the MS level at which quantitation is
# inferred, namely MS1 or MS2.

# In terms of raw data quantitation, most efforts have been devoted to MS2-level
# quantitation. Label-free XIC quantitation has however been addressed in the
# frame of metabolomics data processing by the xcms infrastructure.

# An MSnExp is converted to an MSnSet by the quantitation method. Below, we use
# the iTRAQ 4-plex isobaric tagging strategy (defined by the iTRAQ4 parameter;
# other tags are available) and the trapezoidation method to calculate the area
# under the isobaric reporter peaks.
plot(msexp[[1]], full=TRUE, reporters = iTRAQ4)
dev.print(pdf, 'mzxml_spectra_itraq_reporters_2020-04-01.pdf')

msset <- quantify(msexp, method = "trap", reporters = iTRAQ4, verbose=FALSE)
# I assume this gives the relative expression of a peptide(or protein?) in the spectra?
head(exprs(msset), 20) # there are some NAs in this.

processingData(msset)
msset
# Other MS2 quantitation methods available in quantify include the (normalised)
# spectral index SI and (normalised) spectral abundance factor SAF or simply a
# simple count method.
exprs(si <- quantify(msexp, method = "SIn")) # gives error 
exprs(saf <- quantify(msexp, method = "NSAF")) # same error 

### Data processing and analysis 
# For raw data processing look at MSnbase’s clean, smooth, pickPeaks,
# removePeaks and trimMz for MSnExp and spectra processing methods.
?MSnbase::clean # this seems to remove 0-intensity peaks - but why are these peaks even there in the first place?
?MSnbase::smooth 
?MSnbase::pickPeaks # creates cenroided spectra?
?MSnbase::removePeaks # removes peaks below a specified intensity
?MSnbase::trimMz
#The MALDIquantand xcms packages also features a wide range of raw data
#processing methods on their own ad hoc data instance types.


## processing and normalisation
impurities <- makeImpuritiesMatrix(4)

qnt.crct <- purityCorrect(msset, impurities)
processingData(qnt.crct)
# Various normalisation methods can be applied the MSnSet instances using the
# normalise method: variance stabilisation (vsn), quantile (quantiles), median
# or mean centring (center.media or center.mean), …
qnt.crct.nrm <- normalise(qnt.crct, "quantiles")
qnt.crct.nrm


# The combineFeatures method combines spectra/peptides quantitation values into
# protein data. The grouping is defined by the groupBy parameter, which is
# generally taken from the feature metadata (protein accessions, for example).
## arbitraty grouping
g <- factor(c(rep(1, 25), rep(2, 15), rep(3, 15)))
g

prt <- combineFeatures(qnt.crct.nrm, groupBy = g, fun = "sum") # this wont work unless the groupBy variable 
                                                               # has the same length as the nrow of the data
processingData(prt)

# Dealing with NAs
image(qnt.crct.nrm)
dim(qnt.crct.nrm)
table(is.na(qnt.crct.nrm))
## remove features with missing values
qnt00 <- filterNA(qnt.crct.nrm)
dim(qnt00)
table(is.na(qnt00))
## impute missing values using knn imputation
qnt.imp <- impute(qnt.crct.nrm, method = "knn")
dim(qnt.imp)
table(is.na(qnt.imp))

head(fData(qnt00))
nrow(fData(qnt00))


### Get protein names and function from msnid
# code to extract string between two characters 

test_proteins <- unique(fData(qnt00)$DatabaseAccess)
test_proteins
library(qdapRegex)
try <- rm_between(test_proteins, "|", "|", extract = TRUE)
protein_codes_test <- unlist(try)
head(protein_codes_test)
# remove the NA
protein_codes_test <- protein_codes_test[2:length(protein_codes_test)]
head(protein_codes_test)

# get protein name and description from uniprot
library(purrr)
# source uniprot search functions
source("../current-thesis-template/thesis/scripts/functions.R")
protein_names_test <- map(protein_codes_test, ~find_protein_name(.x))
head(protein_names_test)
protein_descriptions_test <- map(protein_codes_test, ~find_protein_description(.x))
head(protein_descriptions_test)
protein_df_test <- data.frame(name = unlist(protein_names), description = unlist(protein_descriptions_test))

# filter out keratin
protein_df_filtered_test <- dplyr::filter(protein_df_test, !grepl("Keratin", name))
# save dataframe 
save(protein_df_filtered_test, file = "../proteomics/data/mass_spec_protein_test.Rda")

# to do: check is that filtered msid file could be used with this high-level interface

## NOTES
# the mzid file has the accession meta data in it, but using the
# addIdentificationData function doesn't seem to add it in.





### Machine learning 
# The MLInterfaces package provides a unified interface to a wide range of
# machine learning algorithms. Initially developed for microarray and
# ExpressionSet instances, the pRoloc package enables application of these
# algorithms to MSnSet data.

## Classification 
# The example below uses knn with the 5 closest neighbours as an illustration to
# classify proteins of unknown sub-cellular localisation to one of 9 possible
# organelles.
# NOTE: our data doesn't seem to have the "markers" column
data(dunkley2006)
traininds <- which(fData(dunkley2006)$markers != "unknown")
ans <- MLearn(markers ~ ., data = t(dunkley2006), knnI(k = 5), traininds)
ans

## Clustering 
# kmeans
kcl <- MLearn( ~ ., data = dunkley2006, kmeansI, centers = 12)
kcl
plot(kcl, exprs(dunkley2006))
# A wide range of classification and clustering algorithms are also available,
# as described in the ?MLearn documentation page. The pRoloc package also uses
# MSnSet instances as input and ,while being conceived with the analysis of
# spatial/organelle proteomics data in mind, is applicable many use cases.


library("rols")
res <- OlsSearch(q = "ESI", ontology = "MS", exact = TRUE)
res
