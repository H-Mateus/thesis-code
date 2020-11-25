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
library("RforProteomics")
library("mzR")
library("IPPD")
library("biomaRt")
library("MSGFplus")
library("MSnbase")
library("purrr")
library("BiocParallel")

# main guide used:
# https://bioconductor.org/help/course-materials/2017/CSAMA/labs/4-thursday/lab-04-Mass_spec_proteomics_and_metabolomics/01-proteomics/lab.html


# save file path
fl = "/run/media/mateus/4TB-Seagate/proteomic_mass_spec_data/run1_mzml_files/Mateus_iTRAQ_run1_Fr03_4ul.mzML"
# get list of all mzML files - the full names argument includes the file path instead of just the file name
mzMl_files <- list.files("/run/media/mateus/4TB-Seagate/proteomic_mass_spec_data/run1_mzml_files", pattern = "mzML", full.names = TRUE)
mzMl_files
mzMl_files_run2 <- list.files("/run/media/mateus/4TB-Seagate/proteomic_mass_spec_data/run2_mzml_files", pattern = "mzML", full.names = TRUE)
mzMl_files_run2




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

mzid_file = "/run/media/mateus/4TB-Seagate/proteomic_mass_spec_data/run1_mzml_files/Mateus_iTRAQ_run1_Fr03_4ul.mzid"

iddf <- readMzIdData(mzid_file)
names(iddf)

# get identification files for all raw data
#idres_list <- map(mzMl_files, ~runMSGF(msgfpar, .x, memory = 15000))
#idres_list_run2 <- map(mzMl_files_run2, ~runMSGF(msgfpar, .x, memory = 15000))

# many lines of the following are printed when running the search:
# Ignoring spectrum sample=1 period=1 cycle=4777 experiment=3: spectrum is not centroided.

## identification file 
mzid_file <- mzID::files(idres)$id
mzid_file









##### Test with our data #####
# Read in data - This takes a minute or so to run
raw <- readMSData(fl, mode = "onDisk")
raw
raw_try <- raw


# read a single spectrum
sp <- raw[[4]]
sp
# view details of specta
peaksCount(sp)
head(peaksCount(raw))
# retention time
rtime(sp)
# get the mz values and the intensity
mz(sp)
intensity(sp)
peaks(raw)
# use mzid file to identify 
msexp <- addIdentificationData(raw, mzid_file)
fvarLabels(msexp)

## access featureData
fData(msexp)
idSummary(msexp)


# do the same on all run1 data
run1_mzml_data <- map(mzMl_files, ~readMSData(.x, mode = "onDisk"))
run2_mzml_data <- map(mzMl_files_run2, ~readMSData(.x, mode = "onDisk"))

# m/z delta plot
mzp <- plotMzDelta(run1_mzml_data[[1]], reporters = iTRAQ4, verbose = FALSE) + 
  ggtitle("")

mzp

mzp <- plotMzDelta(msexp, reporters = iTRAQ4, verbose = FALSE) + 
  ggtitle("")


# Exercise: Are the spectra in the msexp object centroided? This time, instead of
# looking at the spectra, we can query (and set) this information directly using
# the centroided function (and centroided<- replacement function). The
# isCentroided function can be used to assess the mode from the data.
centroided(msexp)
isCentroided(msexp)
centroided(msexp) <- isCentroided(msexp)


# Spectra and (parts of) experiments can be extracted and plotted.
msexp[[1]]
plot(msexp[[1]], full=TRUE)
# As this data was labeled with iTRAQ4 isobaric tags, we can highlight these four peaks of interest with
plot(msexp[[1]], full=TRUE, reporters = iTRAQ4)

msexp[1:5]
plot(msexp[1:5], full=TRUE)

# I think this adds the iTRAQ labels to the data? - It gives an error but it
# does run anyway, though it takes a minute
test <- quantify(raw, reporters = iTRAQ4, method = "max")
test
# this seems important. Not sure how to use it though 
?addIdentificationData
# is this the spectra by label?
head(exprs(test))

# There are NAs in the data - see how the nrow drops
nrow(test)
test2 <- filterNA(test)
nrow(test2)
# this seems to show how many rows have been removed and some other meta-data
processingData(test2)
head(fData(test2))
protein_translation <- combineFeatures(test2)






# example below - seems like you need a indentification file
## find path to a mzXML file
quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzXML$")
## find path to a mzIdentML file
identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "dummyiTRAQ.mzid")

## create basic MSnExp
msexp <- readMSData(quantFile)

## add identification information
msexp <- addIdentificationData(msexp, identFile)

## access featureData
fData(msexp)
idSummary(msexp)

test_plot <- plotMzDelta(raw_try, reporters = iTRAQ4, verbose = FALSE) + ggtitle("")
test_plot




# view mata-data
head(fData(raw))

# Reporter ions are defined with the Repoterions class - use ?ReporterIons to
# see options and load the relevent class
?ReporterIons
data("iTRAQ4")

reporterNames(iTRAQ4)
sampleNames(raw)


# This took all 32GB of RAM on my laptop and nearly
# filled the SWAP as well. It also didn't finish after running for over 10
# minutes
# plot(raw, full = TRUE) 
plot(test_data, full = TRUE)
plot(test_data[[2]], full = TRUE, reporters = iTRAQ4)

plot(raw[[3]], full = TRUE, reporters = iTRAQ4)

# read chromatographic and/or transition data - again this works with multiple
# files
raw_chrom <- readSRMData(fl)


























## reads the raw data into and MSnExp instance - Note: readMSData can take as
## many files as you like as once whereas openMSfile can only take one at a time
## - the mode = "onDisk" is very important as it determins how the data is
## accessed. The "onDisk" argument prevents the file from being loaded into
## memory

## uses a simple dummy test included in the package
mzXML <- dir(system.file(package="MSnbase",dir="extdata"),
             full.name=TRUE,
             pattern="mzXML$")
basename(mzXML)

test_data <- readMSData(mzXML, verbose = FALSE, centroided = TRUE)
test_data


## Another example dataset:
## Experiment information
library("rpx")
px1 <- PXDataset("PXD000001")
px1
pxfiles(px1)
## Downloading the mzTab data
mztab <- pxget(px1, "PXD000001_mztab.txt")
mztab
### NOTE: this data doesn't appear to be raw - i.e. it already has protein labels included
## Load mzTab peptide data
qnt <- readMzTabData(mztab, what = "PEP", version = "0.9")
sampleNames(qnt) <- reporterNames(TMT6)
head(exprs(qnt))
## remove missing values
qnt <- filterNA(qnt)
processingData(qnt)

## This is an example with raw mass spec data
mzxml <- pxget(px1, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML")
rawms <- readMSData(mzxml, centroided = TRUE, verbose = FALSE)
qntms <- quantify(rawms, reporters = TMT7, method = "max")
qntms
exprs(qntms)

# I think below is taking labels from the earlier labeled example
d <- data.frame(Signal = rowSums(exprs(qntms)[, 1:6]),
                Incomplete = exprs(qntms)[, 7])
d <- log(d)
cls <- rep("#00000050", nrow(qnt))
pch <- rep(1, nrow(qnt))
cls[grep("P02769", fData(qnt)$accession)] <- "gold4" ## BSA
cls[grep("P00924", fData(qnt)$accession)] <- "dodgerblue" ## ENO
cls[grep("P62894", fData(qnt)$accession)] <- "springgreen4" ## CYT
cls[grep("P00489", fData(qnt)$accession)] <- "darkorchid2" ## PHO
pch[grep("P02769", fData(qnt)$accession)] <- 19
pch[grep("P00924", fData(qnt)$accession)] <- 19
pch[grep("P62894", fData(qnt)$accession)] <- 19
pch[grep("P00489", fData(qnt)$accession)] <- 19

plotMzDelta(rawms, reporters = TMT6, verbose = FALSE) + ggtitle("")

plot(Signal ~ Incomplete, data = d,
     xlab = expression(Incomplete~dissociation),
     ylab = expression(Sum~of~reporters~intensities),
     pch = 19,
     col = "#4582B380")
grid()
abline(0, 1, lty = "dotted")
abline(lm(Signal ~ Incomplete, data = d), col = "darkblue")

MAplot(qnt[, c(4, 2)], cex = .9, col = cls, pch = pch, show.statistics = TRUE)


