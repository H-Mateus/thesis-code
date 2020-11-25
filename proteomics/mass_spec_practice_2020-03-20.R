# install packages
# devtools::install_github("sneumann/mzR") # this installs the newest version of
# the package which solves the cvParams error - note you need to restart the R
# session after installing this if you already loaded the normal version of mzR

# BiocManager::install("mzR") # This requiers package: libnetcdf-dev in the
# terminal on Debian 10 - needs ncdf package in r - install netcdf in arch
# BiocManager::install("MSnbase") # note this need gcc-fortran package on arch (use pacman)
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

## Below is use of mzR package to read in files and examine them 
## creates a commection to the mzML file
mz <- openMSfile(fl)
## demonstraction of data access
mz

# get instrament info 
instrumentInfo(mz)

# look at raw spectra by passing the file handle and the index of the spectra of interest
sp1_file2 <- spectra(mz, 1)
# this returns a matrix where the first column is the m/z value and the second column is the intensity of the m/z value
head(sp1_file2)

# plot spectra
sp1_file2 = as.data.frame(sp1_file2)
ggplot(sp1_file2, aes(x = sp1_file2[,1], y = sp1_file2[,2])) +
  geom_point()

# you can pull multiple spectra as a list like so
spl_file2 <- spectra(mz, 1:10)

class(spl_file2)
map(spl_file2, ~head(.x))

# to get the annotations for the individual spectra - this can take some time as
# it creates a dataframe with as many rows as there are spectra
hd <- header(mz)
head(hd)

# once finished, it is good to explicitely close the connection
close(mz)