# Machine learning with random forests in R
# MDCS42270 Bioinformatics analysis of high throughput data
# Lecturer: Dr Carolina Correia
# 22nd March 2021

# Exercises were extracted from the DaMiRseq R package manual (Section 3)
# https://www.bioconductor.org/packages/release/bioc/vignettes/DaMiRseq/inst/doc/DaMiRseq.pdf

# The goal of this script is to prepare data to be used during class on 01/04/2021

#### 01 Setup ####

# You only need to install the packages once, then load the libraries when
# it's time to use them
install.packages("BiocManager") # You need BiocManager to install packages from Bioconductor
BiocManager::install(version = "3.12") # We are using Bioconductor version 3.12
# You might be asked the following question:
# Do you want to install from sources the package which needs compilation? (Yes/no/cancel)
# Answer Yes in the console

BiocManager::install("DaMiRseq") # This is how you install a package from Bioconductor
# You might be asked the following question:
# Update all/some/none? [a/s/n]: a
# Answer all in the console

# You might get an error saying that rJava was not installed or had a problem
# For macOS users: go to terminal and run sudo R CMD javareconf
# You will be asked to type your computer's password
# If you encounter problems in Linux or Windows, send me an email with the
# error and I'll find out how to fix it.

BiocManager::install("HDF5Array") # This is how you install a package from Bioconductor

# Load libraries
library(HDF5Array)
library(DaMiRseq)

#### 02 Look at the dataset ####

# Load the raw count data available from DaMiRseq
# SE stands for Summarized Experiment
data(SE)

# Have a look at a subset of the raw count data
assay(SE)[1:5, c(1:5, 21:25)]

# Have a look at a subset of the samples associated with the raw count data
colData(SE)

#### 03 Filter and normalise data ####

# Using this function, DaMiRseq will apply a filtering by expression values and
# coefficient of Variation (CV). A VST normalisation is also applied by default
# It means that genes (aka features) with very low expression and those that
# could be outliers will be removed from the data
data_norm <- DaMiR.normalization(SE, minCounts = 10, fSample = 0.7,
                                 hyper = "yes", th.cv = 3)

# Using assay() function, we can see that the VST transformation produces data on
# the log2 scale normalized with respect to the library size
assay(data_norm)[c(1:5), c(1:5, 21:25)]

# Now, we need to check for any samples that show inconsistent behaviour
data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.9)

#### 04 Adjust data ####

# Let's test for the presence of surrogate variables (sv) in order to remove
# the effect of putative confounding factors from the expression data
sv <- DaMiR.SV(data_filt)

# Check the correlation between sv and known covariates
DaMiR.corrplot(sv, colData(data_filt), sig.level = 0.01)

# After sv identification, we adjust our expression data by removing unwanted
# effects from the expression data matrix
data_adjust<-DaMiR.SVadjust(data_filt, sv, n.sv = 4)

#### 05 Explore data ####

# Generate 5 plots: heatmap, MultiDimensional Scaling (MDS),
# Relative Log Expression (RLE), Sample-by-Sample expression distribution, and
# Average expression distribution by class

# First, generate all 5 plots with the filtered data (before surrogate variable adjustment)
DaMiR.Allplot(data_filt, colData(data_filt))

# Now generate all 5 plots with the adjusted data (after surrogate variable adjustment)
DaMiR.Allplot(data_adjust, colData(data_adjust))

# You can compare all plots before and after sv adjustment and see how the data have improved

#### 06 Save summarized experiment ####

# Using the HDF5Array R package, we can save all the adjusted data in our
# summarized experiment to use during class next week
saveHDF5SummarizedExperiment(data_adjust, dir = "data_adjust_se",
                             prefix = "", replace = FALSE,
                             chunkdim = NULL, level = NULL,
                             as.sparse = NA, verbose = NA)
