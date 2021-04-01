# Machine learning with random forest in R
# MDCS42270 Bioinformatics analysis of high throughput data
# Lecturer: Dr Carolina Correia
# 01 April 2021

# Exercises were extracted from the DaMiRseq R package manual (Section 3)
# https://www.bioconductor.org/packages/release/bioc/vignettes/DaMiRseq/inst/doc/DaMiRseq.pdf

#### 01 Setup ####

# Required packages were installed when running the prep-data-ML.R script

# Load libraries
library(HDF5Array)
library(DaMiRseq)

# Load summarized experiment with adjusted data
data_adjust <- loadHDF5SummarizedExperiment(dir="data_adjust_se", prefix="")

#### 02 Look at the dataset ####

# Have a look at the adjusted count data
assay(data_adjust)

# Have a look at a subset of the samples associated with the adjusted count data
colData(data_adjust)

#### 03 Feature selection ####

# Let's find a small set of genes (features) that can represent our data
# We currently have 19,343 genes

# Transpose normalised and adjusted expression data, as well as removing
# unwanted characters in feature labels (gene symbols in this case)
set.seed(12345) # This function ensures that randomly created objects can be reproduced
data_clean <- DaMiR.transpose(as.matrix(assay(data_adjust)))

# Variable selection in Partial Least Squares (PLS)
# Here we exclude all non-informative class-related genes
df <- colData(data_adjust)
data_reduced <- DaMiR.FSelect(data_clean, df, th.corr = 0.4)

# Remove highly correlated genes to prevent the inclusion of redundant features that may
# decrease the model performance during the classification step
data_reduced <- DaMiR.FReduct(data_reduced$data)

# Visualise how the most informative samples, obtained after feature selection,
# are clustered in a MultiDimentional Scaling (MDS) plot
DaMiR.MDSplot(data_reduced, df)

# 228 genes is a good reduction from 19,343, but we can go further and
# rank them by importance
df.importance <- DaMiR.FSort(data_reduced, df)

# After the importance score is calculated, a smaller subset of features can be
# selected and used as predictors for classification purpose
# Here, we selected the first 5 genes (default) ranked by importance, check
# the manual to find more information on how to set the number of selected
# features to 'automatic'
selected_features <- DaMiR.FBest(data_reduced, ranking = df.importance,
                                 n.pred = 5)
# Let's visualise the expression and sample data for the 5 selected genes
DaMiR.Clustplot(selected_features$data, df)

#### 04 Classification ####

# All the steps executed so far allowed the reduction of the original expression
# matrix; the objective was to capture a subset of original data as informative
# as possible, in order to carry out a classification analysis

# Let's perform a Bootstrap re-sampling strategy with 30 iterations (iter).
# Ideally we would hundreds of iterations but this means a longer time for
# the code to run.
Classification_res <- DaMiR.EnsembleLearning(selected_features$data,
                                             classes = df$class, fSample.tr = 0.5,
                                             fSample.tr.w = 0.5, iter = 30)

# The violin plot that was automatically generated after classification shows that
# the five selected features ensured moderate and reproducible performance,
# whatever the classifier applied.

# Building a machine learning classification model is only the beginning,
# a lot of effort needs to be put into fine tuning it so that you can obtain
# the optimal prediction model. We are not going to model optimisation
# today, but you can find more information in the package manual.




