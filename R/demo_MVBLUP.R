# demo_MVBLUP.R
# The script was used to build the MVBLUP model, obtain optimal dataset weights, and test the accuracy of the MVBLUP model.

# to run this MVBLUP algrithm, nine files are required to be prepared as the format of the demo data:
# 1. Multi_view1.txt: demonstrative data from the 1st view, consisting of 100 individuals and 3000 characteristics, requires transformation into numeric values. The individual identifier needs to match the Phenotype data.
# 2. Multi_view2.txt: demonstrative data from the 2nd view, consisting of 100 individuals and 1000 characteristics, requires transformation into numeric values. The individual identifier needs to match the Phenotype data.
# 3. Multi_view3.txt: demonstrative data from the 3rd view, consisting of 100 individuals and 1000 characteristics, requires transformation into numeric values. The individual identifier needs to match the Phenotype data.
# 4. Multi_view4.txt: demonstrative data from the 4th view, consisting of 100 individuals and 1000 characteristics, requires transformation into numeric values. The individual identifier needs to match the Phenotype data.
# 5. Multi_view5.txt: demonstrative data from the 5th view, consisting of 100 individuals and 1000 characteristics, requires transformation into numeric values. The individual identifier needs to match the Phenotype data.
# 6. Multi_view6.txt: demonstrative data from the 6th view, consisting of 100 individuals and 1000 characteristics, requires transformation into numeric values. The individual identifier needs to match the Phenotype data.
# 7. Multi_view7.txt: demonstrative data from the 7th view, consisting of 100 individuals and 1000 characteristics, requires transformation into numeric values. The individual identifier needs to match the Phenotype data.
# 8. Multi_view8.txt: demonstrative data from the 8th view, consisting of 100 individuals and 1000 characteristics, requires transformation into numeric values. The individual identifier needs to match the Phenotype data.
# 9. Phenotype.txt: demonstrative phenotypic data with 100 individuals and 2 traits. 
# 10. train_test_id.txt: an information file with 100 rows and 2 columns. The first column contains individual identifiers, and the second column indicates the use of the training or test set.

# Several parameters are demonstrated as followed:

# NP: the size of population
# n: the number of multi-view datasets
# CR: crossover factor
# Mu: scaling factor
# s0: the lower bounds of the search space for weights
# s1: the upper bounds of the search space for weights
# thre: threshold of the early stopping mechanism
# iter: maximum iteration
# cv: the number of folds for cross validation
# trid: which trait
# output_path: directory for storing the output results
rm(list=ls())
setwd("D:/maize/MVBLUP/")
# load the script to obtain the kinship matrices
source("./R/similarity_matrix.R")
# load the script to obtain the trait predictions of training set and test set
source("./R/GP.R")
# load the script to build the MVBLUP model and obtain the accuracy of predicting the test set
source("./R/MVBLUP.R")
# load the script to visualize the learning process of the MVBLUP model
source("./R/plot.R")

library(rrBLUP)
library(ggplot2)
library(RColorBrewer)

NP <- 40 
n <- 8
CR <- 0.5
Mu <- 0.5
s0 <- 0
s1 <- 1
thre <- 0.0001
iter <- 50
cv <- 5
trid <- 1
output_path <- "./res/"

# read the multi-view dataset
Multi_view1 <- read.table("./data/Multi_view1.txt", header=TRUE, sep="", stringsAsFactors=0, check.names=FALSE)
Multi_view2 <- read.table("./data/Multi_view2.txt", header=TRUE, sep="", stringsAsFactors=0, check.names=FALSE)
Multi_view3 <- read.table("./data/Multi_view3.txt", header=TRUE, sep="", stringsAsFactors=0, check.names=FALSE)
Multi_view4 <- read.table("./data/Multi_view4.txt", header=TRUE, sep="", stringsAsFactors=0, check.names=FALSE)
Multi_view5 <- read.table("./data/Multi_view5.txt", header=TRUE, sep="", stringsAsFactors=0, check.names=FALSE)
Multi_view6 <- read.table("./data/Multi_view6.txt", header=TRUE, sep="", stringsAsFactors=0, check.names=FALSE)
Multi_view7 <- read.table("./data/Multi_view7.txt", header=TRUE, sep="", stringsAsFactors=0, check.names=FALSE)
Multi_view8 <- read.table("./data/Multi_view8.txt", header=TRUE, sep="", stringsAsFactors=0, check.names=FALSE)
Multi_dataset <- list(Multi_view1, Multi_view2, Multi_view3, Multi_view4, Multi_view5, Multi_view6, Multi_view7, Multi_view8)

# read phenotype
Phenotype <- read.table("./data/Phenotype.txt", header = TRUE, sep = "", stringsAsFactors = 0, check.names = FALSE)
train_test_id <- read.table("./data/train_test_id.txt", header = TRUE, sep = "", stringsAsFactors = 0, check.names = FALSE)

# obtain single-view inner product for calculating multi-view kinship matrix
A_N <- similarity_inner_product(Multi_dataset, n)

# optimize the weights of various single-view datasets by employing the MVBLUP model
MVBLUP_results <- MVBLUP_weights(NP, n, CR, Mu, s0, s1, thre, iter, cv, trid, A_N,
                                 Phenotype, train_test_id)

# multi-view kinship matrix obtained from the MVBLUP model
kinship_test <- MVBLUP_results$kinship_test
# testop1_ex: true phenotypes of the test set
testop1_ex <- MVBLUP_results$testop1_ex
# trainop1_ex: true phenotypes of the train set
trainop1_ex <- MVBLUP_results$trainop1_ex

# We output a demo result of MVBLUP model, which includes the optimal weight values for different single-view datasets and the number of iterations required for convergence.
write.table(MVBLUP_results$MVBLUP_information, paste0(output_path, "trait", trid, "-MVBLUP_information.txt"), row.names = FALSE)

# test the MVBLUP model and store the results
MVBLUP_accuracy <- MVBLUP_PreTest(kinship_test, testop1_ex, trainop1_ex)
write.table(MVBLUP_accuracy, paste0(output_path, "trait", trid, "-MVBLUP_accuracy.txt"), row.names = TRUE)

# Visualize the learning process of the MVBLUP model and save it to the output directory
plot_learning(MVBLUP_results$Training_Accuracy,
              MVBLUP_results$Training_Weight,
              MVBLUP_results$MVBLUP_information,
              output_path)

# MVBLUP_results$Training_Accuracy: the average accuracy obtained from each iteration of five-fold cross-validation on the training set.
# MVBLUP_results$Training_Weight: the optimal weights obtained in each iteration.
# MVBLUP_results$MVBLUP_information: the optimal weight values and the number of iterations required for convergence
# output_path: the path for storing the output results.
