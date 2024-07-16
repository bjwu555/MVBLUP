# The script was used to construct the MVBLUP model and determine the accuracy of test set predictions
# Notes:
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
# A_N: the kinship matrix calculated using a custom similarity function
# Phenotype: demonstrative phenotypic data with 100 individuals and 2 traits
# train_test_id: an information file with 100 rows and 2 columns. The first column contains individual identifiers, and the second column indicates the use of the training or test set
# A_test: multi-view kinship matrix obtained from the MVBLUP model
# testop1_ex: true phenotypes of the test set
# trainop1_ex: true phenotypes of the train set

set.seed(121)
#========Learning the optimal weights for datasets from different perspectives based on differential evolution algorithm================
MVBLUP_weights <- function(NP, n, CR, Mu, s0, s1, thre, iter, cv, trid, A_N, Phenotype, train_test_id) {
  
  num_ex <- floor(nrow(Phenotype) / cv)
  index_sample <- match(train_test_id[, 1], Phenotype[, 1])
  dt_ex <- Phenotype[index_sample, -1]
  samplenames <- train_test_id[, 1]
  dt_in <- dt_ex[, trid]
  dt_in[is.na(dt_in)] <- mean(dt_in, na.rm = TRUE)
  trainop1_ex <- dt_in[-c(1:num_ex)]
  testop1_ex <- dt_in[1:num_ex]
  
  AA <- array(0, dim = c(nrow(train_test_id), nrow(train_test_id), n))
  for (n0 in 1:n) {
    AA[,, n0] <- A_N[,, n0][index_sample, index_sample]
  }
  
  Training_Accuracy <- data.frame(iter = rep(1:iter, each = NP), group = rep(1:NP, iter), Accuracy = NA)
  Training_Weight <- data.frame(iter = rep(1:iter, each = n), Parameter = rep(1:n, iter), Weight = NA)
  num_test <- floor(length(trainop1_ex) / cv)
  
  BM <- matrix(runif(n * NP, 0, 1), nrow = n)
  F0 <- NULL
  for (k in 1:NP) {
    b <- BM[1:n, k]
    A_M <- Reduce("+", lapply(1:n, function(n0) (b[n0]^2) * AA[,,n0]))
    rownames(A_M) <- samplenames
    colnames(A_M) <- samplenames
    A <- similarity_kinship(A_M)
    A_train <- A[-c(1:num_ex), -c(1:num_ex)]
    res_M1 <- Pre_ER5(trainop1_ex, A_train, num_test)
    F0[k] <- res_M1[length(res_M1)]
  }
  Training_Accuracy[1:NP, ncol(Training_Accuracy)] <- F0
  Training_Weight[1:n, ncol(Training_Weight)] <- BM[, which.max(F0)]
  F_N0 <- matrix(0, nrow = NP, ncol = iter)
  F_N0[, 1] <- F0
  
  for (j in 2:iter) {
    V <- matrix(0, nrow = n, ncol = NP)
    U <- matrix(0, nrow = n, ncol = NP)
    U <- BM
    for (l in 1:NP) {
      Index_NP <- c(1:NP)
      r <- sample(Index_NP[-l])[1:3]
      V[, l] <- BM[, r[1]] + Mu * (BM[, r[2]] - BM[, r[3]])
    }
    randrow <- runif(n)
    if (min(randrow) > CR) {
      rid <- sample(nrow(BM))[1]
      U[rid, ] <- V[rid, ]
    } else {
      U[which(randrow <= CR), ] <- V[which(randrow <= CR), ]
    }
    
    for (k in 1:NP) {
      U[U[, k] <= s0, k] <- s0
      U[U[, k] > s1, k] <- s1
      if (max(U[, k]) == 0) {
        U[, k] <- BM[, k]
      }
      b <- U[1:n, k]
      A_M <- Reduce("+", lapply(1:n, function(n0) (b[n0]^2) * AA[,,n0]))
      rownames(A_M) <- samplenames
      colnames(A_M) <- samplenames
      A <- similarity_kinship(A_M)
      A_train <- A[-c(1:num_ex), -c(1:num_ex)]
      res_M1 <- Pre_ER5(trainop1_ex, A_train, num_test)
      F_N0[k, j] <- res_M1[length(res_M1)]
      if (F_N0[k, j] < F0[k]) {
        U[, k] <- BM[, k]
        F_N0[k, j] <- F0[k]
      }
    }
    BM <- U
    if (max(F_N0[, j] - F_N0[, (j - 1)]) < thre) {
      F0 <- F_N0[, j]
      Training_Accuracy[((j - 1) * NP + 1):(j * NP), ncol(Training_Accuracy)] <- F0
      Training_Weight[((j - 1) * n + 1):(j * n), ncol(Training_Weight)] <- U[, which.max(F0)]
      break
    }
    F0 <- F_N0[, j]
    Training_Accuracy[((j - 1) * NP + 1):(j * NP), ncol(Training_Accuracy)] <- F0
    Training_Weight[((j - 1) * n + 1):(j * n), ncol(Training_Weight)] <- U[, which.max(F0)]
  }
  MVBLUP_information <- data.frame(parameters = c(paste0("View", 1:n, "_W"), "iter", "Trait"),
                                   values = c(BM[, which.max(F0)], j, colnames(dt_ex)[trid]))
  
  results <- list(Training_Accuracy = Training_Accuracy, 
                  Training_Weight = Training_Weight, 
                  MVBLUP_information = MVBLUP_information, 
                  kinship_test = A,
                  trainop1_ex = trainop1_ex,
                  testop1_ex = testop1_ex)
  return(results)
}

#=====Using the training set, we derive the ideal weight configuration and apply it to the test set to forecast phenotypic precision=======
MVBLUP_PreTest <- function(A_test, testop1_ex, trainop1_ex) {
  Tres_M1 <- T_ER(trainop1_ex, testop1_ex, A_test)
  Coe<-data.frame(name = c(names(Tres_M1)[-length(Tres_M1)],"Accuracy"), value = Tres_M1)
  rownames(Coe) <- 1:nrow(Coe)
  return(Coe)
}

