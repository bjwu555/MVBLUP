# The script functioned to derive the kinship matrix by defining a similarity function
# Notes:
# data: the initial single-view dataset multiplied by its transpose
# Multis: a list composed of multiple viewpoint datasets
# n: the number of multi-view datasets
#============Calculate the kinship matrix using the similarity function=========
similarity_kinship <- function(A_M) {
  A_N0 <- A_M
  for (i in 1:nrow(A_N0)) {
    for (j in 1:nrow(A_N0)) {
      A_N0[i, j] <- A_M[i, j] / sqrt(A_M[i, i] * A_M[j, j])
    }
  }
  return(A_N0)
}
#=============Computing the inner product of a single-view dataset==============
similarity_inner_product <- function(Multis,n) {
  A_N <- array(0, dim = c(nrow(Phenotype), nrow(Phenotype), n))
  for (n0 in 1:n) {
    data_S <- apply(Multis[[n0]][, -1], 2, function(x) {
      (as.numeric(x) - mean(as.numeric(x))) / sd(as.numeric(x))
    })
    A_M <- data_S %*% t(data_S)
    A_M[diag(A_M) == 0] <- 1
    A_N[, , n0] <- A_M
  }
  return(A_N)
}