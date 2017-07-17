# Feature Specific Quantile Normalization (FSQN)
# Author: Jennifer Franks, jennifer.m.franks.gr@dartmouth.edu
#
#---------quantileNormalizeByFeature----------
#INPUT
#   to_normalize = matrix (m x n) with m samples as rows, and n features as
#                   columns
#   target_dist = matrix (m2 x n) with m2 samples as rows, and n features as
#                   columns to use as the target distribution
#
# OUTPUT
#   data.qn = normalized matrix (m x n)
#
#
#---------coordinateMatrices----------
#
#INPUT
#   matrix1 = matrix (m1 x n1) with rownames and oclnames
#   matrix2 = matrix (m2 x n2) with rownames and colnames
#
#OUTPUT
#   list (matrix1, matrix2), with ordered corresponding columns
#
# ------------------------------------------------

library(preprocessCore)

quantileNormalizeByFeature <- function(matrix_to_normalize,
                                       target_distribution_matrix){

    if (ncol(matrix_to_normalize) != ncol(target_distribution_matrix)){
        cat("ERROR: Data matrices are not compatible - column lengths differ!")
    }
    else{

        data.qn <- matrix(0, nrow = nrow(matrix_to_normalize),
                          ncol = ncol(matrix_to_normalize))

        for (i in 1:ncol(matrix_to_normalize)){
            feature.to.normalize <- matrix_to_normalize[,i]
            target.feature.dist <- target_distribution_matrix[,i]
            result <- normalize.quantiles.use.target(
                x = as.matrix(feature.to.normalize),
                target = target.feature.dist,
                copy = TRUE)
            data.qn[,i] <- result
        }
        rownames(data.qn) = rownames(matrix_to_normalize)
        colnames(data.qn) = colnames(matrix_to_normalize)
        return(data.qn)
    }
}



coordinateMatrices <- function(matrix1, matrix2){
    if (length(colnames(matrix1)) == 0 | length(colnames(matrix2))==0){
        cat("ERROR: At least one matrix is missing colnames!")
    }
    else{
        comm_genes <- list()
        comm_genes$matrix1 <- colnames(matrix1)
        comm_genes$matrix2 <- colnames(matrix2)
        common.symbols <- Reduce(intersect, comm_genes)

        # need to get the matrices to match by colnames using what is in common
        matrix1 <- matrix1[,which(colnames(matrix1) %in% common.symbols)]
        matrix1 <- matrix1[,order(colnames(matrix1))]

        matrix2 <- matrix2[,which(colnames(matrix2) %in% common.symbols)]
        matrix2 <- matrix2[,order(colnames(matrix2))]

        return(list(V1=matrix1,V2=matrix2))
    }

}
