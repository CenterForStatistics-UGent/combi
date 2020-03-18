#' Prepare an empty Jacobian matrix, with useful entries prefilled.
#' In case of distribution "gaussian", it returns the lhs matrix of the linear system for finding the feature paramters
#' @param n the number of parameters
#' @param m the dimension
#' @param lower the current parameter estimates
#' @param distribution A character string, the distributional assumption for the data
#' @param normal a boolean, are normalization restrictions in place?
#' @param nLambda1s The number of centering restrictions
#' @param centMat The centering matrix
#' @param weights Vector of feature weights
#'
#' @importFrom Matrix sparseMatrix
#'
#' @return an empty jacobian matrix, or the lhs of the system of estimating equations
buildEmptyJac = function(n, m, lower, distribution = "quasi", normal = FALSE,
                         nLambda1s = 1, centMat = NULL, weights = 1){
    Dim = n + m + normal + nLambda1s - 1
    if(distribution=="quasi"){
    i = rep(seq_len(n), m + nLambda1s - 1)
    j = c(rep(c(n + seq_len(nLambda1s),
        if(m>1) seq(n+1+normal+nLambda1s, n+m+normal+nLambda1s-1) else NULL), each = n))
    sparseMatrix(i = i, j = j, dims = c(Dim, Dim), symmetric = TRUE,
                 x = c(if(nLambda1s==1) {if(length(weights)==1) rep(weights, n) else
                     weights} else t(centMat), lower[, seq_len(m-1)]*weights))
    } else if(distribution == "gaussian"){
        #Not a real jacobian, but lhs of the system of linear equations
    i1 = rep(seq_len(n), m + nLambda1s - 1)
    j1 = c(rep(c(n + seq_len(nLambda1s),
                    if(m>1) seq(n+1+nLambda1s, n+m+nLambda1s-1) else NULL), each = n))
    x1 = -weights*cbind(1, lower[, seq_len(m-1)])
    i2 = rep(seq(n+1, n + m + nLambda1s - 1), each =  n)
    x2 = weights *cbind(1, lower[, seq_len(m-1)])
        sparseMatrix(i = c(i1,i2), j = c(j1,i1), dims = c(Dim, Dim), symmetric = FALSE,
                     x = c(x1,x2))
    }
}