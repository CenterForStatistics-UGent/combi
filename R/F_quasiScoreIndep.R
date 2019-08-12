#' Quasi score equations for column offset parameters of sequence count data
#'
#' @param x the initial guess for the current margin
#' @param data the data matrix
#' @param meanVarTrend the function describing the mean-variance trend
#' @param otherMargin The other margin
#' @param col A logical, is the column being estimated?
#' @param libSizes The library sizes
#' @param ... passed on to prepareJacMat
#'
#' @return the evaluated estimating equation
quasiScoreIndep = function(x, data, otherMargin, meanVarTrend, col, libSizes, ...){
    mu = exp(buildMuMargins(x = x, otherMargin = otherMargin, col = col))
    if(anyNA(mu) | any(is.infinite(mu)) ){
        return(rep(1e16, length(x)))
    }
    (if(col) colSums else rowSums)(prepareScoreMat(data = data, mu = mu,
                                      meanVarTrend = meanVarTrend,
                                      CompMat = exp(if(col) x else otherMargin),
                                      libSizes = libSizes),
                                   na.rm = TRUE)
}