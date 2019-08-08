#' prepare the jacobian matrix
#'
#' @param mu the mean matrix
#' @param data the count matrix
#' @param meanVarTrend The mean variance trend
#' @param CompMat The compisition matrix
#' @param libSizes The library sizes
#'
#' @return the matrix which can be summed over
prepareJacMat = function(mu, data, meanVarTrend, CompMat, libSizes){
    trendEval = meanVarTrend(c(CompMat), libSizes = libSizes)
    derivEval = meanVarTrend(c(CompMat), libSizes = libSizes, deriv = 1L)
    mu*((data-2*mu)*trendEval -
        rowMultiply(mu*(data-mu), derivEval))/
    trendEval^2
}