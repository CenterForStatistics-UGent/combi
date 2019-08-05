#' prepare the jacobian matrix
#'
#' @param mu the mean matrix
#' @param data the count matrix
#' @param meanVarTrend
#' @param CompMat
#' @param libSizes
#'
#' @return the matrix which can be summed over
prepareJacMat = function(mu, data, meanVarTrend, CompMat, libSizes){
    trendEval = meanVarTrend(c(CompMat), libSizes = libSizes)
    derivEval = meanVarTrend(c(CompMat), libSizes = libSizes, deriv = 1L)
    mu*((data-2*mu)*trendEval -
        rowMultiply(mu*(data-mu), derivEval))/
    trendEval^2
}