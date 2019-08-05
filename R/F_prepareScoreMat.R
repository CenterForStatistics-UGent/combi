#' Prepare a helper matrix for score function evaluation under quasi-likelihood
#'
#' @param data
#' @param mu
#' @param meanVarTrend
#' @param CompMat
#' @param libSizes
#'
#' @return The helper matrix
prepareScoreMat = function(data, mu, meanVarTrend, CompMat, libSizes){
    (data - mu)*mu/meanVarTrend(c(CompMat), libSizes = libSizes)
}