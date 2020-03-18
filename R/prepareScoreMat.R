#' Prepare a helper matrix for score function evaluation under quasi-likelihood
#'
#' @inheritParams prepareJacMat
#'
#' @return The helper matrix
prepareScoreMat = function(data, mu, meanVarTrend, CompMat, libSizes){
    (data - mu)*mu/meanVarTrend(c(CompMat), libSizes = libSizes)
}