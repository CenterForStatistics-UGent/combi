#' prepare the jacobian for the latent variabels compostional
#'
#' @param mu the mean matrix
#' @param data the count matrix
#' @param meanVarTrend The mean variance trend
#' @param CompMat0 The compisition matrix
#' @param libSizes The library sizesv
#' @param paramEsts Current parameter estimates
#'
#' @return The empty jacobian matrix with entries maximally filled out
prepareJacMatComp = function(mu, paramEsts, CompMat0, meanVarTrend, data, libSizes){
    Sum = rowSums(CompMat0)
    CompMat = CompMat0/Sum
    CompMat[CompMat==0] = .Machine$double.eps
    meanEval = meanVarTrend(CompMat, libSizes, outerProd = FALSE)*libSizes
    CC = c(CompMat0 %*% paramEsts)
    foo = mu*CC/Sum
    tmp = -((rowMultiply(mu, paramEsts) - foo)^2*
                ((data - mu)*meanVarTrend(CompMat, deriv = 1L)/meanEval + 1) -
                (rowMultiply(mu, paramEsts^2) - 2*rowMultiply(foo, paramEsts) +
                     2*mu*(CC/Sum)^2 - mu*c(CompMat0 %*% paramEsts^2)/Sum)*
                (data - mu))/meanEval
    tmp
}