#' Score equations for conditioning under compositionality
#'
#' @param x Confounder parameter estimates
#' @param confMat confounder matrix
#' @param data data
#' @param meanVarTrend mean variance trend
#' @param marginModel marginal models
#' @param biasReduction A boolean, should a bias reduced estimation be applied?
#' @param allowMissingness A boolean, are missing values allowed
#' @param subtractMax A boolean, shuold the maximum be subtracted before softmax
#' transformation? Recommended for numerical stability
#'
#' @return The evaluation of the estimating equations
scoreConfoundersComp = function(x, confMat, data, meanVarTrend, marginModel,
                                biasReduction, allowMissingness,
                                subtractMax = TRUE){
    parMat = matrix(x, ncol(confMat), ncol(data)-1)
    parMat = cbind(0, parMat)# The additive log-ratio transform
    logCompMat = matrix(marginModel$colOff, nrow(data), ncol(data), byrow = TRUE) +
        + confMat %*% parMat
    if(subtractMax) logCompMat = logCompMat - apply(logCompMat, 1, max)
    CompMat0 = exp(logCompMat)
    Sum = rowSums(CompMat0)
    CompMat = CompMat0/Sum
    mu = CompMat*exp(marginModel$rowOff)
    if(allowMissingness){
        isNA = is.na(data)
        data[isNA] = mu[isNA]
    }
    tmp = crossprod(confMat/Sum,
              ((CompMat*(Sum-CompMat))*(data-mu)/meanVarTrend(CompMat, outerProd = FALSE))[,-1])#
    if(biasReduction) tmp = tmp -
        c((diag(confMat %*% tcrossprod(solve(crossprod(confMat)), confMat))/2) %*%
             confMat)
    return(c(tmp))
}