#' Score equations for conditioning under compositionality
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
    # if(anyNA(mu)){
    #     return(rep(1e16, length(x)))
    # }
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