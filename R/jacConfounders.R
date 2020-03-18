#' Jacobian when estimating confounder variables
#' @inheritParams scoreConfounders
#' @param distribution,offSet distribution and offset of the view
#' @param libSizes,CompMat Library sizes and relative abunance
#' @return the jacobian matrix
jacConfounders = function(confMat, data, distribution, x, meanVarTrend,
                       offSet, CompMat, libSizes, allowMissingness){
    if(distribution == "gaussian"){
        -crossprod(confMat)
    } else if(distribution == "quasi"){
        mu = offSet*exp(confMat %*% x)
        if(allowMissingness){
            isNA = is.na(data)
            data[isNA] = mu[isNA]
        }
        crossprod(confMat * c(prepareJacMat(data = data, mu = mu,
                                    meanVarTrend = meanVarTrend,
                                    CompMat = CompMat, libSizes = libSizes)),
                  confMat)
    }
}