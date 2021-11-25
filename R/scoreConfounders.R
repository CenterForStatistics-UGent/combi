#' Score functions for confounder variables
#'
#' @param data,distribution,offSet,confMat,meanVarTrend
#' Characteristics of the views
#' @param x the parameter estimates
#' @param libSizes,CompMat Library sizes and relative abunance
#' @param allowMissingness a boolean, should missing values be allowed
#' @return The evaluation of the estimating equations
scoreConfounders = function(x, data, distribution, offSet, confMat,
                            meanVarTrend, allowMissingness, libSizes, CompMat){
    if(distribution == "gaussian"){
        mu = offSet + confMat %*% x
        if(allowMissingness){
            isNA = is.na(data)
            data[isNA] = mu[isNA]
        }
        crossprod(confMat, (data - mu))
    } else if(distribution == "quasi"){
        mu = offSet*exp(confMat %*% x)
        if(allowMissingness){
            isNA = is.na(data)
            data[isNA] = mu[isNA]
        }
        crossprod(confMat, (data - mu)*mu/meanVarTrend(CompMat,
                                                       libSizes = libSizes))
    }
}