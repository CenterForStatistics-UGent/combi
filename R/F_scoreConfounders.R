#' Score functions for confounder variables
#'
#' @param x
#' @param data
#' @param distribution
#' @param offSet
#' @param confMat
#' @param meanVarTrend
#' @param allowMissingness
#'
#' @return The evaluation of the estimating equations
scoreConfounders = function(x, data, distribution, offSet, confMat,
                              meanVarTrend, allowMissingness){
    if(distribution == "gaussian"){
        mu = offSet + confMat %*% x
        if(allowMissingness){
            isNA = is.na(data)
            data[isNA] = mu[isNA]
        }
        crossprod(confMat, (data - mu))
    } else if(distribution == "quasi"){
        mu = offSet*exp(confMat %*% x)
        crossprod(confMat, prepareScoreMat(mu = mu, data = data,
                                      meanVarTrend = meanVarTrend))
    } else if(distribution == "binomial"){

    }
}