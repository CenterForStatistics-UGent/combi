#' Jacobian when estimating confounder variables
jacConfounders = function(confMat, data, distribution, x, meanVarTrend,
                       offSet){
    if(distribution == "gaussian"){
        -crossprod(confMat)
    } else if(distribution == "quasi"){
        mu = offSet*exp(confMat %*% x)
        crossprod(confMat * c(prepareJacMat(data = data, mu = mu,
                                    meanVarTrend = meanVarTrend)), confMat)
    } else if(distribution == "binomial"){

    }
}