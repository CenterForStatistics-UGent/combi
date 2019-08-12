#' Jacobian when estimating confounder variables
#' @inheritParams scoreConfounders
#' @param distribution,offSet distribution and offset of the view
#' @return the jacobian matrix
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