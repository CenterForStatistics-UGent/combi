#' Evaluate the score functions for the estimation of the feature parameters for a
#' single dataset
#' @inheritParams estFeatureParameters
#' @return A vector with evaluated score function
#'@param data,distribution,offSet,meanVarTrend,indepModel,compositional,paramEstsLower
#' Characteristics of the views
#' @param x the parameter estimates
#' @param latentVar the latent variables
#' @param mm the dimension
#' @param allowMissingness a boolean, should missing values be allowed
scoreFeatureParams = function(x, data, distribution, offSet, latentVar,
                           meanVarTrend, mm, indepModel, compositional,
                           paramEstsLower, allowMissingness, ...){
    if(!compositional){
        mu = buildMu(offSet, latentVar[,mm], x,
                     distribution)
        if(allowMissingness){
            isNA = is.na(data)
            data[isNA] = mu[isNA]
        }
    }
    if(distribution == "gaussian"){
        latentVar[,mm] %*% (data - mu)
    } else if(distribution == "quasi"){
        if(compositional){
        CompMat = buildCompMat(indepModel$colMat, rbind(paramEstsLower, x),
                               latentVar, m = mm, norm = TRUE, subtractMax = TRUE)
        mu = CompMat*indepModel$libSizes
        if(allowMissingness){
            isNA = is.na(data)
            data[isNA] = mu[isNA]
        }
        CompMatVar = CompMat/meanVarTrend(CompMat, outerProd = FALSE)
        crossprod(latentVar[,mm],
                  (1-CompMat)*(data-mu)*CompMatVar)
        } else {
        latentVar[,mm] %*% prepareScoreMat(mu = mu, data = data,
                        meanVarTrend = meanVarTrend)
       }
    } else if(distribution == "binomial"){

    }
}