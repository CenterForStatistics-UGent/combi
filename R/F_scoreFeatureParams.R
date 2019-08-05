#' Evaluate the score functions for the estimation of the feature parameters for a
#' single dataset
#'
#' @inherit estFeatureParams param
#'
#' @return A vector with evaluated score function
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
    # if(anyNA(mu) ){
    #     return(rep(1e16, length(x)))
    # }# Put up a big "Keep Out" sign when drifting towards extreme means
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
        # CompMatVar[is.na(CompMatVar)] = 1
        crossprod(latentVar[,mm],
                  (1-CompMat)*(data-mu)*CompMatVar)
        } else {
        latentVar[,mm] %*% prepareScoreMat(mu = mu, data = data,
                        meanVarTrend = meanVarTrend)
       }
    } else if(distribution == "binomial"){

    }
}