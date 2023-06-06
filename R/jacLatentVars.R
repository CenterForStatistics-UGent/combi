#' Evaluate the jacobian for estimating the latent variable for one view
#' @inheritParams scoreLatentVars
#' @param distribution,data,varPosts,compositional,meanVarTrend,offSet,paramEsts,paramMats,indepModel Characteristics of each view
#' @return The diagonal of the jacobian matrix
#' @param n the number of samples
jacLatentVars = function(latentVar, data, distribution, paramEsts, paramMats, offSet,
                         meanVarTrend, n, varPosts, mm, indepModel,
                         latentVarsLower, compositional, allowMissingness, ...){
    if(!compositional){
        mu = buildMu(offSet = offSet, latentVar = latentVar,
                     paramEsts = paramMats, distribution, paramMatrix = TRUE)
        if(allowMissingness){
            isNA = is.na(data)
            data[isNA] = mu[isNA]
        }
    }
    if(distribution == "gaussian"){
        rep(-sum(paramEsts[mm, ]^2/varPosts), n)
    } else if(distribution == "quasi"){
        if(compositional){
            CompMat0 = buildCompMat(indepModel$colMat, paramEsts,
                                    cbind(latentVarsLower, latentVar), m = mm, norm = FALSE)
        mu = CompMat0/rowSums(CompMat0)*indepModel$libSizes
        if(allowMissingness){
            isNA = is.na(data)
            data[isNA] = mu[isNA]
        }
        rowSums(prepareJacMatComp(mu = mu, CompMat0 = CompMat0,
            paramEsts = paramEsts[mm,], data = data, meanVarTrend = meanVarTrend,
            libSizes = indepModel$libSizes)
            )
        } else {
        prepareJacMat(data = data, mu = mu, meanVarTrend = meanVarTrend) %*%
            paramEsts[mm,]^2
        }
    } else if(distribution == "binomial"){

    }
}