#' Evaluate the score functions for the estimation of the latent variables for a
#' single dataset
#'
#' @inheritParams estLatentVars
#'
#' @return A vector of length n, with evaluated score function
#' @param data,distribution,offSet,meanVarTrend,indepModel,varPosts,paramEsts,paramMats,compositional
#' Characteristics of the views
#' @param latentVar the latent variable estimates
#' @param covMat a matrix of constraining covariates
#' @param constrained a boolean, is this a constrained analysis
#' @param mm the current dimension
#' @param latentVarsLower the lower dimensional latent variables
#' @param allowMissingness a boolean, should missing values be allowed
scoreLatentVars = function(data, distribution, paramEsts, paramMats, offSet, latentVar,
                           meanVarTrend, constrained = FALSE, covMat = NULL,
                           varPosts, compositional, indepModel, mm,
                           latentVarsLower, allowMissingness, ...){
    if(!compositional){
    mu = buildMu(offSet = offSet, latentVar = if(constrained)
        c(covMat %*% latentVar) else latentVar, paramEsts = paramMats,
        distribution, paramMatrix = TRUE)
    # if(anyNA(mu)){
    #     return(rep(1e16, length(latentVar)))
    # }
    if(allowMissingness){
        isNA = is.na(data)
        data[isNA] = mu[isNA]
    }
    }
# | any(is.infinite(mu))
if(distribution == "gaussian"){
    (if(constrained) crossprod(covMat, data - mu) else
      data - mu) %*% (paramEsts[mm,]/varPosts)
} else if(distribution == "quasi"){
    if(compositional){
    CompMat = buildCompMat(indepModel$colMat, paramEsts,
                           latentVar = if(constrained){
                               covMat %*% cbind(latentVarsLower, latentVar)} else{
                                   cbind(latentVarsLower, latentVar)}, m = mm,
                               norm = TRUE)
    mu = CompMat*indepModel$libSizes
    # if(anyNA(mu)){
    #     return(rep(1e16, length(latentVar)))
    # }
    if(allowMissingness){
        isNA = is.na(data)
        data[isNA] = mu[isNA]
    }
    CompMatVar = CompMat/meanVarTrend(CompMat, outerProd = FALSE)
    # CompMatVar[is.na(CompMatVar)] = 1
    #This can be improved
    tmpMat  = CompMatVar*(data-mu)*
        (paramMats - c(CompMat %*% paramEsts[mm,]))
    if(constrained){
    tmpMat = crossprod(covMat, tmpMat)
    }
    rowSums(tmpMat)
    } else {
    prepMat = prepareScoreMat(mu = mu, data = data, meanVarTrend = meanVarTrend)
    (if(constrained) crossprod(covMat, prepMat) else prepMat) %*% paramEsts[mm,]
    }
} else if(distribution == "binomial"){

}
}
