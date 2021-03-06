#' Evaluate the jacobian for estimating the latent variable for one view for constrained ordination
#' @inheritParams estLatentVars
#' @importFrom tensor tensor
#' @param distribution,compositional,meanVarTrend,offSet,paramEsts,indepModel,varPosts,data
#' Characteristics of each view
#' @return The jacobian matrix
#' @param covMat the covariates matrix
#' @param mm the dimension
#' @param latentVar current latent variable estimates
#' @param latentVarsLower latent variable estimates of lower dimensions
#' @param numCov the number of covariates
jacLatentVarsConstr = function(latentVar, data, distribution, paramEsts, offSet,
                         meanVarTrend, numCov, covMat, varPosts, compositional, mm,
                         indepModel, latentVarsLower, ...){
    if(distribution == "gaussian"){
        -colSums(tensor(outer(paramEsts[mm,]^2/varPosts, covMat), covMat, 2, 1))
    } else if(distribution == "quasi"){
        if(compositional){
        CompMat0 =  buildCompMat(indepModel$colMat, paramEsts,
                        cbind(latentVarsLower, latentVar), m = mm, norm = FALSE)
        mu = CompMat0*indepModel$libSizes/rowSums(CompMat0)
        tmpMat = prepareJacMatComp(mu = mu, CompMat0 = CompMat0,
                                   paramEsts = paramEsts[mm,],
                  data = data, meanVarTrend = meanVarTrend,
                  libSizes = indepModel$libSizes)
        } else {

        mu = buildMu(offSet, cbind(latentVarsLower, latentVar), paramEsts,
                         distribution, paramMatrix = TRUE)
        tmpMat = rowMultiply(prepareJacMat(data = data, mu = mu,
                                meanVarTrend = meanVarTrend), paramEsts[mm,]^2)
        }
        colSums(tensor(vapply(seq_len(numCov),
                              FUN.VALUE = tmpMat, function(x) {
                                  covMat[, x] * tmpMat
                              }), covMat, 1, 1))
    } else if(distribution == "binomial"){

    }
}