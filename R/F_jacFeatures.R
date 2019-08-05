#' Evaluate the jacobian for estimating the feature parameters for one view
#'
#' @inherit estFeatureParams param
#'
#' @return The jacobian matrix
jacFeatures = function(latentVars, data, distribution, paramEsts, meanVarTrend,
                       offSet, compositional, indepModel, m, paramEstsLower,
                       minNumber = .Machine$double.eps^10, allowMissingness, ...){
    if(distribution == "gaussian"){
        -sum(latentVars[,m]^2)
    } else if(distribution == "quasi"){
        if(compositional) {
            CompMat0 = buildCompMat(indepModel$colMat, rbind(paramEstsLower, paramEsts),
                                latentVars, m = m, norm = FALSE, subtractMax = TRUE)
            Sum = rowSums(CompMat0)
            CompMat = CompMat0/Sum
            CompMat[CompMat==0] = .Machine$double.eps
            mu = CompMat*indepModel$libSizes
            if(allowMissingness){
                isNA = is.na(data)
                data[isNA] = mu[isNA]
            }
            BB = Sum - CompMat0
            meanEval = meanVarTrend(CompMat, outerProd = FALSE)*indepModel$libSizes
            # meanEvalSq = meanEval^2
            # meanEvalSq[meanEvalSq < minNumber] = minNumber
            dMudBeta =  mu*BB*latentVars[,m]/Sum
            #Prepare some matrices
            foo1 = (CompMat0-BB)*mu*
                latentVars[,m]^2/Sum^2*(data-mu)/meanEval
            foo2  = dMudBeta*(1 + (data-mu)*meanVarTrend(CompMat, deriv = 1L)/meanEval)/
                meanEval

            tmpMat = crossprod(foo1, CompMat0)
            diag(tmpMat) = -colSums(BB*foo1)

            tmpMat2 = -crossprod((mu*latentVars[,m]/Sum*foo2), CompMat0)
            diag(tmpMat2) = colSums(foo2*dMudBeta)
            tmpMat - tmpMat2

            } else {
            (latentVars[,m]^2) %*% prepareJacMat(data = data,
                mu = buildMu(offSet, latentVars[,m], paramEsts,
                             distribution), meanVarTrend = meanVarTrend)}
    } else if(distribution == "binomial"){

    }
}
