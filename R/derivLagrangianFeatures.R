#' The score function to estimate the feature parameters
#' @inheritParams estFeatureParameters
#' @param ... arguments to the jacobian function, currently ignored
#' @param distribution,compositional,meanVarTrend,offSet,numVar,indepModel,paramEstsLower
#' Characteristics of the view
#' @param mm the current dimension
#' @param x current parameter estimates
#'
#' @return A vector with the evaluation of the score functions of the feature parameters
derivLagrangianFeatures = function(x, data, distribution, offSet, latentVars,
                                numVar, paramEstsLower, mm, indepModel,
                                meanVarTrend, weights, compositional, ...){
    #Extract the parameter estimates
    param = x[seq_len(numVar)]
    score = scoreFeatureParams(data = data, distribution = distribution,
                            x = param, offSet = offSet,
                            latentVar = latentVars, meanVarTrend = meanVarTrend,
                            mm = mm, compositional = compositional,
                            indepModel = indepModel, paramEstsLower = paramEstsLower, ...) +
        weights*(x[numVar+1] + (2*param*x[numVar+2]) + (if(mm==1) 0 else (x[(numVar+3):(numVar+mm+1)] %*% paramEstsLower)))
    #Extract lagrange multipliers immediately
    centering = sum(param*weights)
    norma = sum(param^2*weights) - 1
    orthogonality = if(mm==1) NULL else (paramEstsLower %*% (param*weights))
        return(c(score, centering, norma,  orthogonality))
}