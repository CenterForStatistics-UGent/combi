#' The score function to estimate the latent variables
#' @inheritParams estFeatureParameters
#' @param distribution,compositional,meanVarTrend,offSet,numVar
#' Characteristics of the view
#' @param x parameter estimates
#' @param paramEstsLower lower dimension estimates
#' @param indepModel the independence model
#'
#' @return A vector of length n, the evaluation of the score functions of the latent variables
deriv2LagrangianFeatures = function(x, data, distribution, offSet, latentVars,
                                    numVar, paramEstsLower, mm, Jac,
                                    meanVarTrend, weights, compositional,
                                    indepModel,...){
    Seq = seq_len(numVar)
    Jac[Seq, numVar + 2] = Jac[numVar + 2, Seq] = 2*x[Seq]*weights
    JacTmp = jacFeatures(data = data, distribution = distribution,
                      paramEsts = x[Seq], offSet = offSet,
                      latentVars = latentVars, meanVarTrend = meanVarTrend,
                      m = mm, compositional = compositional,
                      paramEstsLower = paramEstsLower, indepModel = indepModel, ...)
    if(compositional){
    Jac[Seq, Seq] = JacTmp
    } else {
    diag(Jac)[Seq] = JacTmp
    }
    diag(Jac)[Seq] = diag(Jac)[Seq] + 2*x[numVar+2]*weights
    return(Jac)
}