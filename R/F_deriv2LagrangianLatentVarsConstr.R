#' The score function to estimate the latent variables
#' @inheritParams deriv2LagrangianLatentVars
#' @param Jac an empty jacobian matrix
#' @param ... arguments to the jacobian function, currently ignored
#'
#' @return A vector of length nn, the evaluation of the score functions of the latent variables
deriv2LagrangianLatentVarsConstr = function(x, data, distributions, offsets, paramEsts, paramMats,
                                     numVars, latentVarsLower, nn, m, Jac,
                                     numSets, meanVarTrends, links, numCov,
                                     covMat, nLambda1s, varPosts,
                                     compositional, indepModels,...){
    latentVar = c(covMat %*% x[seq_len(numCov)])
    latentVarsLower = covMat %*% latentVarsLower
    sepJacs = vapply(seq_len(numSets),
                     FUN.VALUE = matrix(0, nrow = numCov, ncol = numCov),
                     function(i){
            jacLatentVarsConstr(data = data[[i]], distribution = distributions[[i]],
                        paramEsts = paramEsts[[i]], offSet = offsets[[i]], paramMats = paramMats[[i]],
                        latentVar = latentVar, latentVarsLower = latentVarsLower,
                        meanVarTrend = meanVarTrends[[i]], n = numCov,
                        covMat = covMat, varPosts = varPosts[[i]], compositional = compositional[[i]],
                        indepModel = indepModels[[i]], mm = m)
        })
    Jac[seq_len(numCov), seq_len(numCov)] = rowSums(sepJacs, dims = 2)
    diag(Jac)[seq_len(numCov)] = diag(Jac)[seq_len(numCov)] + 2*x[numCov + nLambda1s + 1]
    Jac[seq_len(numCov), numCov+1+nLambda1s] = Jac[numCov+1+nLambda1s, seq_len(numCov)] =
        2*x[seq_len(numCov)]
    #Extract Lagrange multipliers immediately
    return(Jac)
}