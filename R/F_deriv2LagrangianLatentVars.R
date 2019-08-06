#' The jacobian function to estimate the latent variables
#' @inheritParams derivLagrangianLatentVars
#' @param ... arguments to the jacobian function, currently ignored
#'
#' @return A vector of length n, the evaluation of the score functions of the latent variables
deriv2LagrangianLatentVars = function(x, data, distributions, offsets, paramEsts, paramMats,
                                     numVars, latentVarsLower, n, m, Jac,
                                     numSets, meanVarTrends, links, varPosts, indepModels, compositional,...){
    sepJacs = vapply(seq_len(numSets), FUN.VALUE = numeric(n),function(i){
            jacLatentVars(data = data[[i]], distribution = distributions[[i]],
                        paramEsts = paramEsts[[i]], offSet = offsets[[i]], paramMats = paramMats[[i]],
                        latentVar = x[seq_len(n)], meanVarTrend = meanVarTrends[[i]],
                        n = n, varPosts = varPosts[[i]], indepModel = indepModels[[i]],
                        latentVarsLower = latentVarsLower, compositional = compositional[[i]],
                        mm =m, ...)})
    diag(Jac)[seq_len(n)] = rowSums(sepJacs)
    #Extract lagrange multipliers immediately
    return(Jac)
}