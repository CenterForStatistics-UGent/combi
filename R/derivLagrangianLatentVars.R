#' The score function to estimate the latent variables
#' @param ... arguments to the jacobian function, currently ignored
#' @param x The current estimates of the latent variables
#' @param n The number of samples
#' @param m The dimensions
#' @param numSets The number of views
#' @param latentVarsLower The parameter estimates of the lower dimensions
#' @param compositional,links,indepModels,meanVarTrends,numVars,distributions,data,offsets,varPosts,paramMats,paramEsts
#' Lists of inforamtion on all the views
#'
#' @return A vector of length n, the evaluation of the score functions of the latent variables
derivLagrangianLatentVars = function(x, data, distributions, offsets, paramEsts,
                                     paramMats, numVars, n, m, numSets,
                                     meanVarTrends, links, varPosts,
                                     latentVarsLower, compositional, indepModels, ...){
    #Extract the latent variable estimates
    latent = x[seq_len(n)]
    score = rowSums(
    vapply(seq_len(numSets), FUN.VALUE = latent, function(i){
     scoreLatentVars(data = data[[i]], distribution = distributions[[i]],
                     paramEsts = paramEsts[[i]], offSet = offsets[[i]], paramMats = paramMats[[i]],
                     latentVar = latent, meanVarTrend = meanVarTrends[[i]],
                     varPosts = varPosts[[i]], mm = m,
                     latentVarsLower = latentVarsLower, compositional = compositional[[i]],
                     indepModel = indepModels[[i]], constrained = FALSE, ...)
    })) + x[n+1] + if(m==1) 0 else (latentVarsLower %*% x[(n+2):(n+m)])
    #Extract lagrange multipliers immediately
    if(m==1){
        return(c(score, sum(latent)))
    } else {
        return(c(score, sum(latent), crossprod(latentVarsLower, latent)))
    }
}