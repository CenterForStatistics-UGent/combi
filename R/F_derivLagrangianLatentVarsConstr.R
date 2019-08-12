#' The score function to estimate the latent variables
#' @param ... arguments to the jacobian function, currently ignored
#' @param x The current estimates of the latent variables
#' @param n The number of samples
#' @param m The dimensions
#' @param numSets The number of views
#' @param latentVarsLower The parameter estimates of the lower dimensions
#' @param compositional,links,indepModels,meanVarTrends,numVars,distributions,data,offsets,varPosts,paramMats,paramEsts
#' Lists of information on all the views
#' @param covMat The covariance matrix
#' @param centMat A centering matrix
#' @param nLambda1s The number of dummy variables
#' @param numCov The number of covariates
#'
#' @return A vector of length n, the evaluation of the score functions of the latent variables
derivLagrangianLatentVarsConstr = function(x, data, distributions, offsets, paramEsts,
                                     numVars, latentVarsLower, n, m,
                                     numSets, meanVarTrends, links, covMat,
                                     numCov, centMat, nLambda1s, varPosts,
                                     compositional, indepModels, paramMats,
                                     ...){
    #Extract the latent variable estimates
    alpha = x[seq_len(numCov)]
    # latent = c(covMat %*% alpha)
    # latentVarsLower = covMat %*% latentVarsLower
score = rowSums(
vapply(seq_len(numSets), FUN.VALUE = alpha, function(i){
 c(scoreLatentVars(data = data[[i]], distribution = distributions[[i]],
                 paramEsts = paramEsts[[i]], offSet = offsets[[i]],
                 paramMats = paramMats[[i]],
                 latentVar = alpha, meanVarTrend = meanVarTrends[[i]],
                 nLambda1s = nLambda1s, constrained = TRUE, covMat = covMat,
                 varPosts = varPosts[[i]], compositional = compositional[[i]],
                 indepModel = indepModels[[i]], latentVarsLower = latentVarsLower,
                 mm = m, ...))
})) + c(x[numCov + seq_len(nLambda1s)] %*% centMat) +
    2*x[numCov + nLambda1s + 1]*alpha +
    if(m==1) 0 else
        c(latentVarsLower %*% x[(numCov+2+nLambda1s):(numCov+nLambda1s+m)])
#Extract lagrange multipliers immediately
centerFactors = centMat %*% alpha
if(m==1){
    return(c(score, centerFactors, sum(alpha^2)-1))
} else {
    return(c(score, centerFactors, sum(alpha^2)-1,
             crossprod(latentVarsLower, alpha)))
}
}