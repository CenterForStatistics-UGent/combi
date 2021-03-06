#' Evaluate the influence function
#'
#' @param modelObj The model object
#' @param samples A boolean, should we look at sample variables? Throws an error otherwise
#' @param Dim,View Integers, the dimension and views required
#'
#' @method influence combi
#' @details Especially the influence of the different views on the latent
#' variable or gradient estimation may be of interest. The influence values are
#' not all calculated. Rather, the score values and inverse jacobian are returned
#' so they can easily be calculated
#' @importFrom Matrix diag solve
#' @return A list with components
#' \item{score}{The evaluation of the score function}
#' \item{InvJac}{The inverted jacobian matrix}
influence.combi = function(modelObj, samples = is.null(View), Dim = 1, View = NULL){
    with(modelObj, {
    if(samples){
        lambdaLatent = lambdasLatent[seqM(Dim, normal = FALSE)]
        constrained = !is.null(covariates)
        if(constrained) covMat = buildCovMat(covariates)$covModelMat
        #Find the normalization lagrange mutliplier?
        score = Reduce(f = cbind, lapply(seq_along(data), function(i){
            if(distributions[[i]] == "gaussian"){
                mu = buildMu(offSet = buildOffsetModel(modelObj, i),
                             latentVar = latentVars[,Dim],
                    paramEsts = paramEsts[[i]][Dim,], distributions[[i]])
                rowMultiply(if(constrained) crossprod(covMat, data[[i]] - mu) else
                    data[[i]] - mu, paramEsts[[i]][Dim,]/varPosts[[i]])
            } else if(distributions[[i]] == "quasi"){
                if(compositional[[i]]){
                    CompMat = buildCompMat(indepModels[[i]]$colMat,
                                           paramEsts[[i]],
                                           latentVar = latentVars, m = Dim,
                                           norm = TRUE)
                    mu = CompMat * indepModels[[i]]$libSizes
                    tmpMat = CompMat*(data[[i]]-mu)/meanVarTrends[[Dim]][[i]](CompMat, outerProd = FALSE)*
                        (matrix(paramEsts[[i]][Dim,], byrow = TRUE, nrow(mu), ncol(mu)) -
                             c(CompMat %*% paramEsts[[i]][Dim,]))
                    if(constrained){
                        tmpMat = crossprod(covMat, tmpMat)
                    }
                    return(tmpMat)
                } else {
                    prepMat = prepareScoreMat(mu = mu, data = data[[i]], meanVarTrend = meanVarTrends[[Dim]][[i]])
                    rowMultiply(if(constrained) crossprod(covMat, prepMat) else prepMat, paramEsts[[i]][Dim,])
                }
           }})) + (lambdaLatent[1] +
                       if(Dim==1) 0 else c(latentVars[, seq_len(Dim-1), drop = FALSE] %*% lambdaLatent[-1]))
    # Inverse Jacobian
        n = if(constrained) ncol(covMat) else nrow(data[[1]])
        nLambda1s = if(constrained) nrow(centMat) else 1
        Jacobian = buildEmptyJac(n = n, m = Dim,
                                 lower = if(constrained) alphas else latentVars,
                                 nLambda1s = nLambda1s,
                                 normal = constrained, centMat = centMat)
        if(constrained){
            Jacobian[seq_len(n),seq_len(n)] = rowSums(dims = 2, vapply(seq_along(data), FUN.VALUE = diag(double(n)), function(i){
                jacLatentVarsConstr(data = data[[i]], distribution = distributions[[i]],
                                    paramEsts = paramEsts[[i]], offSet = buildOffsetModel(modelObj, i),
                                    latentVar = latentVars[, Dim], meanVarTrend = meanVarTrends[[Dim]][[i]],
                                    varPosts = varPosts[[i]], mm = Dim, covMat = covMat,
                                    latentVarsLower = latentVars[, seq_len(Dim-1)], compositional = compositional[[i]],
                                    indepModel = indepModels[[i]], n = n,
                                    allowMissingness = modelObj$allowMissingness,
                                    paramMats = matrix(paramEsts[[i]][Dim,], byrow = TRUE, nrow(data[[i]]), ncol(data[[i]])))
            }))
            numCov = ncol(covMat)
            Jacobian[seq_len(numCov), numCov+1+nLambda1s] = Jacobian[numCov+1+nLambda1s, seq_len(numCov)] =
                2*alphas[seq_len(numCov), Dim]
        } else {
            diag(Jacobian)[seq_len(n)] = rowSums(vapply(seq_along(data), FUN.VALUE = double(n), function(i){
                jacLatentVars(data = data[[i]], distribution = distributions[[i]],
                              paramEsts = paramEsts[[i]], offSet = buildOffsetModel(modelObj, i),
                              latentVar = latentVars[, Dim], meanVarTrend = meanVarTrends[[Dim]][[i]],
                              varPosts = varPosts[[i]], mm = Dim,
                              latentVarsLower = latentVars[, seq_len(Dim-1)], compositional = compositional[[i]],
                              indepModel = indepModels[[i]], n = n,
                              allowMissingness = modelObj$allowMissingness,
                              paramMats = matrix(paramEsts[[i]][Dim,], byrow = TRUE, nrow(data[[i]]), ncol(data[[i]])))
            }))
        }

    InvJac = solve(Jacobian)[seq_len(n), seq_len(n)]
    } else {
        stop("Only influence functions of sample variables are implemented currently!\n")
    }
    # Matrix of all influences becomes too large: return score and
    # inverse jacobian and calculate influences on demand
    return(list(score = score, InvJac = InvJac))
    })
}