#' Check for monotonicity in compositional datasets fro given dimensions
#' @param modelObj The combi fit
#' @param Dim The dimensions considered
#'
#'@return A boolean matrix indicating monotonicity for every feature
checkMonotonicity = function(modelObj, Dim){
    #Check (local) monotonicity for compositionality
    with(modelObj, {
    lapply(seq_along(data), function(i){
        if(!compositional[[i]]){return(NULL)}
            sign(latentVars[,Dim, drop = FALSE] %*% paramEsts[[i]][Dim,,drop = FALSE]) ==
                      sign(buildCompMat(colMat = indepModels[[i]]$colMat, paramEsts = paramEsts[[i]],
                                         latentVar = latentVars, m = NA, id = Dim) -
                                exp(indepModels[[i]]$colMat)/rowSums(exp(indepModels[[i]]$colMat)))
        })
    })
}