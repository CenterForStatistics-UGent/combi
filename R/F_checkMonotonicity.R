#' Check for monotonicity in compositional datasets fro given dimensions
checkMonotonicity = function(modelObj, Dim){
    #Check (local) monotonicity for compositionality
    with(modelObj, {
    lapply(seq_along(data), function(i){
        if(!compositional[[i]]){return(NULL)}
            # colMat = matrix(indepModels[[i]]$colOff,
            #      byrow = TRUE, nrow = nrow(data[[i]]), ncol = ncol(data[[i]])) +
            # if(is.null(modelObj$confMats[[i]])) 0 else confMats[[i]]$confModelMat %*% confVars[[i]]$paramEsts
            sign(latentVars[,Dim, drop = FALSE] %*% paramEsts[[i]][Dim,,drop = FALSE]) ==
                      sign(buildCompMat(colMat = indepModels[[i]]$colMat, paramEsts = paramEsts[[i]],
                                         latentVar = latentVars, m = NA, id = Dim) -
                                exp(indepModels[[i]]$colMat)/rowSums(exp(indepModels[[i]]$colMat)))
        })
    })
}