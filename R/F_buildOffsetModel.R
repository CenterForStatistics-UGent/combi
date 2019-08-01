#' Build a marginal offset matrix given a model
#' @param modelObj a modelDI object
#' @param The view for which to build the offset
#'
#' @return The offset matrix
buildOffsetModel = function(modelObj, View){
with(modelObj, {
    if(distributions[[View]]=="quasi"){
        expColMat = exp(indepModels[[View]]$colMat)
        return(indepModels[[View]]$libSizes *
            if(compositional[[View]]) expColMat/rowSums(expColMat) else expColMat)
    } else if(distributions[[View]]=="gaussian"){
        return(indepModels[[View]]$colMat + indepModels[[View]]$libSizes)
    }
})
}