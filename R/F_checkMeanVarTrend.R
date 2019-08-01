#' Quickly check if the mean variance trend provides a good fit
#'
#' @param data The n-by-p data matrix
#' @param meanVarFit The type of mean variance fit
#'
#' @return A plot object
checkMeanVarTrend = function(data, meanVarFit = "spline", returnTrend = FALSE, ...){
    libSizes = rowSums(data)
    baseAbundances = colSums(data)/sum(data)
    data = data[libSizes>0, baseAbundances>0]
    libSizes = libSizes[libSizes>0]
    baseAbundances = baseAbundances[baseAbundances>0]
    foo = estMeanVarTrend(data, meanMat = outer(libSizes, baseAbundances),
                    plot = TRUE, meanVarFit = meanVarFit, baseAbundances = baseAbundances,
                    libSizes = libSizes, ...)
    if(returnTrend) return(foo) else return(invisible())
}