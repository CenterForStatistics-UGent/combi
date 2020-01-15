#' Quickly check if the mean variance trend provides a good fit
#'
#' @param data The n-by-p data matrix
#' @param meanVarFit The type of mean variance fit
#' @param returnTrend A boolean, should the estimated trend be returned or
#' only plotted?
#' @param ... passed on to the estMeanVarTrend() function
#'
#' @return A plot object
#' @export
#' @examples
#' data(Zhang)
#' par(mfrow = c(1,2))
#' lapply(extractData(list("microbiome" = zhangMicrobio, "metabolome" = zhangMetabo)),
#'  checkMeanVarTrend)
#'  par(mfrow = c(1,1))
checkMeanVarTrend = function(data, meanVarFit = "spline", returnTrend = FALSE,
                             ...){
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