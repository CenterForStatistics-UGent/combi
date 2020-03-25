#' Quickly check if the mean variance trend provides a good fit
#'
#' @param data Data in any acceptable format (see details ?combi)
#' @param meanVarFit The type of mean variance fit, either "cubic" or "spline"
#' @param returnTrend A boolean, should the estimated trend be returned (TRUE) or
#' only plotted (FALSE)?
#' @param ... passed on to the estMeanVarTrend() function
#'
#' @return A plot object
#' @export
#' @examples
#' data(Zhang)
#' par(mfrow = c(1,2))
#' lapply(list("microbiome" = zhangMicrobio, "metabolome" = zhangMetabo),
#'  checkMeanVarTrend)
#'  par(mfrow = c(1,1))
checkMeanVarTrend = function(data, meanVarFit = "spline", returnTrend = FALSE,
                             ...){
    data = extractData(data)
    if(!meanVarFit %in% c("cubic", "spline"))
        stop("Mean-variance trend must be either cubic or spline!\n")
    libSizes = rowSums(data)
    baseAbundances = colSums(data)/sum(data)
    data = data[libSizes>0, baseAbundances>0]
    libSizes = libSizes[libSizes>0]
    baseAbundances = baseAbundances[baseAbundances>0]
    mvt = estMeanVarTrend(data, meanMat = outer(libSizes, baseAbundances),
                    plot = TRUE, meanVarFit = meanVarFit,
                    baseAbundances = baseAbundances, libSizes = libSizes, ...)
    if(returnTrend) return(mvt) else return(invisible())
}