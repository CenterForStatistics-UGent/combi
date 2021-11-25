#' Extract coordinates from fitted object
#' @param modelObj The fitted model
#' @param Dim the required dimensions
#'
#' @return A list with components (matrices with two columns)
#' \item{latentData}{The latent variables}
#' \item{featureData}{The feature parameters}
#' \item{varData}{The variables}
#' @export
#' @examples
#' data(Zhang)
#' \dontrun{
#' #Unconstrained
#' microMetaboInt = combi(
#' list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
#' distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
#' logTransformGaussian = FALSE, verbose = TRUE)
#' }
#'     #Load the fits
#' load(system.file("extdata", "zhangFits.RData", package = "combi"))
#' extractCoords(microMetaboInt, Dim = c(1,2))
extractCoords = function(modelObj, Dim){
    dimNames = paste0("Dim", Dim)
    latentData = data.frame(modelObj$latentVars[, Dim])
    colnames(latentData) = dimNames
    featureData = lapply(modelObj$paramEsts, function(x){
        tmp  = data.frame(t(x[Dim,]))
        colnames(tmp) = dimNames
        tmp[["featNames"]] = colnames(x)
        tmp
    })
    if(!is.null(modelObj$covariates)) {
        varData = data.frame(modelObj$alphas[,Dim], "varNames" = rownames(modelObj$alphas))
    } else {varData = NULL}
    list(latentData = latentData, featureData = featureData,
         varData = varData)
}