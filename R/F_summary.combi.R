#' Show an overview of a fitted combi object
#' @param object a fitted combi object
#'
#' @return A ggplot object containing the plot
#' @method summary combi
#'
#' @export
#' @examples
#' data(Zhang)
#' \donttest{
#' #Unconstrained
#' microMetaboInt = combi(
#' list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
#' distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
#' logTransformGaussian = FALSE, verbose = TRUE)
#' #Constrained
#' microMetaboIntConstr = combi(
#'     list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
#'     distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
#'     logTransformGaussian = FALSE, covariates = zhangMetavars, verbose = TRUE)}
#'     #Load the fits
#' load(system.file("extdata", "zhangFits.RData", package = "combi"))
#' summary(microMetaboInt)
#' summary(microMetaboIntConstr)
summary.combi = function(object, ...){
    constr = if(is.null(object$covariates)) "Unconstrained" else "Constrained"
    dim = length(object$iter)
    datSets = length(object$data)
    nSam = nrow(object$latentVars)
    viewsString = vapply(seq_along(object$data), FUN.VALUE = character(1),
                         function(i){
        paste0(names(object$data[i]),": ", ncol(object$paramEsts[[i]]), "\n")
    })
    cat(constr, "combi ordination of", dim, "dimensions on",
        datSets, "views with", nSam,
        "samples.\nViews and number of features were:\n", viewsString)
}