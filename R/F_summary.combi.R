#' Show an overview of a fitted combi object
#' @param comb a fitted combi object
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
summary.combi = function(comb){
    constr = if(is.null(comb$covariates)) "Unconstrained" else "Constrained"
    dim = length(comb$iter)
    datSets = length(comb$data)
    nSam = nrow(comb$latentVars)
    viewsString = sapply(seq_along(comb$data), function(i){
        paste0(names(comb$data[i]),": ", ncol(comb$paramEsts[[i]]), "\n")
    })
    cat(constr, "combi ordination of", dim, "dimensions on",
        datSets, "views with", nSam,
        "samples.\nViews and number of features were:\n", viewsString)
}