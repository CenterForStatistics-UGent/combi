#' print an overview of a fitted combi x
#' @param x a fitted combi x
#' @param ... Further arguments, currently ignored
#'
#' @return A ggplot x containing the plot
#' @method print combi
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
#' print(microMetaboInt)
#' print(microMetaboIntConstr)
print.combi = function(x, ...){
    constr = if(is.null(x$covariates)) "Unconstrained" else "Constrained"
    dim = length(x$iter)
    datSets = length(x$data)
    nSam = nrow(x$latentVars)
    viewsString = vapply(seq_along(x$data), FUN.VALUE = character(1),
                         function(i){
        paste0(names(x$data[i]),": ", ncol(x$paramEsts[[i]]), "\n")
    })
    cat(constr, "combi ordination of", dim, "dimensions on",
        datSets, "views with", nSam,
        "samples.\nViews and number of features were:\n", viewsString,
        if(!is.null(x$covariates)) {
            paste0("Number of sample variables included was ",ncol(x$covariates),
                  ",\nfor which ", ncol(x$covMat),
                  " parameters were estimated per dimension.")
            })
}