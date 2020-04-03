#' Helper function to extract data matrix from phyloseq, expressionset objects etc. Also filers out all zero rows
#' @param logTransformGaussian A boolean, should array data be logtransformed
#' @param data The list of data objects, either matrix, phyloseq or ExpressionSet objects
#' @return the raw data matrices, samples in the rows
#' @export
#' @examples
#' data(Zhang)
#' matrixList = extractData(list("microbiome" = zhangMicrobio,
#' "metabolome" = zhangMetabo))
extractData = function(data, logTransformGaussian = TRUE){
    datNames = names(data)
    data = lapply(seq_along(data), function(View){
        dat  = extractMat(data[[View]],
                          logTransformGaussian = logTransformGaussian)
        if(is.null(colnames(dat))){
            colnames(dat) = paste0("feat", seq_len(ncol(dat)), "View", View)
        }
        return(dat)
    })
    names(data) = datNames
    data
}
#' A function to extract a data matrix from a number of objects
#' @param Y a phyloseq or eSet object, or another object, or a raw data matrix
#' @param logTransformGaussian A boolean, should array data be logtransformed
#' @param ... additional arguments for the extractor function
#' @return A data matrix with samples in the rows and features in the columns
#' @rdname extractMat
setGeneric("extractMat", function(Y, ...) standardGeneric("extractMat"))
#' @import methods
#' @importFrom phyloseq t otu_table taxa_are_rows
setMethod("extractMat", "phyloseq", function(Y,...){
    if(taxa_are_rows(Y)) t(as(otu_table(Y), "matrix")) else
        as(otu_table(Y), "matrix")
    })
#' @import methods
#' @importFrom phyloseq t taxa_are_rows
setMethod("extractMat", "otu_table", function(Y,...){
    if(taxa_are_rows(Y)) t(as(Y, "matrix")) else as(Y, "matrix")
})
#' @rdname extractMat
#' @import methods
#' @importFrom Biobase exprs
setMethod("extractMat", "ExpressionSet", function(Y, logTransformGaussian,...){
    dat = t(exprs(Y))
    if(logTransformGaussian) dat = log(dat)
    return(dat)
})
#' @rdname extractMat
#' @import methods
#' @importFrom SummarizedExperiment assay
setMethod("extractMat", "SummarizedExperiment", function(Y, ...){
    t(assay(Y))
})
#' @rdname extractMat
setMethod("extractMat", "matrix", function(Y, ...){
    Y
})