#' Helper function to extract data matrix from phyloseq, expressionset objects etc. Also filers out all zero rows
#'
#' @param logTransformGaussian A boolean, should array data be logtransformed
#' @param data The list of data objects, either matrix, phyloseq or ExpressionSet objects
#'
#' @return the raw data matrices, samples in the rows
#'
#' @importFrom phyloseq taxa_are_rows otu_table t
#' @importFrom Biobase exprs
#' @export
#' @examples
#' data(hmp2)
#' matrixList = extractData(list("microbiome" = microPruneVir,
#' "virome" = virPrune))
extractData = function(data, logTransformGaussian = TRUE){
    datNames = names(data)
    data = lapply(seq_along(data), function(View){
        dat  = data[[View]]
        if(is(dat, "phyloseq") | is(dat, "otu_table")){
    if(taxa_are_rows(dat)){dat = t(dat)}
    dat = as(otu_table(dat), "matrix")
    } else if(is(dat, "ExpressionSet")){
    dat = t(exprs(dat))
    if(logTransformGaussian) dat = log(dat)
    }
        if(is.null(colnames(dat))) colnames(dat) = paste0("feat", seq_len(ncol(dat)), "View", View)
        return(dat)
    })
    names(data) = datNames
    data
}