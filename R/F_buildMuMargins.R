#' Build the marginal mu matrix
#'
#' @param x The marginal parameters begin estimated
#' @param otherMargin The parameters of the other margin
#' @param col A logical, are the column parameters being estimated?
#'
#' @return a matrix of means
buildMuMargins = function(x, otherMargin, col){
    if(col) outer(otherMargin, x, "+") else outer(x, otherMargin, "+")
    # -log(sum(exp(x)))  -log(sum(exp(otherMargin)))
}