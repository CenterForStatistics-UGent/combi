#' Build an offset matrix from an marginal model object
#'
#' @param indepModel The fitted marginal model, a list
#' @param invLink The inverse link function
#'
#' @return an offset matrix of the size of the data
buildMarginalOffset = function(indepModel, invLink){
    invLink(outer(indepModel$rowOff, indepModel$colOff, "+"))
}