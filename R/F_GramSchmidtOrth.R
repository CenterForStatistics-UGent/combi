#' Gram schimdt orhtogonalize a with respect to b, and normalize
#'
#' @param a the vector to be orthogonalized
#' @param weights weights vector
#' @param norm a boolean, should the result be normalized?
#' @param b the vector to be orthogonalized to
#'
#' @return The orthogonalized vector
GramSchmidtOrth = function(a, b, weights = 1, norm = TRUE){
    tmp = a-c(crossprod(a*weights, b))/sum(b^2*weights)*b
    if(norm) tmp/sqrt(sum(tmp^2*weights)) else tmp
}