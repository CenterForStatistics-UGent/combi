#' A small auxiliary function for the indices of the lagrange multipliers
#'
#' @param m an integer, the current dimension
#' @param nLambda1s the number of centering restrictions
#' @param normal a logical, is there a normalization restriction?
#'
#' @return a vector containing the ranks of the current lagrangian multipliers
seq_m = function(y, normal = TRUE, nLambda1s = 1) {
    (y - 1) * (normal + nLambda1s + (y - 2)/2) +
        seq_len(y + nLambda1s - 1 + normal)
}