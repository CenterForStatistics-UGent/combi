#' Build the composition matrix for a certain dimension m dimensions
#' @param colMat The nxp independence model composition matrix
#' @param paramEsts The matrix of feature parameter estimates
#' @param latentVar The matrix of latent variables
#' @param m the required dimension
#' @param norm a boolean, should the composition matrix be normalized?
#' @param id The vector of dimensions to consider
#' @param subtractMax A boolean, should the maximum be substracted from every
#' composition prior to exponentiation? Recommended for numerical stability
#'
#' @return A matrix with compositions in the rows
buildCompMat = function(colMat, paramEsts, latentVar, m, norm = TRUE,
                        id = seq_len(m), subtractMax = TRUE){
    tmp0 = colMat + (latentVar[,id, drop = FALSE] %*%
                     paramEsts[id,, drop = FALSE])
    #tmp[is.infinite(tmp)] = 1e16
    if(subtractMax) tmp0 = tmp0 - apply(tmp0, 1, max)
    tmp = exp(tmp0)
    if(norm) tmp/rowSums(tmp) else tmp
}