#' Array multiplication
#' @param centralMat an nxp matrix
#' @param outerMat an nxd matrix
#' @param ncols an integer, the number of columns of outerMat
#'
#' @return an nxpxd matrix, the stacked matrices of centralMat multiplied to every column of outermat
arrayMult = function(centralMat, outerMat, ncols = ncol(outerMat)){
    vapply(seq_len(ncols),
           FUN.VALUE = centralMat, function(x) {
               outerMat[, x] * centralMat
           })
}