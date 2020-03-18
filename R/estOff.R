#' Estimate the row/column parameters of the independence model
#'
#' @param data a list of data matrices with the same number of samples n
#' in the rows. Also phyloseq objects are acceptable
#' @param distribution a character string describing which distributional
#' assumption should be used.
#' @param meanVarTrend The estimated mean-variance trend
#' @param col A logical, should column offsets be estimated
#' @param ... passed onto nleqslv
#' @param rowOff,colOff current row and column offset estimates
#' @param newtonRaphson A boolean, should Newton-Raphson be used to solve the
#' estimating equations
#' @param libSizes The library sizes, used to evaluate the mean-variance trend
#'
#' @return The estimated marginal parameters
#'
#' @importFrom BB dfsane
estOff = function(data, distribution, rowOff, colOff, meanVarTrend, col,
                 newtonRaphson, libSizes,...){
    libSizes = exp(rowOff)
    if(distribution == "gaussian"){
        if(col) colMeans(data-rowOff, na.rm = TRUE) else colMeans(t(data)-colOff, na.rm = TRUE)
    } else if(distribution == "quasi"){
        if(col){
            if(newtonRaphson){
                nleqslv(x = if(col) colOff else rowOff, fn = quasiScoreIndep,
                           jac = quasiJacIndep, data = data, otherMargin = if(col) rowOff else colOff,
                            meanVarTrend = meanVarTrend, col = col, libSizes = libSizes, ...)$x
            } else {
                dfsane(par = colOff, fn = quasiScoreIndep,
                           data = data, otherMargin = rowOff,
                           meanVarTrend = meanVarTrend, libSizes = libSizes,
                       col = col)$par #Memory efficient
            }
        } else {
            nleqslv(x = if(col) colOff else rowOff, fn = quasiScoreIndep,
                jac = quasiJacIndep, data = data, otherMargin = if(col) rowOff else colOff,
                meanVarTrend = meanVarTrend, col = col, libSizes = libSizes,...)$x
          }
    } else if(distribution == "binomial"){
    }
}