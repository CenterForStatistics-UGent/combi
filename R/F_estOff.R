#' Estimate the column parameters of the independence model
#'
#' @param data a list of data matrices with the same number of samples n in the rows.
#' Also phyloseq objects are acceptable
#' @param distribution a character string describing which distributional assumption should be used.
#' @param rowOff,colOff current values for row and column offsets
#' @param meanVarTrend The estimated mean-variance trend
#' @param col A logical, should column offsets be estimated
#' @param ... passed onto nleqslv
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
            # CompMat = exp(colOff)
            # #Prevent trivial solutions
            # obsLibSizes = rowSums(data, na.rm = TRUE)
            # lowers  = log(obsLibSizes/10)
            # highers = log(obsLibSizes*10)
            # Jacobian = matrix(1,2)
            # # Converting to minimization problem squares the condition number, which is painfully felt.
            # # Still no other solution seems feasible
            # sapply(seq_len(nrow(data)), function(i){
            #     constrOptim.nl(par = rowOff[i], fn = scoreLibsConstr,
            #                    gr = jacLibsConstr,
            #                    Data = data[i,], CompMat = CompMat,
            #                    meanVarTrend = meanVarTrend,
            #                    hin = function(x,...) c(x-lowers[i], highers[i] - x),
            #                    control.outer= list(trace = FALSE),
            #                    hin.jac = function(...) Jacobian)$par
            #     # if(inherits(tmp, "try-error")){
            #     #     foo = 2
            #     # }
            #     # tmp
            #     # spg(par = rowOff[i], fn = scoreLibsConstr, gr = jacLibsConstr,
            #     #     Data = data[i,], CompMat = CompMat,
            #     #     meanVarTrend = meanVarTrend, lower = lowers[i],
            #     #     control = list(trace = FALSE, checkGrad = TRUE))$par
            # })
       }
    } else if(distribution == "binomial"){
    }
}