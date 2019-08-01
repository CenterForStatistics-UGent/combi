#' Estimate the independence model belonging to one view
#'
#' @param data a list of data matrices with the same number of samples n in the rows.
#' Also phyloseq objects are acceptable
#' @param distribution a character string describing which distributional assumption should be used.
#' @param compositional A logical indicating if the dataset should be treated as compositional
#' @param maxIt an integer, the maximum number of iterations
#' @param tol A small scalar, the convergence tolerance
#'
#' @return A list with elements
#' \item{rowOff}{The row offsets}
#' \item{colOff}{The column offsets}
#' \item{converged}{A logical flag, indicating whether the fit converged}
#' \item{iter}{An integer, the number of iterations}
#' @importFrom nleqslv nleqslv
estIndepModel = function(data, distribution, compositional, maxIt, tol, link,
                         invLink, meanVarFit, newtonRaphson, dispFreq,...){
    converged = FALSE
    iter = 1L
    #Disregard NA's
    NArows = apply(data, 1, function(x) all(is.na(x)))
    #Create starting values
    rowOff = link(rowSums(data, na.rm = TRUE))
    colOff = link(colSums(data, na.rm = TRUE)/sum(data, na.rm = TRUE))
    while(iter < maxIt && !converged){
        #Estimate the mean-variance trend if needed
        if(((iter-1) %% dispFreq == 0) && distribution == "quasi"){
            #Mean variance trend
            meanVarTrend = suppressWarnings(estMeanVarTrend(data = data,
                                           meanMat = invLink(outer(rowOff,
                                                                   colOff, "+")),
                                           meanVarFit = meanVarFit,
                                           libSizes = invLink(rowOff),
                                           baseAbundances = invLink(colOff)))
        }
        rowOffOld = rowOff; colOffOld = colOff

        rowOff[!NArows] = estOff(colOff = colOff, data = data[!NArows,], distribution = distribution,
                         col = FALSE, rowOff = rowOff[!NArows], meanVarTrend = meanVarTrend,
                        newtonRaphson = newtonRaphson, libSizes = invLink(rowOff),...)
        if(distribution == "quasi"){
            rowOff[!NArows][rowOff[!NArows] < link(1)] = link(1)
        }
        colOff = estOff(colOff = colOff, data = data[!NArows,], distribution = distribution,
                        rowOff = rowOff[!NArows], meanVarTrend = meanVarTrend, col = TRUE,
                        newtonRaphson = newtonRaphson, libSizes = invLink(rowOff),...)
        if(distribution == "quasi"){
            colOff[colOff < link(1e-8)] = link(1e-8)
        }
        #Enforce compositionality here, avoids trouble with the jacobian,
        #amongst others
        if(compositional) {
            if(distribution == "quasi"){
                colOff = colOff - log(sum(exp(colOff)))
            } else if(distribution == "gaussian"){
                colOff = colOff/sum(colOff)
            }
        }
        iter = iter + 1L
        converged = sqrt(mean((rowOff[!NArows] - rowOffOld[!NArows])^2)) < tol &&
            sqrt(mean((colOff - colOffOld)^2)) < tol
    }
    if(!converged){
        warning(immediate. = TRUE,
                "Fit of independence model did not converge. Please investigate cause.\n")
    }
    rowOff[NArows] = 0
    return(list(rowOff = rowOff, colOff = colOff, converged = converged, iter = iter))
    }