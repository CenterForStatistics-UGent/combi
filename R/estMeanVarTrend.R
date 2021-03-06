#' Estimate a column-wise mean-variance trend
#'
#' @param data the data matrix with n rows
#' @param meanMat the estimated mean matrix
#' @param plot A boolean, should the trend be plotted?
#' @param baseAbundances The baseline abundances
#' @param libSizes Library sizes
#' @param meanVarFit A character string describing the type of trend to be
#'  fitted: either "spline" or "cubic"
#' @param degree The degree of the spline
#' @param constraint Constraint to the spline
#' @param ... additional arguments passed on to the plot() function
#'
#' @import cobs stats graphics
#' @importFrom alabama constrOptim.nl
#'
#' @return A list with components
#' \item{meanVarTrend}{An smoothed trend function, that can map a mean on a variance}
#' \item{meanVarTrendDeriv}{A derivative function of this}
estMeanVarTrend = function(data, meanMat, baseAbundances, libSizes,
                           plot = FALSE, meanVarFit, degree = 2L,
                           constraint = "none",...){
    # A vector of variances
    varVec = colSums((meanMat - data)^2/libSizes, na.rm = TRUE)/(nrow(data)-1)
    #Logarithms of mean vector
    logMeanVec = log(baseAbundances); logVarVec = log(varVec)
    #The smallest value
    minFit = min(logMeanVec)
    if(meanVarFit == "spline"){
        #Slope of linear fit
        slopeLin = lm(logVarVec ~ logMeanVec)$coef[2]
        #Find linearization point by minimizing residuals.
        #Is quite computation intensive, so go for a rough approximation
        minPos = optimize(f = function(minPos){
            suppressWarnings(sum(cobs(x = logMeanVec, y = logVarVec, constraint = constraint,
                       pointwise = rbind(c(0, minPos, minPos),
                                         c(2,minPos, 1)),
                       lambda = -1,
                    print.mesg = FALSE, print.warn = FALSE, maxiter = 25,
                    lambda.length = 10, keep.x.ps = FALSE, degree = degree)$resid^2))
        }, interval = c(log(1/sum(meanMat))-10, minFit + 1),
        tol = 0.1)$minimum
    #Use optimal value
    pointwiseMat = rbind(c(0, minPos, minPos),
                         c(2,minPos, 1))
    #Enforce slope 1 and y equal to x at small value + slope 1 at upper end
    mvFit = suppressWarnings(cobs(x = logMeanVec, y = logVarVec, constraint = constraint,
         pointwise = pointwiseMat, lambda = -1, degree = degree,
         print.mesg = FALSE, print.warn = FALSE, keep.x.ps = FALSE))
    linX = minPos
    if (length(mvFit$coef) > (nk1 = length(mvFit$knots) + 1L))
        length(mvFit$coef) = nk1
    new.knots <- c(rep.int(mvFit$knots[1], degree), mvFit$knots,
                   rep.int(mvFit$knots[length(mvFit$knots)], degree))
    mvFit = mvFit[["coef"]]
    } else if(meanVarFit == "cubic"){
        #Enforce slope >1 and function value > x at smallest observed mean
        mvFit = lm(logVarVec ~ logMeanVec + I(logMeanVec^2) + I(logMeanVec^3))$coef
        lowerEval = c(c(1, minFit, minFit^2, minFit^3) %*% mvFit)
        lowerSlope = c(c(1, 2*minFit, 3*minFit^2) %*% mvFit[-1])
        mvFit = mvFit[4:1] #Invert order for Horner's rule
    if(((lowerSlope-1) * (lowerEval-minFit)) >0){
    #If curve points towards diagonal Poisson line, find quadratic fit that has
    #1)The same slope
    #2)The same value at the end of the spline, and
    #3) has the diagonal line as a tangent
    quadFunc = function(x){
        a = x[1]; b = x[2]; cc = x[3]
        c(deriv = minFit*2*a+ b - lowerSlope,
          cont = minFit^2*a + minFit*b + cc - lowerEval,
          discrim = (b-1)^2-4*a*cc)
    }
    coefsQuad = nleqslv(rep(1,3), fn = quadFunc)$x
    linX = (1-coefsQuad[2])/(2*coefsQuad[1]) #The x where parabole has slope 1 and touches the diagonal
    } else {
    # If curve points away from the diagonal, backtrack until the slope was equal to 1 and linearize from there (still not optimal)
    minFit = nleqslv(minFit, fn = function(x){
        switch(meanVarFit,
           "cubic" = c(1, 2*x, 3*x^2) %*% mvFit[-1],
           "spline" = predict(mvFit, x, deriv = 1L)$y)-1})$x
    linX = minFit
    coefsQuad = integer(3)
    }
    } else if(meanVarFit == "quadratic") {
    coefsQuad = constrOptim.nl(par = rep.int(1L,3),
               heq = function(pars){(pars[2]-1)^2-4*pars[1]*pars[3]},
               fn = function(pars){sum((polyHorner(pars, logMeanVec)-logVarVec)^2)},
               control.outer = list(trace = FALSE))$par
    linX = (1-coefsQuad[2])/(2*coefsQuad[1])
    minFit = Inf; mvFit = NULL
    } else {}
    dispFunction = function(means, libSizes, deriv = 0L, outerProd = !deriv){
        logMeans = log(means)
        idSmaller = logMeans < linX
        baa = means
        if(!deriv){
            baa[!idSmaller] = exp(predictSpline(mvFit, logMeans[!idSmaller],
                              linX, coefsQuad, meanVarFit = meanVarFit,
                              minFit = minFit, new.knots = new.knots,
                              degree = degree))
        } else{
            baa[!idSmaller] = predictSpline(mvFit, logMeans[!idSmaller], linX,
                         coefsQuad, deriv = deriv, meanVarFit = meanVarFit,
                         minFit = minFit, new.knots = new.knots, degree = degree)*
                exp(predictSpline(mvFit, logMeans[!idSmaller],
                                  linX, coefsQuad, meanVarFit = meanVarFit,
                                  minFit = minFit, new.knots = new.knots,
                                  degree = degree))/means[!idSmaller]
            baa[idSmaller] = 1
        }
        if(outerProd) outer(libSizes, baa) else baa
    }
    if(plot){
    plot(x = logMeanVec, y = logVarVec, xlab = "log(relative abundance)",
         ylab = "log(variance/library size)", xlim = c(linX-1, max(logMeanVec)),
         ylim = c(linX-1, max(logVarVec)), col = "grey", cex = 0.5,...)
    abline(0, 1, lty ="dotted", col = "red")
    Seq = seq(linX-2, max(logMeanVec), 0.01)
    idNum = (Seq < minFit) + (Seq < linX)
    if(meanVarFit == "spline") idNum[idNum==1] = 0
    eval = log(dispFunction(exp(Seq), outerProd = FALSE))
    lines(Seq[idNum==0], eval[idNum==0], col = "black")
    lines(Seq[idNum==1], eval[idNum==1], col = "blue")
    lines(Seq[idNum==2], eval[idNum==2], col = "orange")
    legend("topleft", legend = c(switch(meanVarFit, "spline"= "Spline",
                                        "cubic" = c("Cubic fit", "Quadratic part")),
                                 "Linear part"),
           lty = 1, col  = c("black",
                             switch(meanVarFit, "spline"= NULL, "cubic" = "blue"), "orange"))
    }
    return(dispFunction)
}