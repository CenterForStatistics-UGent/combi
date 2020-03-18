#' A custom spline prediction function, extending linearly with a slope such that prediction never drops below first bisectant
#' @param fit The existing spline fit
#' @param newdata points in which the spline needs to be evaluated
#' @param linX The x at which the fit becomes linear and intersects the diagonal
#'  line
#' @param coefsQuad parameters of a quadratic fit
#' @param deriv An integer. Which derivative is required?
#' @param meanVarFit A character string, indicating which type of mean variance
#' fit is being used
#' @param minFit The lower bound of the cubic fit
#' @param new.knots The knots at which the spline is to be evaluated
#' @param degree The degree of the polynomial fit
#'
#'@return The evaluation of the spline, i.e. the predicted variance
predictSpline = function(fit, newdata, linX, coefsQuad, deriv = 0L,
                         meanVarFit, minFit, new.knots, degree){
    n = length(newdata)
    if(meanVarFit == "spline"){
        orderZ = sort.list(newdata)
        #Call C routine directly for speed
        .C("spline_value", PACKAGE = "cobs", new.knots, fit,
                 length(fit), degree + 1L,
                 newdata[orderZ], n, deriv,
                 y = double(n))$y[sort.list(orderZ)]
    } else {
        out = newdata
        idLinear = newdata < linX
        idLower = newdata < minFit
        idLowerNotLinear = idLower & !idLinear
        out[idLowerNotLinear] = if(deriv) coefsQuad[2] +
            coefsQuad[1]*2*newdata[idLowerNotLinear] else
            polyHorner(coefsQuad, newdata[idLowerNotLinear])
        out[!idLower] = if(deriv) polyHorner(fit[-4]*c(3,2,1), newdata[!idLower])  else
            polyHorner(fit, newdata[!idLower])
        return(out)
    }

}