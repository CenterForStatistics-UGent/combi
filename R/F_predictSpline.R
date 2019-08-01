#' A custom spline prediction function, extending linearly with a slope such that prediction never drops below first bisectant
predictSpline = function(fit, newdata, lowerSlope, lowerEval, linX,
                         coefsQuad, deriv = 0L, meanVarFit, minFit, new.knots,
                         degree){
    out = newdata
    if(meanVarFit == "spline"){
        orderZ = sort.list(newdata)
    } else {
        idLinear = newdata < linX
        idLower = newdata <  minFit
        idLowerNotLinear = idLower & !idLinear
        out[idLowerNotLinear] = if(deriv) coefsQuad[2] + coefsQuad[1]*2*newdata[idLowerNotLinear] else
            polyHorner(coefsQuad, newdata[idLowerNotLinear])
    }
    n = length(newdata)
    #Call C routine directly for speed
 switch(meanVarFit,
                           "spline" = .C("spline_value", new.knots, fit,
                                         length(fit), degree + 1L,
                                         newdata[orderZ], n, deriv,
                                         y = double(n))$y[sort.list(orderZ)],
                            if(deriv) polyHorner(fit[-4]*c(3,2,1), newdata[!idLower])  else
                               polyHorner(fit, newdata[!idLower]))
    # predict(fit, newdata[!idLower], deriv = deriv)[,2]
# tmp2 = cobs:::.splValue(degree = 2, fit$knots, fit$coef,
    # newdata[!idLower][orderZ], deriv = deriv)[order(orderZ)]
}