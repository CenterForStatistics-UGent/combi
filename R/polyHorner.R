#' Horner's method to evaluate a polynomial, copied from the polynom package. the most efficient way
#' @param coefs the polynomial coefficients
#' @param x the input values for the polynomial function
#'
#' @return the evaluated polynomial
polyHorner = function(coefs, x){
    out = 0
    for(coefj in coefs) out = x*out + coefj
    out
}