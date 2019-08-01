#' Estimate the latent variables
#'
#' @param data A list of data matrices
#' @param distributions A character vector describing the distributions
#' @param offsets A list of offset matrices
#' @param paramEsts A list of parameter estimates for the different views
#' @param numVars The number of variables
#' @param latentVars A vector of latent variables
#' @param latentVarsLower The latent variables of the lower dimensions
#' @param n,m The number of samples and dimensions, respectively
#' @param Jac a jacobian matrix, pre-filled where possible
#'
#' @return A vector of length n, the estimates of the latent variables
#' @importFrom nleqslv nleqslv
estLatentVars = function(latentVars, lambdasLatent, constrained, fTol,...){
    out = if(constrained){
        nleqslv(x = c(latentVars, lambdasLatent), fn = derivLagrangianLatentVarsConstr,
                jac = deriv2LagrangianLatentVarsConstr, ...)
    } else {
        nleqslv(x = c(latentVars, lambdasLatent), fn = derivLagrangianLatentVars,
                jac = deriv2LagrangianLatentVars, ...)
    }
    out = list(x = out$x, conv = if(max(abs(out$fvec)) < fTol) 1 else out$termcd)
}