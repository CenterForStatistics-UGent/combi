#' Estimate the latent variables
#'
#' @param lambdasLatent A vector of Lagrange multipliers
#' @param constrained A boolean, is the ordination constrained?
#' @param fTol The convergence tolerance
#' @param ...
#' @param latentVars A vector of latent variables
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