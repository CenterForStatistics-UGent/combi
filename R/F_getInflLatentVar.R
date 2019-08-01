#' Extract the influence on the estimation of the latent variable
#' @param score The score matrix
#' @param InvJac The inverse Jacobian
#' @param i the sample index
#'
#' @return The influence of all observations on the i-th latent variable
getInflLatentVar = function(score, InvJac, i){
    score*InvJac[,i]
}