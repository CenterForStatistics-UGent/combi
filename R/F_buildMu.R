#' A function to build the mu matrix
#' @param offSet the offset matrix
#' @param latentVar,paramEsts,distribution Latent variables,
#' parameter estimates and distribution type
#' @param paramMatrix A boolean, are feature parameters provided as matrix
buildMu = function(offSet, latentVar, paramEsts, distribution,
                   paramMatrix = FALSE){
    centralMat  = if(paramMatrix)
        latentVar*paramEsts else
            outer(latentVar, paramEsts)
    switch(distribution,
           "quasi" = offSet*exp(centralMat),
            "gaussian" = offSet + centralMat
    )
}
#