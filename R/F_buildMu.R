#' A function to build the mu matrix
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