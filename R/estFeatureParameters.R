#' Estimate the feature parameters
#'
#' @param data A list of data matrices
#' @param distributions A character vector describing the distributions
#' @param offsets A list of offset matrices
#' @param paramEsts Current list of parameter estimates for the different views
#' @param numVars The number of variables
#' @param lambdasParams The lagrange multipliers
#' @param seqSets A vector with view indices
#' @param nCores The number of cores to use in multithreading
#' @param m The dimension
#' @param JacFeatures An empty Jacobian matrix
#' @param meanVarTrends The mean-variance trends of the different views
#' @param control A list of control arguments for the nleqslv function
#' @param weights The normalization weights
#' @param compositional A list of booleans indicating compositionality
#' @param indepModels A list of independence model
#' @param fTol A convergence tolerance
#' @param allowMissingness A boolean indicating whether missing values are
#' allowed
#' @param maxItFeat An integer, the maximum number of iterations
#' @param ... Additional arguments passed on to the score and jacobian functions
#' @param latentVars A vector of latent variables
#'
#' @importFrom nleqslv nleqslv
#' @importFrom BB BBsolve
#' @importFrom stats rnorm
#'
#' @return A vector with estimates of the feature parameters
#'
#' @details If forking is available on the OS and the number of cores specified
#'   is larger than 1, this function will multithread the estimation of the
#'   feature parameters for the different views.
estFeatureParameters = function(paramEsts, lambdasParams, seqSets, data,
                                distributions, offsets, nCores, m, JacFeatures,
                                meanVarTrends, latentVars, numVars, control,
                                weights, compositional, indepModels, fTol,
                                allowMissingness, maxItFeat,...){
    bplapply(seqSets, function(i){
        if(distributions[[i]] == "gaussian"){
            #Solve system of linear equations, variances still not needed
            diag(JacFeatures[[i]])[seq_len(numVars[[i]])] = sum(latentVars[,m]^2)
            x = c(colSums(na.rm = TRUE, latentVars[,m]*(data[[i]] - offsets[[i]])) ,integer(m))
            sol = solve(JacFeatures[[i]], x)
            return(list(x = sol@x, conv = 1))
        } else {
        out = try(nleqslv(x = c(paramEsts[[i]][m,], lambdasParams[[i]][
            seqM(m, normal = compositional[[i]])]),
            fn = derivLagrangianFeatures, jac = deriv2LagrangianFeatures,
            distribution = distributions[[i]], numVar = numVars[[i]],
            paramEstsLower = paramEsts[[i]][seq_len(m-1),,drop = FALSE],
            Jac = as.matrix(JacFeatures[[i]]), meanVarTrend = meanVarTrends[[i]],
            data = data[[i]], latentVars = latentVars, mm = m,
            control = control, weights = weights[[i]], offSet = offsets[[i]],
            compositional = compositional[[i]], indepModel = indepModels[[i]],
            allowMissingness = allowMissingness, ...),
           silent = TRUE)
        if(!inherits(out, "try-error")){
        conv = if(max(abs(out$fvec)) < fTol) 1 else out$termcd
        #Introduce some randomness to avoid saddle points, local extrema. Not ideal but happens so often (e.g. MCMC)
        sol = out$x
        solf = sum(out$fvec^2)
        iter = 0
        id  = seq_along(paramEsts[[i]][m,])
        while((conv != 1L) && ((iter = iter + 1) <= maxItFeat)){
            boo = rnorm(length(out$x), sd = 1.1^iter)
            boo[id] = boo[id] - sum(boo[id]*weights[[i]])
            for(jj in seq_len(m-1)){
                boo[id] =
                gramSchmidtOrth(boo[id], paramEsts[[i]][jj,],
                                weights[[i]], norm = FALSE)
            }
            boo[id] = boo[id]/sqrt(sum(boo[id]^2*weights[[i]]))
            out = nleqslv(x = boo,
                fn = derivLagrangianFeatures, jac = deriv2LagrangianFeatures,
                distribution = distributions[[i]], numVar = numVars[[i]],
                paramEstsLower = paramEsts[[i]][seq_len(m-1),,drop = FALSE],
                Jac = as.matrix(JacFeatures[[i]]), meanVarTrend = meanVarTrends[[i]],
                data = data[[i]], latentVars = latentVars, mm = m,
                control = control, weights = weights[[i]], offSet = offsets[[i]],
                compositional = compositional[[i]], indepModel = indepModels[[i]],
                allowMissingness = allowMissingness, ...)
            conv = if(max(abs(out$fvec)) < fTol) 1L else out$termcd
            if(sum(out$fvec^2) < solf){
                sol = out$x; solf = sum(out$fvec^2)
            }
        }
        return(list(x = sol, conv = conv))
        } else {
        out = BBsolve(par = c(paramEsts[[i]][m,], lambdasParams[[i]][
            seqM(m, normal = compositional[[i]])]),
            fn = derivLagrangianFeatures, distribution = distributions[[i]], numVar = numVars[[i]],
            paramEstsLower = paramEsts[[i]][seq_len(m-1),,drop = FALSE],
             meanVarTrend = meanVarTrends[[i]], Jac =NULL,
            data = data[[i]], latentVars = latentVars, mm = m,
            control = list(), weights = weights[[i]], offSet = offsets[[i]],
            compositional = compositional[[i]], indepModel = indepModels[[i]],
            allowMissingness = allowMissingness,...)
        conv = switch(as.character(out$convergence),
                      "0" = 1, "1" = 4, "2" = 3, 4) # Same error codes as nleqslv
        return(list(x = out$par, conv = conv))
}

        }
    })
}