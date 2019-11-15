#' Perform model-based data integration
#'
#' @param data a list of data matrices with the same number of samples n in the rows.
#' Also phyloseq objects are acceptable
#' @param M the required dimension of the fit, a non-negative integer
#' @param covariates a dataframe of n samples with sample-specific variables.
#' @param distributions a character vector describing which distributional assumption should be used. See details.
#' @param compositional A logical vector with the same length as "data",
#'  indicating if the datasets should be treated as compositional
#' @param maxIt an integer, the maximum number of iterations
#' @param tol A small scalar, the convergence tolerance
#' @param verbose Logical. Should verbose output be printed to the console?
#' @param nCores The number of cores to be used in estimating the feature parameters of each view. See details.
#' @param confounders A dataframe or a list of dataframes with the same length as data.
#'  In the former case the same dataframe is used for conditioning,
#'  In the latter case each view has its own conditioning variables (or NULL).
#' @param minFraction a scalar, each taxon's total abundance
#' should equal at least the number of samples n times minFraction,
#'   otherwise it is trimmed.
#' @param prevCutOff a scalar, the prevalance cutoff for the trimming.
#' @param record A boolean, should intermediate estimates be stored? Can be useful to check convergence
#' @param fTol The tolerance for solving the estimating equations
#' @param logTransformGaussian A boolean, should the array data be logtransformed?
#' @param nleq.control A list of arguments to the nleqslv function
#' @param weights A character string, either 'marginal' or 'uniform', indicating
#' rrhow the feature parameters should be weighted in the normalization
#' @param meanVarFit The type of mean variance fit, see details
#' @param maxFeats The maximal number of features for a Newton-Raphson procedure
#'  to be feasible
#' @param dispFreq An integer, the period after which the variances should be
#' reestimated
#' @param allowMissingness A boolean, should NA values be allowed?
#' @param biasReduction A boolean, should bias reduction be applied to allow for
#' confounder correction in groups with all zeroes? Not guaranteed to work
#' @param maxItFeat,maxItFilt Integers, the maximum allowed number of iterations
#' in the estimation of the feature parametes and confounder parameters
#' respectively
#'
#' @return An object of the "compInt" class, containing all information on the
#' data integration and fitting procedure
#'
#' @details Using more than one core is only implemented on Unix systems.
#' Setting nCores > 1 on Windows will use a single core, with a warning.
#' When the number of cores specified is larger than the number of views,
#' nCores is silently set to the number of views.
#' meanVarFit = "spline" yields a cubic spline fit for the abundance-variance
#'  trend, "cubic" gives a third degree polynomial. Both converge to the
#'  diagonal line with slope 1 for small means.
#'  Distribution can be either "quasi" for quasi likelihood or "gaussian" for Gaussian data
#' @aliases compIntegrate
#' @importFrom limma squeezeVar
#' @importFrom vegan rda
#' @importFrom parallel mclapply
#' @importFrom stats sd
#' @export
#' @examples
#' data(hmp2)
#' microVirDI = compInt(data = list("microbiome" = microPruneVir,
#' "virome" = virPrune), distributions = c("quasi", "quasi"),
#' compositional = c(TRUE, TRUE), verbose = TRUE)
#' microVirDIconstr = compInt(data = list("microbiome" = microPruneVir,
#' "virome" = virPrune), distributions = c("quasi", "quasi"),
#' compositional = c(TRUE, TRUE),
#' covariates = hmp2samVar[, c("diagnosis","biopsy_location","sex")],
#' verbose = TRUE)
compInt = function(data, M = 2L, covariates = NULL, distributions,
                   compositional, maxIt = 3e2L, tol = 1e-3, verbose = FALSE,
                   prevCutOff = 0.95, minFraction = 0.1, logTransformGaussian = TRUE,
                   confounders = NULL, nleq.control = list(maxit = 1e3L, cndtol = 1e-16),
                   record = TRUE, weights = NULL, fTol = 1e-5, nCores = 1,
                   meanVarFit = "spline", maxFeats = 2e3, dispFreq = 10L,
                   allowMissingness = FALSE, biasReduction = TRUE, maxItFeat = 2e1L,
                   maxItFilt = 50L){
    #Switch off multithreading on windows
    if(.Platform$OS.type == "windows" && nCores>1){
        message("Forking not supported on Windows machine!\nUsing only one core.")
        nCores = 1
    }
    nCores = min(nCores, length(data)) # No point in using more cores than views
    #Perform checks
    if(M %in% c(0L,1L) | (as.integer(M)!=M)){
        stop("Please supply non-negative integer dimension of at least 2!")
    }
    #if(length(distributions)==1) distributions = rep(distributions, length(data))
    #if(length(compositional)==1) compositional = rep(compositional, length(data))
    if(!is.logical(compositional)) stop("'compositional' should be a logical vector!")
    if(!is.character(distributions)) stop("'distributions' should be a character vector!")
    if(!all(vapply(FUN.VALUE = integer(1),
                   c(length(data), length(distributions)),
                   identical, length(compositional)))){
        stop("Make sure data, distribution, links and compositional have the same length")
    }
    if(length(data)==1){
        warning("Please provide at least two views for data integration!
                Now fitting an ordination model.",
                immediate. = TRUE)}
    if(is.null(names(data))) names(data) = paste0("View", seq_along(data))
    namesData = names(data)
    #Extract otu table from phyloseq objects6
    data = extractData(data, logTransformGaussian = logTransformGaussian)
    if(any(vapply(FUN.VALUE = TRUE, data, function(x){is.null(rownames(x))}))){
        stop("Make sure to provide sample names for all views!")
    }
    if(any(vapply(FUN.VALUE = TRUE, data, anyNA)) & !allowMissingness){
        stop("Missing data present. To allow fit with missing data, set allowMissingness to TRUE")
    }
    rowNames = lapply(data, function(x){sort(rownames(x))})
    if(!all(vapply(FUN.VALUE = logical(1), rowNames[-1], FUN = identical,
                   rowNames[[1]]))){
        if(allowMissingness){
        # Make sure at least some samples are shared
            if(length(Reduce(rowNames, f = intersect))==0){
                stop("No samples common to all views have been found!
                     Data integration is questionable in this case. Check rownames of objects provided.")
            } else {
                #Fill in NA's for missing samples. These observations do not contribute to the estimating equations
                allRownames = Reduce(rowNames, f = union)
                data  = lapply(data, function(dat){
                    rownamesMissing = allRownames[!allRownames %in% rownames(dat)]
                    naMatrix = matrix(NA, length(rownamesMissing), ncol(dat))
                    rownames(naMatrix) = rownamesMissing
                    dat = rbind(dat, naMatrix)
                    dat[allRownames,]
                })
            }
        } else {
            stop("Rownames do not match!
                 Fix rownames first or set 'allowMissingness' to TRUE to allow for missing data")
            }
    }
    n = nrow(data[[1]]) #Total number of samples
    zeroRows = apply(vapply(data, FUN.VALUE = logical(nrow(data[[1]])),
                            function(x){rowSums(x, na.rm = TRUE)==0}), 1, any)
    if(any(zeroRows)){
        warning("Zero rows\n", paste(which(zeroRows), collapse = " ") ,
                "\nfiltered out prior to fit", immediate. = TRUE)
        data = lapply(data, function(x){x[!zeroRows,]})
        covariates = covariates[!zeroRows,]
        confounders = confounders[!zeroRows,]
    }
    if(!is.null(confounders) && !length(confounders) %in% c(1,length(data))){
        stop("Please provide a single confounder matrix or as many as you provide views!\n")
    }
    #Assign link functions
    links = lapply(distributions, function(distr){
        switch(distr, "gaussian" = "identity", "quasi" = "log",
        "binomial" = "logit")
    })
    #Get the inverse link functions
    invLinks = lapply(links, function(link){
        if(is.function(link)){link = deparse(substitute(link))}
        #Convert to string if needed
        switch(link, "identity" = identity, "log" = exp,
               "logit" = function(x){exp(x)/(1+exp(x))})
    })
    links = lapply(links, match.fun) #Match the link functions
    #Assign weighting schemes
    weights = if(is.null(weights)) ifelse(distributions %in%
                                              c("quasi"), "marginal", "uniform") else weights
    #Define some handy quantities
    numSets = length(data) #Number of views
    if(numSets < nCores){
        message("It has no use supplying more cores than views, only ", numSets, " cores used!")
        nCores = numSets
    }
    seqSets = seq_len(numSets) #A sequence along the views
    #Build confounder matrices
    oneConfMat = length(confounders) == 1 #Single confounder matrix for all views?
    confounders = lapply(seqSets, function(i){if(is.null(confounders[[i]])) NULL else
        data.frame(confounders[[i]])[!zeroRows,,drop = FALSE]}) #Coerce to data frames
    confMats = lapply(confounders, buildConfMat)
    #Build covariate matrix
    constrained = !is.null(covariates)
    # Mean variance trends
    meanVarFit = if(length(meanVarFit)==1){
        rep(meanVarFit, numSets)
    } else if(length(meanVarFit)==numSets){
        meanVarFit
    } else {
        stop("Provide mean-variance trend vector of length one or equal to number of datasets!\n")
    }
    # Prune count data, also using the confounders if needed
    data = lapply(seq_along(data), function(i){
        if(distributions[[i]] == "quasi"){
                data[[i]] = data[[i]][, colMeans(data[[i]]==0, na.rm = TRUE) <= prevCutOff &
                                        colSums(data[[i]], na.rm = TRUE) >= (n * minFraction)]
            if(!oneConfMat && !is.null(confounders[[i]]) && !biasReduction){
            data[[i]] = trimOnConfounders(data = data[[i]], prevCutOff = prevCutOff,
                                       minFraction = minFraction, n = n,
confounders = confMats[[if(length(confounders)>1) i else 1]]$confModelMatTrim)
            }
        }
        return(data[[i]])
    })
    #All zero rows, but not all NAs
    n = nrow(data[[1]])
    zeroRowsIndep = apply(
        vapply(data, FUN.VALUE = logical(n),
               function(x){rowSums(x, na.rm = TRUE)==0 &
                       !apply(x,1, function(x) all(is.na(x)))}), 1, any)
    if(any(zeroRowsIndep)){
        message("Zero rows\n", paste(which(zeroRowsIndep), collapse = " ") ,
                "\nfiltered out after filtering features")
    }
    data = lapply(data, function(x){x[!zeroRowsIndep,]})
    if (constrained) {
        if(!is.data.frame(covariates)){
            stop("Please provide covariate matrix as a dataframe")
        }
        covariates = droplevels(covariates[!zeroRowsIndep,])#Drop unused levels
        tmp = buildCovMat(covariates)
        covMat = tmp$covModelMat
        numCov = ncol(covMat)
        # Already prepare the matrix that defines the
        # equations for centering the coefficients of the
        # dummy variables
        tmp = buildCentMat(tmp$datFrame)
        centMat = tmp$centMat
        covariates = tmp$datFrame
        rm(tmp)

        # Remove rows with NA's, we might want to find
        # something better or leave it to the end user
        if (anyNA(covariates)) {
            NArows = apply(covariates, 1, anyNA)
            if (all(NArows)) {
                stop("All samples have missing covariates")
            }
            data = lapply(data, function(X) X[!NArows, ])
            covariates = covariates[!NArows, ]
            if (!is.null(confounders)) {
                confMats = lapply(confMats,function(confMat){
                    confMat$confModelMatTrim = confMat$confModelMatTrim[!NArows,]
                    confMat$confModelMat = confMat$confModelMat[!NArows,]
                    confMat
                })
            }
            warning(paste("Some covariates contain missing values.
                    We removed samples \n",
                          paste(which(NArows), collapse = ", "),
                          "\n prior to analysis. Imputation of missing values in predictors is left to the user"),
                    immediate. = TRUE)
        }
    } else {
        covModelMat = centMat = NULL
    }
    nLambda1s = if(constrained) nrow(centMat) else 1
    names(data) = namesData
    n = nrow(data[[1]])
    numVars = vapply(FUN.VALUE = integer(1), data, ncol) #Number of variables per view

    IDs = lapply(seqSets, function(i){
        (sum(numVars[seq_len(i-1)])+1):sum(numVars[seq_len(i)])
    }) # The indices of the features of each view
    newtonRaphson = numVars <= maxFeats #Indicates if too many features for newton raphson

    #### MARGINAL INDEPENDENCE MODELS ####
    if(verbose){cat("Estimating independence models...\n")}
    marginModels = mclapply(seqSets, mc.cores = nCores, function(i){
        estIndepModel(data = data[[i]], distribution = distributions[[i]],
                       link = links[[i]], compositional = compositional[[i]],
                       invLink = invLinks[[i]], tol = tol, maxIt = maxIt,
                      meanVarFit = meanVarFit[[i]], newtonRaphson = newtonRaphson[[i]],
                      control = nleq.control, dispFreq = dispFreq)
        })
    #Corresponding offsets
    offsetsMargin = lapply(seqSets, function(i){buildMarginalOffset(indepModel = marginModels[[i]],
                                                      invLink = invLinks[[i]])})
    #Assign weights
    weightsList = lapply(seqSets, function(i){
switch(weights[[i]],
       "uniform" = rep(1, numVars[[i]]),#/numVars[[i]]
       "marginal" = invLinks[[i]](marginModels[[i]]$colOff))
    })
    #Estimate the mean-variance trends for conditioning and starting value calculation
    meanVarTrends = lapply(seqSets, function(i){
        if(distributions[[i]] == "quasi"){
            estMeanVarTrend(data = data[[i]], meanMat = offsetsMargin[[i]],
                            meanVarFit = meanVarFit[[i]],
                            libSizes = invLinks[[i]](marginModels[[i]]$rowOff),
                            baseAbundances = invLinks[[i]](marginModels[[i]]$colOff))
        } else {NULL}
    })

    #### CONDITIONING ####
    if(verbose && !all(vapply(FUN.VALUE = logical(1), confounders, is.null))){
        cat("Conditioning on known confounders ...\n")}
    confVars = lapply(seqSets, function(i){
        if(!oneConfMat && is.null(confounders[[i]])){return(NULL)
            } else {filterConfounders(offSet = offsets[[i]], data = data[[i]],
                              distribution = distributions[[i]], link = links[[i]],
                              compositional = compositional[[i]], invLink = invLinks[[i]],
                              confMat = confMats[[if(length(confounders)>1) i else 1]]$confModelMat,
                              meanVarTrend = meanVarTrends[[i]], numVar = numVars[[i]],
                              control = nleq.control, marginModel = marginModels[[i]],
                              allowMissingness = allowMissingness, biasReduction = biasReduction,
                              maxItFilt = maxItFilt)
            }
        })
    #Prepare a list of independence models
    indepModels = lapply(seqSets, function(i){
        libSizes = invLinks[[i]](marginModels[[i]]$rowOff)
        colMat = matrix(marginModels[[i]]$colOff, n, numVars[[i]], byrow = TRUE) +
                    if(!is.null(confounders[[i]])){confMats[[i]]$confModelMat %*% confVars[[i]]} else 0
        return(list(colMat = colMat, libSizes = libSizes))
            })
    offsets = lapply(seqSets, function(i){
        if(distributions[[i]]=="quasi"){
            expColMat = exp(indepModels[[i]]$colMat)
            return(indepModels[[i]]$libSizes *
                       if(compositional[[i]]) expColMat/rowSums(expColMat) else expColMat)
        } else if(distributions[[i]]=="gaussian"){
            return(indepModels[[i]]$colMat + indepModels[[i]]$libSizes)
        }
    })
    numVars = lapply(offsets, ncol)
    n = nrow(offsets[[1]])
    #### MAIN MODEL ####
    # Pre-allocate some quantities
    iter = rep.int(1L, M)
    converged = logical(M)
    # Lagrange multipliers
    lambdasLatent = if(constrained) numeric((nLambda1s+1)*M + M*(M-1)/2) else
        numeric(M + M*(M-1)/2)
    lambdasParams = lapply(compositional, function(comp){
        integer(M*(2 + (M-1)/2))})
    #M more restrictions for the normalisation
    #Prepare to record the iterations if necessary
    if(record){
        if(constrained) alphaRec = array(0, dim = c(numCov ,M, maxIt)) else
            latentRec = array(0, dim = c(n ,M, maxIt))
        paramRec = lapply(numVars, function(numVar){
            array(0, dim = c(M, numVar, maxIt))
        })
        #Convergence
        latentConv = matrix(nrow = maxIt, ncol =  M)
        paramConv = array(dim = c(maxIt, numSets, M))
    } else {
        latentConv = paramConv = alphaRec = paramRec = latentRec = NULL
    }
    if(allowMissingness){
        naIdList = lapply(data, is.na)
    }

    ### STARTING VALUES ###
    # Find starting values through svd on concatenated datasets, normalized as well as possible
    concat = Reduce(lapply(seqSets, function(i){
        rawResiduals = data[[i]]-offsets[[i]]
        out = switch(distributions[[i]], "gaussian" = rowMultiply(rawResiduals, 1/apply(rawResiduals, 2, sd)),
               "quasi" = rawResiduals/sqrt(meanVarTrends[[i]](invLinks[[i]](marginModels[[i]]$colOff),
                                                              invLinks[[i]](marginModels[[i]]$rowOff))))
        if(allowMissingness){
            out[naIdList[[i]]] = 0
        }
        out
    }), f = cbind)
    if(constrained){
        CCA = rda(X = concat, Y = covMat)$CCA
        # Redundancy analysis for starting values
        if (sum(!colnames(covMat) %in% CCA$alias) < M) {
            M = sum(!colnames(covMat) %in% CCA$alias)
            warning(immediate. = TRUE, paste("Can only fit an ordination with",
                                             M, "dimensions with so few covariates!"))
        }
        alphas = matrix(0, numCov, M)
        alphas[!colnames(covMat) %in% CCA$alias,
              ] = CCA$biplot[, seq_len(M)]
        alphas = t(t(alphas) - colMeans(alphas))
        #alphas = t(t(alphas)/sqrt(colSums(alphas^2)))
        latentVars = covMat %*% alphas
        paramEsts = lapply(seqSets, function(i){
            t(CCA$v[IDs[[i]],seq_len(M)])
        })
    } else {
    Svd = svd(concat)
    #Re-center parameter estimates
    tmpMats = lapply(seqSets, function(i){
        tmpMat = Svd$v[IDs[[i]], seq_len(M)]
        tmpMat = t(tmpMat) - colMeans(tmpMat*weightsList[[i]])
        tmpMat
      })
    #Transfer some weight to fit the constraints
    tmpWeights = vapply(seqSets, FUN.VALUE = double(M), function(i){
        sqrt(colSums(t(tmpMats[[i]])^2*weightsList[[i]]))
    })
    latentVars = scale((Svd$u)[,seq_len(M)], scale = FALSE,
                       center = TRUE)
    # %*% diag(Svd$d)
    #latentVars = t(t(latentVars)*exp(rowMeans(log(tmpWeights))))
    paramEsts =  lapply(seqSets, function(i){
        tmpMats[[i]]/tmpWeights[,i]
    })
    rm(Svd)
    }
    rm(concat)
    meanVarTrends = vector("list", M)
    #Estimate the mean-variance trend once, given the complete indendence model

    for(m in seq_len(M)){
        if(verbose) cat("Estimating dimension", m, "\n")
        #Modify the offset if needed
        if(m>1){
            offsets = lapply(seqSets, function(i){
                if(compositional[[i]]){
                  indepModels[[i]]$libSizes*
                        buildCompMat(indepModels[[i]]$colMat, paramEsts[[i]], latentVars, m-1)
                } else {
                buildMu(offSet = offsets[[i]], latentVar = latentVars[,m-1],
                        paramEsts = paramEsts[[i]][m-1,], distribution = distributions[[i]])
                }
            })
        }
        #Prepare the jacobians
                ## Latent variables
        emptyJacLatent = buildEmptyJac(n = if(constrained) numCov else n, m = m,
                          lower = if(constrained) alphas else latentVars, nLambda1s = nLambda1s,
                          normal = constrained, centMat = centMat)
                ## Feature parameters
        emptyFeatureJacs = lapply(seqSets, function(i){
            buildEmptyJac(numVars[[i]], m = m, lower = t(paramEsts[[i]]),
                          normal = compositional[[i]], weights = weightsList[[i]],
                          distribution = distributions[[i]])})
        while(iter[[m]] <= maxIt && !converged[[m]]){

            #Estimate mean-variance trend
            if(((iter[[m]]-1) %% dispFreq)==0){
                if(verbose) cat("Estimating mean variance trends\n")
            meanVarTrends[[m]] = lapply(seqSets, function(i){
                if(distributions[[i]] == "quasi"){
                    estMeanVarTrend(data = data[[i]], meanMat = offsetsMargin[[i]],
                                    meanVarFit = meanVarFit[[i]],
                                    libSizes = indepModels[[i]]$libSizes,
                                    baseAbundances = invLinks[[i]](marginModels[[i]]$colOff))
                } else {NULL}
            })
            }
            #Store old estimates
            if(constrained){
                alphaOld = alphas[,m]
            } else {latentVarsOld = latentVars[,m]}
        paramEstsOld = lapply(paramEsts, function(x){x[m,]})

         #Estimate parameters
        if(verbose) cat("Estimating feature parameters of dimension", m,
                        ": iteration", iter[[m]], "\n")
        paramEstsTmp = estFeatureParameters(data = data, distributions = distributions,
                                            offsets = offsets, paramEsts = paramEsts,
                                            latentVars = latentVars, m = m,
                                            numVars = numVars, meanVarTrends = meanVarTrends[[m]],
                                            lambdasParams = lambdasParams,
                                            JacFeatures = emptyFeatureJacs, seqSets = seqSets,
                                            control = nleq.control, nCores = nCores,
                                            weights = weightsList, compositional = compositional,
                                            indepModels = indepModels, fTol = fTol,
                                            newtonRaphson = newtonRaphson, allowMissingness = allowMissingness,
                                            maxItFeat = maxItFeat)
        for (i in seqSets){
            #Enforce restrictions to stabilize algorithm
            tmp = paramEstsTmp[[i]]$x[seq_len(numVars[[i]])]
            tmp = tmp - sum(tmp*weightsList[[i]]) #Center
            if(m>1){
                paramEsts[[i]][m,] = GramSchmidtOrth(paramEsts[[i]][m,],
                                                     paramEsts[[i]][m-1,],
                                                     weights = weightsList[[i]],
                                                     norm = compositional[[i]])
            } #Orthogonalize
            paramEsts[[i]][m,] = tmp/sqrt(sum(tmp^2*weightsList[[i]])) #Scale
            lambdasParams[[i]][seq_m(m, nLambda1s = 1, normal = compositional[[i]] )] =
                paramEstsTmp[[i]]$x[-seq_len(numVars[[i]])] #+ compositional[[i]]
            #Record convergence
            if(record) paramConv[iter[[m]],i,m] = paramEstsTmp[[i]]$conv
        };rm(tmp)
        ### Estimate nuissance parameters ###
        #if(verbose) cat("Estimating nuissance parameters ...\n")
        #For each dimension, re-estimate the posterior standard deviation using limma, if applicable
        varPosts = lapply(seqSets, function(i){
            if(distributions[[i]] == "gaussian"){
                varEst = colSums(na.rm = TRUE, (data[[i]] - offsets[[i]] -
                                      outer(latentVars[,m], paramEsts[[i]][m,]))^2)/
                    (n-m-1)
                squeezeVar(var = varEst, df = n-m-1)$var.post
            } else {NULL}
        })
        ### Estimate latent variables ###
        if(verbose) cat("Estimating latent variables of dimension", m,
                        ": iteration", iter[[m]], "\n")
        #Prefab the parameter matrices
        paramMats = lapply(seqSets, function(i){
            matrix(paramEsts[[i]][m,], byrow = TRUE, n, numVars[[i]])
            })
        #nleq.control$trace = TRUE
        latentVarsTmp = estLatentVars(data = data,
                            distributions = distributions,
                            offsets = offsets, paramEsts = paramEsts, paramMats = paramMats,
                            latentVars = (if(constrained) alphas else latentVars)[, m],
                            latentVarsLower = (if(constrained) alphas else latentVars)[,
                                                seq_len(m-1), drop = FALSE],
                            n = n, m = m, numSets = numSets, numVars = numVars,
                            meanVarTrends = meanVarTrends[[m]], links = links,
                            lambdasLatent = lambdasLatent[
                                          seq_m(m, normal = constrained, nLambda1s = nLambda1s)],
                            Jac = as.matrix(emptyJacLatent),
                            control = nleq.control, covMat = covMat,
                            numCov = numCov, centMat = centMat, nLambda1s = nLambda1s,
                            varPosts = varPosts, constrained = constrained,
                            compositional = compositional, indepModels = indepModels,
                            fTol = fTol, allowMissingness = allowMissingness)
        #nleq.control$trace = FALSE
        lambdasLatent[seq_m(m, normal = constrained, nLambda1s = nLambda1s)] =
            latentVarsTmp$x[-seq_len(if(constrained) numCov else n)]
        if(constrained){
            alphas[,m] = latentVarsTmp$x[seq_len(if(constrained) numCov else n)]
            if(m>1){
            alphas[,m] = GramSchmidtOrth(alphas[,m], alphas[,m-1], norm = TRUE)
            }
            latentVars[,m] = covMat %*% alphas[,m]
        } else {
            latentVars[, m] = latentVarsTmp$x[seq_len(n)] -
                mean(latentVarsTmp$x[seq_len(n)]) #Enforce centering for stability
            if(m>1){
                latentVars[, m] = GramSchmidtOrth(latentVars[, m], latentVars[, m-1], norm = FALSE)
            } #Orthogonalize
        }
        if(record) latentConv[iter[[m]],m] = latentVarsTmp$conv #Store convergence
        # Store intermediates if necessary
        if(record){
            if(constrained){
        alphaRec[,m,iter[[m]]] = alphas[,m]
            } else {
        latentRec[,m,iter[[m]]] = latentVars[,m]
            }
        for (i in seqSets){
            paramRec[[i]][m,,iter[[m]]] = paramEsts[[i]][m,]
        }
        }
        # Change iterator
        iter[[m]] = iter[[m]] + 1
        # Check for convergence
        converged[[m]] = all(vapply(seqSets, FUN.VALUE = logical(1),
                                    function(i){
            sqrt(mean((1-paramEsts[[i]][m,]/paramEstsOld[[i]])^2)) < tol})) &&
            (if(constrained) sqrt(mean((1-alphas[,m]/alphaOld)^2)) else
                sqrt(mean((1-latentVars[,m]/latentVarsOld)^2))) < tol
        }
    iter[[m]] = iter[[m]]-1
    }
    if(!all(converged)){
        warning("Some dimensions did not converge, please investigate cause!")
    }
    ### Assign names
    rownames(latentVars) =  rownames(data[[1]])
    colnames(latentVars) = paste0("Dim", seq_len(M))
    if(constrained) {
        colnames(alphas) = colnames(latentVars)
        rownames(alphas) = colnames(covMat)
    }
    for(i in seqSets){
        colnames(paramEsts[[i]]) = colnames(data[[i]])
        rownames(paramEsts[[i]]) = colnames(latentVars)
    }
    #### RETURN RESULT ####
    out = list(latentVars = latentVars, paramEsts = paramEsts,
               indepModels = indepModels, lambdasLatent = lambdasLatent,
               lambdasParams = lambdasParams, iter = iter, converged = converged,
               alphas = if(constrained) alphas else NULL,
               covMat = if(constrained) covMat else NULL, covariates = covariates,
               weights = weightsList, varPosts = varPosts, latentConv = latentConv,
               paramConv = paramConv, data = data, compositional = compositional,
               meanVarTrends = meanVarTrends, distributions = distributions,
               confVars = confVars, confMats = confMats, marginModels = marginModels,
               newtonRaphson = newtonRaphson, allowMissingness = allowMissingness,
               centMat = centMat)
    if(record){
        out$paramRec = lapply(seqSets, function(i){
            tmp = paramRec[[i]][,,seq_len(max(iter))]
            colnames(tmp) = colnames(paramEsts[[i]])
            tmp})
        out$paramConv = out$paramConv[seq_len(max(iter)),,]
        out$latentConv = out$latentConv[seq_len(max(iter)),]
        if(constrained){
        out$alphaRec = alphaRec[,,seq_len(max(iter))]
        rownames(out$alphaRec) = rownames(alphas)
        } else {
        out$latentRec = latentRec[,,seq_len(max(iter))]
        rownames(latentRec) = rownames(latentVars)
        }
    }
    class(out) = "compInt"
    return(out)
}