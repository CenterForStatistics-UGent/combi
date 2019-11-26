#' A ggplot line plot showing the influences
#' @importFrom reshape2 melt
#' @importFrom stats aggregate
#'
#' @param modelObj The fitted data integration
#' @param plotType The type of plot requested, see details
#' @param pointFun The function to calculate the summary measure to be plotted
#' @param lineSize The line size
#' @param Dim The dimension required
#' @param samples index vector of which samples to be plotted
#' @param ... additional arguments passed on to the influence() function
#'
#' @details The options for plotType are: "pointPlot": Dot plot of total influence
#'  per view and sample, "boxplot": plot boxplot of influence of all
#'   observations per view and sample, "boxplotSingle": boxplot of log absolute
#'    total influence per view, "lineplot": line plot of total influence
#'  per view and sample. In the pointplot, dots crosses represent parameter estimates
#'
#' @return A ggplot object
#' @export
#' @examples
#' data(Zhang)
#' microMetaboInt = combi(
#' list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
#' distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
#' logTransformGaussian = FALSE, verbose = TRUE)
#' inlfPlot(microMetaboInt)
#' #Constrained
#' microMetaboIntConstr = combi(
#'     list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
#'     distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
#'     logTransformGaussian = FALSE, covariates = zhangMetavars, verbose = TRUE)
#' inflPlot(microMetaboIntConstr)
inflPlot = function(modelObj, plotType = ifelse(length(modelObj$data) <= 2,
                                                "pointplot", "boxplot"),
                    pointFun = "sum", lineSize = 0.07, Dim = 1,
                    samples = seq_len(nrow(if(is.null(modelObj$covariates))
                        modelObj$latentVars else modelObj$alphas)), ...){
    #if(length(modelObj$data) >2) "sum" else "mean"
    if(!plotType %in% c("lineplot", "pointplot", "boxplot", "boxplotSingle")){
        stop("plotType not recognized, see details!")
    }
    constrained = !is.null(modelObj$covariates)
    inflObj = influence(modelObj, Dim = Dim,...)
    inflMat = as.matrix(with(inflObj, InvJac %*% score))[samples,,drop = FALSE]
    numVars = vapply(FUN.VALUE = integer(1), modelObj$data, ncol)
    IDs = lapply(seq_along(numVars), function(i){
        (sum(numVars[seq_len(i-1)])+1):sum(numVars[seq_len(i)])
    })
    rownames(inflMat) = (if(constrained) rownames(modelObj$alphas) else
        rownames(modelObj$latentVars))[samples]
    moltInflMat = melt(as.matrix(inflMat))
    names(moltInflMat) = c("LatentVariable", "Features", "Influence")
    #A vector with the views
    Views = unlist(lapply(seq_along(modelObj$data), function(i){
        rep(names(modelObj$data)[i], numVars[[i]])
    }))
    names(Views) = colnames(inflMat)
    moltInflMat$View = factor(Views[as.character(moltInflMat$Features)])
    moltInflMat$LatentVariable = factor(moltInflMat$LatentVariable)
    if(grepl(plotType, pattern = "boxplot")){
        moltInflMat$logAbsoluteInfluence  = log(abs(moltInflMat$Influence))
    }

    if(plotType == "lineplot"){
        Plot = ggplot(data = moltInflMat,
                      aes_string(x = "LatentVariable", y = "Influence", group = "Features", colour = "View")) +
            geom_line(size = lineSize) + ylab("Influence per feature")
    } else if(plotType == "boxplot"){
        Plot = ggplot(data = moltInflMat,
                      aes_string(y = "logAbsoluteInfluence", x = "LatentVariable")) +
            geom_boxplot(aes_string(fill = "View"), outlier.shape = NA, outlier.size = 0.75)
        # Means = aggregate(data = moltInflMat, Influence ~ LatentVariable + View, FUN = mean)
        # Plot = Plot + geom_point(data = Means, aes_string(y = "Influence", x = "factor(LatentVariable)"),inherit.aes = FALSE) +
    } else if (plotType == "boxplotSingle"){
        Plot = ggplot(data = moltInflMat,
                      aes_string(y = "logAbsoluteInfluence", x = "View", fill = "View")) + geom_boxplot()
    } else if (plotType == "pointplot"){
    aggMoltInfl = aggregate(data = moltInflMat,
                           Influence ~ View + LatentVariable, FUN = pointFun)
    latentDf = data.frame(LatentVariable = factor(rownames(inflMat)),
                          value = if(constrained) modelObj$alphas[samples,Dim] else
                              modelObj$latentVars[samples,Dim])
    #Scale to fit window
    latentDf$value = latentDf$value*1.1*max(abs(aggMoltInfl$Influence)/max(abs(latentDf$value)))
    Plot = ggplot(data = aggMoltInfl, aes_string(y = "Influence", x = "LatentVariable")) +
        geom_point(aes_string(colour = "View")) +
        ylab(switch(pointFun, "sum" = "Total influence per view",
                    "mean" = "Mean influence per view")) +
        geom_hline(data = data.frame(h = 0), aes(yintercept = h), linetype = "dashed", col ="grey75") +
        geom_point(inherit.aes = FALSE, data = latentDf, aes(x = LatentVariable, y = value),
                   size = 0.5, shape = 3, colour ="grey25")
        #Add reference line and latent variable estimates
    }
    Plot = Plot + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
        ggtitle(paste("Dimension", Dim))
    if(plotType != "boxplotSingle") Plot = Plot +
        xlab(if(constrained) "Environmental variable" else "Latent variable")
    Plot
    }