#' Make triplots or higher dimensional plots of the ordination object
#' @param modelObj the fit
#' @param Dim the dimensions to be plotted
#' @param samDf a dataframe of sample variables
#' @param samCol A variable name from samDf used to colour the samples
#' @param featNum,varNum Integer, the number of features or variables to be plotted
#' @param featCols Colours for the features
#' @param manExpFactorTaxa,manExpFactorVar expansion factors to size all
#' components to the same scale
#' @param featSize,crossSize,varSize,samSize,strokeSize Size of the corresponding objects
#' @param samColValues
#' @param warnMonotonicity A boolean, should a warning be thrown when the
#' feature proportions of compositional views do not all vary monotonically
#' with all latent variables?
#' @param returnCoords A boolean, should coordinates be returned, e.g. for use
#'  in thrird party software
#' @param squarePlot A boolean, should the axes be square? Strongly recommended
#' @param featAlpha Controls the transparacny of the features
#'
#' @return A ggplot object containing the plot
#'
#' @export
#' @import ggplot2
#' @importFrom grDevices terrain.colors
plot.compInt = function(modelObj, Dim = c(1,2), samDf = NULL, samCol = NULL,
                        featNum = 20L, featCols = c("darkblue", "darkgreen",
                                                    "darkred", terrain.colors(5)),
                        manExpFactorTaxa = 0.975, featSize = 2.5, crossSize = 4,
                        manExpFactorVar = 0.975, varNum = nrow(modelObj$alphas),
                        varSize = 2.5, samColValues = NULL, samSize = 1.5,
                        strokeSize = 0.05, warnMonotonicity = FALSE,
                        returnCoords = FALSE, squarePlot = TRUE, featAlpha = 0.5){
    #palette(rainbow(length(unique(samDf[[samCol]]))))
    nViews = length(modelObj$data)
    coords = extractCoords(modelObj, Dim)
    latentData = coords$latentData; featureData = coords$featureData
    varData = coords$varData

    if(!is.null(samDf)){
        latentData = data.frame(latentData, samDf[rownames(modelObj$latentVars),, drop = FALSE])
    }

    DimChar = paste0("Dim", Dim)
    #### Latent variables ####
    Plot = ggplot(aes_string(x = DimChar[1], y = DimChar[2], fill = samCol), data = latentData) +
        geom_point(size = samSize, shape = 21, stroke = strokeSize) + theme_bw() +
        if(!is.null(samColValues)) scale_fill_manual(values = samColValues)

    #### Views ####
    if(length(featNum)==1){featNum = rep(featNum, nViews)}
    checkComp = checkMonotonicity(modelObj, Dim)
    for(i in seq_len(nViews)){
        arrowLengths = rowSums(featureData[[i]][, DimChar]^2)
        featurePlot = arrowLengths >= quantile(arrowLengths, 1-min(1,featNum[i]/nrow(featureData[[i]])))
        featureData[[i]] = featureData[[i]][featurePlot,]
        scalingFactorTmp = apply(latentData[, DimChar], 2, range)/
            apply(featureData[[i]][,DimChar], 2, range)
        scalingFactor = min(scalingFactorTmp[scalingFactorTmp >0]) *
            manExpFactorTaxa
        featureData[[i]][, DimChar] = featureData[[i]][, DimChar]*scalingFactor
        #Warn for non monotonicity in case of compositionality
        if(warnMonotonicity && !all(checkComp[[i]][,featurePlot])){
            checkMate = apply(checkComp[[i]][, featurePlot], c(1,2), all)
warning("Features \n", paste(colnames(modelObj$data[[i]])[featurePlot][!apply(checkMate, 2 ,all)], collapse = "\n"),
        "\nnot monotonous on this plot")
        }
        if(length(featCols[[i]])>1) featCols[[i]] = featCols[[i]][featurePlot]
        Plot = Plot + geom_text(aes_string(x = DimChar[1], y = DimChar[2], label = "featNames"),
                    inherit.aes = FALSE, data = featureData[[i]],
                    col = featCols[[i]], size = featSize, alpha = featAlpha)
    }
    # Views Legend, TO DO!
    #### Gradient ####
    if(!is.null(modelObj$covariates)){

        arrowLengthsVar = rowSums(varData[,DimChar]^2)
        varPlot = arrowLengthsVar >= quantile(arrowLengthsVar, 1-varNum/nrow(modelObj$alphas))
        varData = varData[varPlot,]
        scalingFactorTmpVar = apply(latentData[, DimChar], 2, range)/
            apply(varData[,DimChar], 2, range)
        scalingFactorVar = min(scalingFactorTmpVar[scalingFactorTmpVar >0]) *
            manExpFactorVar
        varData[,DimChar] = varData[,DimChar]*scalingFactorVar
    Plot = Plot + geom_text(aes_string(x = DimChar[1], y = DimChar[2], label = "varNames"),
                            inherit.aes = FALSE, data = varData, size = varSize)
    }
    # Add cross in the centre
    Plot = Plot + geom_point(data = data.frame(x = 0, y = 0),
                             aes_string(x = "x", y = "y"), size = crossSize,
                             inherit.aes = FALSE, shape = 3)
    #Square plot, or throw warning
    if(squarePlot){
        Plot = Plot + coord_fixed()
    } else {
        warning("Plot not square!\n This may be misleading and is not recommended!")
    }
    ## TO DO: add colour code for the views
    if(returnCoords){
        return(list(Plot = Plot, latentData = latentData, featureData = featureData,
               varData = varData))
    } else {
        return(Plot)
    }
}