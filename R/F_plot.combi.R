#' Make multiplots of the data integration object
#'
#' @param x the fit
#' @param ... additional arguments, currently ignored
#' @param Dim the dimensions to be plotted
#' @param samDf a dataframe of sample variables
#' @param samCol A variable name from samDf used to colour the samples
#' @param samShape A variable name from samDf used to shape the samples
#' @param featCols Colours for the features
#' @param samColValues Colours for the samples
#' @param warnMonotonicity A boolean, should a warning be thrown when the
#' feature proportions of compositional views do not all vary monotonically
#' with all latent variables?
#' @param returnCoords A boolean, should coordinates be returned, e.g. for use
#'  in thrird party software
#' @param squarePlot A boolean, should the axes be square? Strongly recommended
#' @param featAlpha Controls the transparacny of the features
#' @param featNum,varNum The number of features and variables to plot
#' @param manExpFactorTaxa,manExpFactorVar Expansion factors for taxa and
#' variables, normally calculated natively
#' @param featSize,crossSize,varSize,samSize,strokeSize Size parameters for
#' the feature dots, central cross, variable labels, sample dots and sample
#' strokes
#' @param xInd,yInd x and y indentations
#' @param checkOverlap A boolean, should overlapping labels be omitted?
#' @param shapeValues the shapes, as numeric values
#'
#' @return A ggplot object containing the plot
#' @method plot combi
#'
#' @export
#' @import ggplot2
#' @importFrom grDevices terrain.colors
#' @importFrom stats quantile
#' @examples
#' data(Zhang)
#' \dontrun{
#' #Unconstrained
#' microMetaboInt = combi(
#' list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
#' distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
#' logTransformGaussian = FALSE, verbose = TRUE)
#' #Constrained
#' microMetaboIntConstr = combi(
#'     list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
#'     distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
#'     logTransformGaussian = FALSE, covariates = zhangMetavars, verbose = TRUE)}
#'     #Load the fits
#' load(system.file("fits", "zhangFits.RData", package = "combi"))
#' plot(microMetaboInt)
#' plot(microMetaboInt, samDf = zhangMetavars, samCol = "ABX")
#' #Constrained
#' plot(microMetaboIntConstr, samDf = zhangMetavars, samCol = "ABX")
plot.combi = function(x, ..., Dim = c(1,2), samDf = NULL, samCol = NULL,
                        samShape = NULL, featNum = 15L,
                        featCols = c("darkblue", "darkgreen", "darkred",
                                     terrain.colors(5)),
                        manExpFactorTaxa = 0.975, featSize = 2.5, crossSize = 4,
                        manExpFactorVar = 0.975, varNum = nrow(x$alphas),
                        varSize = 2.5, samColValues = NULL, samSize = 1.5,
                        strokeSize = 0.05, warnMonotonicity = FALSE,
                        returnCoords = FALSE, squarePlot = TRUE, featAlpha = 0.5,
                        xInd = 0, yInd = 0, checkOverlap = FALSE,
                        shapeValues = (21:(21+length(unique(samDf[[samShape]]))))){
    #palette(rainbow(length(unique(samDf[[samCol]]))))
    nViews = length(x$data)
    coords = extractCoords(x, Dim)
    latentData = coords$latentData; featureData = coords$featureData
    varData = coords$varData

    if(!is.null(samDf)){
        latentData = data.frame(latentData, samDf[rownames(x$latentVars),, drop = FALSE])
    }

    DimChar = paste0("Dim", Dim)
    if(!is.null(samShape)){
        if(!(is.factor(samDf[[samShape]]) | is.character(samDf[[samShape]])))
            stop("Shape must be a discrete variable!\n")
    }
    #### Latent variables ####
    Plot = ggplot(aes_string(x = DimChar[1], y = DimChar[2], fill = samCol, shape = samShape), data = latentData) +
        theme_bw() +
        (if(!is.null(samColValues)) scale_fill_manual(values = samColValues, name = samCol))
    if(!is.null(samShape)) {
        Plot = Plot + scale_shape_manual(name = samShape, values = shapeValues) +
             geom_point(size = samSize, stroke = strokeSize) +
            guides(fill = guide_legend(override.aes = list(shape = 21)))}
    else {
        Plot = Plot+ geom_point(size = samSize, shape = 21, stroke = strokeSize)
}

    #### Views ####
    if(length(featNum)==1){featNum = rep(featNum, nViews)}
    checkComp = checkMonotonicity(x, Dim)
    for(i in seq_len(nViews)){
        if(featNum[i]==0) break
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
warning("Features \n", paste(colnames(x$data[[i]])[featurePlot][!apply(checkMate, 2 ,all)], collapse = "\n"),
        "\nnot monotonous on this plot")
        }
        if(length(featCols[[i]])>1) featCols[[i]] = featCols[[i]][featurePlot]
        Plot = Plot + geom_text(aes_string(x = DimChar[1], y = DimChar[2], label = "featNames"),
                    inherit.aes = FALSE, data = featureData[[i]],
                    col = featCols[[i]], size = featSize, alpha = featAlpha,
                    check_overlap = checkOverlap)
    }
    # Views Legend, TO DO!
    #### Gradient ####
    if(!is.null(x$covariates)){

        arrowLengthsVar = rowSums(varData[,DimChar]^2)
        varPlot = arrowLengthsVar >= quantile(arrowLengthsVar, 1-varNum/nrow(x$alphas))
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
        warning("Plot not square!\n This may be misleading and is not recommended!
                See squarePlot parameter")
    }
    ## TO DO: add colour code for the views
    Plot = indentPlot(Plot, xInd = xInd, yInd = yInd) #Indentation
    if(returnCoords){
        return(list(Plot = Plot, latentData = latentData, featureData = featureData,
               varData = varData))
    } else {
        return(Plot)
    }
}