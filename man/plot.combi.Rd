% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F_plot.combi.R
\name{plot.combi}
\alias{plot.combi}
\title{Make multiplots of the data integration object}
\usage{
\method{plot}{combi}(x, ..., Dim = c(1, 2), samDf = NULL,
  samCol = NULL, samShape = NULL, featNum = 20L,
  featCols = c("darkblue", "darkgreen", "darkred", terrain.colors(5)),
  manExpFactorTaxa = 0.975, featSize = 2.5, crossSize = 4,
  manExpFactorVar = 0.975, varNum = nrow(x$alphas), varSize = 2.5,
  samColValues = NULL, samSize = 1.5, strokeSize = 0.05,
  warnMonotonicity = FALSE, returnCoords = FALSE, squarePlot = TRUE,
  featAlpha = 0.5, xInd = 0, yInd = 0, checkOverlap = FALSE,
  shapeValues = (21:(21 + length(unique(samDf[[samShape]])))))
}
\arguments{
\item{x}{the fit}

\item{...}{additional arguments, currently ignored}

\item{Dim}{the dimensions to be plotted}

\item{samDf}{a dataframe of sample variables}

\item{samCol}{A variable name from samDf used to colour the samples}

\item{samShape}{A variable name from samDf used to shape the samples}

\item{featNum, varNum}{The number of features and variables to plot}

\item{featCols}{Colours for the features}

\item{manExpFactorTaxa, manExpFactorVar}{Expansion factors for taxa and
variables, normally calculated natively}

\item{featSize, crossSize, varSize, samSize, strokeSize}{Size parameters for
the feature dots, central cross, variable labels, sample dots and sample
strokes}

\item{samColValues}{Colours for the samples}

\item{warnMonotonicity}{A boolean, should a warning be thrown when the
feature proportions of compositional views do not all vary monotonically
with all latent variables?}

\item{returnCoords}{A boolean, should coordinates be returned, e.g. for use
in thrird party software}

\item{squarePlot}{A boolean, should the axes be square? Strongly recommended}

\item{featAlpha}{Controls the transparacny of the features}

\item{xInd, yInd}{x and y indentations}

\item{checkOverlap}{A boolean, should overlapping labels be omitted?}

\item{shapeValues}{the shapes, as numeric values}
}
\value{
A ggplot object containing the plot
}
\description{
Make multiplots of the data integration object
}
\examples{
data(Zhang)
\dontrun{
#Unconstrained
microMetaboInt = combi(
list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
logTransformGaussian = FALSE, verbose = TRUE)
#Constrained
microMetaboIntConstr = combi(
    list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
    distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
    logTransformGaussian = FALSE, covariates = zhangMetavars, verbose = TRUE)}
load(system.file("fits", "zhangFits.RData", package = "combi"))
plot(microMetaboInt)
plot(microMetaboInt, samDf = zhangMetavars, samCol = "ABX")
#Constrained
plot(microMetaboIntConstr, samDf = zhangMetavars, samCol = "ABX")
}