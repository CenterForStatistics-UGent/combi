#' Plot the convegrence of the different parameter estimates in a line plot
#' @param model A fitted modelDI object
#' @param latent A boolean, should latent variable trajectory be plotted
#' @param nVars An integer, the number of variables to plot. By default all are plotted
#' @param Dim An integer, the dimension to be plotted
#' @param View An integer or character string, indicating the view
#' to be plotted (if latent = FALSE)
#' @param size The line size (see ?geom_path)
#'
#' @importFrom reshape2 melt
#'
#' @return A ggplot object containing the convergence plot
#' @export
#' @examples
#' \donttest{
#' data(Zhang)
#' #Unconstrained
#' microMetaboInt = combi(
#' list("microbiome" = zhangMicrobio, "metabolomics" = zhangMetabo),
#' distributions = c("quasi", "gaussian"), compositional = c(TRUE, FALSE),
#' logTransformGaussian = FALSE, verbose = TRUE)}
#' load(system.file("extdata", "zhangFits.RData", package = "combi"))
#' convPlot(microMetaboInt)
#' convPlot(microMetaboInt, Dim = 2)
#' convPlot(microMetaboInt, View = "microbiome")
convPlot = function(model, latent = is.null(View), nVars = Inf, Dim = 1L,
                    View = NULL, size = 0.125){
    if(is.null(model$paramRec)){
        stop("Intermediate parameter estimates not recorded.
             Try rerunning with record = TRUE.\n")
    }
    stopifnot(is(model, "combi"), is.logical(latent), is.numeric(Dim),
              length(Dim)==1, is.numeric(nVars), is.numeric(View) |
                  is.character(View) | is.null(View), length(View) %in% c(0,1),
              is.numeric(size))
    df = if(latent){
if(is.null(model$covariates)) model$latentRec[,Dim,seq_len(model$iter[[Dim]])] else
    model$alphaRec[,Dim,seq_len(model$iter[[Dim]])]
    } else {
model$paramRec[[View]][Dim,,seq_len(model$iter[[Dim]])]
    }
    nVars = min(nVars, nrow(df))
    #Only plot variabels with largest relative departure
    whichVars = order(apply(df[,(ncol(df)-1):ncol(df)], 1, function(x){
abs(1-x[2]/x[1])
    })) <= nVars
    moltDf = melt(df[whichVars,])
    names(moltDf) = c("Variable","Iteration", "Value")
    moltDf$Iteration = moltDf$Iteration - 1 #Shift back one to make the lines fit
    #Add colour based on convergence
    Gaussian = !is.null(View) && model$distributions[[View]]=="gaussian"
    moltDf$Colour = if(Gaussian){
        1
    } else {
        rep(if(latent) model$latentConv[seq_len(model$iter[[Dim]]),Dim] else
        model$paramConv[seq_len(model$iter[[Dim]]),View,Dim],
        each = sum(whichVars))}
    moltDf$Colour = factor(labels = c("Yes", "Stalled", "No",
                                      "Jacobian too\n ill-conditioned"),
                           vapply(moltDf$Colour, FUN.VALUE = character(1),
                                  switch, "1" = "Yes", "2" = "Stalled",
                                  "3" ="Stalled",
                                  "5" = "Jacobian too\n ill-conditioned", "No"),
                           levels = c("Yes", "Stalled",
                                      "Jacobian too\n ill-conditioned", "No"))
    gTitle = paste("Trace plots of", if(latent) {
        if(is.null(model$covariates))"latent variable estimates" else
            "gradient estimates"
        } else {
        paste("feature parameter estimates of view", View)},
        "of dimension", Dim)
    ggplot(moltDf, aes(x = Iteration, y = Value, group = Variable, col = Colour)) +
        geom_path(size = size) +
        theme_bw() +
        ylab("Parameter estimate") +
        ggtitle(gTitle) +
        scale_colour_manual(values = c("Yes" = "black", "Stalled" = "blue",
                                       "No" = "red", "Jacobian too\n ill-conditioned" = "darkgreen"),
                            name = paste0(if(latent || model$newtonRaphson[View]) "Newton-Raphson" else "Barilai-Borwein",
                                          "\n convergence"), guide = if(Gaussian) FALSE else guide_legend())
}
