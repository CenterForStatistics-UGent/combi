#' Filter out the effect of known confounders
#'
#' @param data data matrix
#' @param distribution,link,invLink,compositional,meanVarTrend,offSet,numVar,marginModel
#' Characteristics of the view
#' @param confMat A confounder design matrix
#' @param control A list of control elements to the nleqslv function
#' @param biasReduction A boolean, should bias reduction be applied
#' @param allowMissingness A boolean, are missing values allowed?
#'
#' @return Parameter estimates accounting for the effects of the confounders
filterConfounders = function(confMat, data, distribution, link, invLink,
                             compositional, control, meanVarTrend, offSet,
                             numVar, marginModel, biasReduction, allowMissingness){
    if(compositional){
        #Remove parameters for one taxon, set them to zero (additive log-ratio transform)
        #Choice of reference does not mattter, the point is to modify the offset
        paramEstsTmp = nleqslv(x = numeric(ncol(confMat)*(numVar-1)),
                               fn = scoreConfoundersComp,
                               jac = jacConfoundersComp, confMat = confMat,
                               data = data, control = control,
                               meanVarTrend = meanVarTrend,
                            marginModel = marginModel, biasReduction = biasReduction,
                            allowMissingness = allowMissingness)
        if(paramEstsTmp$termcd!=1){
            paramEstsTmp2 = BBsolve(par = paramEstsTmp$x, fn = scoreConfoundersComp,
                                   confMat = confMat, data = data,
                                   meanVarTrend = meanVarTrend,
                                   marginModel = marginModel, biasReduction = biasReduction,
                                   allowMissingness = allowMissingness)
            if(paramEstsTmp2$convergence!=0){
                warning("Estimation of conditioning parameters did not converge!
                    Try stricter filtering using the prevCutOff argument.")
            }
        }
        paramEsts = cbind(0, matrix(if(paramEstsTmp$termcd==1) paramEstsTmp$x else paramEstsTmp2$par,
                                    nrow = ncol(confMat)))

    } else {
        paramEsts = vapply(FUN.VALUE = numeric(ncol(confMat)), seq_len(numVar),
                       function(i){
            offSet = with(marginModel, exp(outer(rowOff, colOff, "+")))
            libSizes = exp(marginModel$rowOff)
            comp = exp(marginModel$colOff[i]) #Use baseline composition, not ideal
            # but at least it may fit
            tmp = try(nleqslv(x = numeric(ncol(confMat)),
            fn = scoreConfounders, jac = jacConfounders, offSet = offSet[,i],
            distribution = distribution, meanVarTrend = meanVarTrend,
            data = data[,i], CompMat = comp,
            libSizes = libSizes, confMat = confMat, control = control,
            allowMissingness = allowMissingness)$x, silent = TRUE)
            #There may be numerical inaccuracy, try with numerical jacobian
            if(inherits(tmp, "try-error")){
            tmp = nleqslv(x = numeric(ncol(confMat)),
                          fn = scoreConfounders, jac = NULL, offSet = offSet[,i],
                          distribution = distribution, meanVarTrend = meanVarTrend,
                          data = data[,i], confMat = confMat,
                          CompMat = comp, allowMissingness = allowMissingness,
                          libSizes = libSizes, control = control)$x
            }
            return(tmp)})
    }
   paramEsts
}