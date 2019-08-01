#' Filter out the effect of known confounders
#' @param indepModel
#' @param data
#' @param distribution
#' @param link
#' @param invLink
#' @param compositional
filterConfounders = function(confMat, data, distribution, link, invLink,
                             compositional, control, meanVarTrend, offSet,
                             numVar, marginModel, biasReduction, allowMissingness,
                             maxItFilt){
    if(compositional){
        #Remove parameters for one taxon, set them to zero (additive log-ratio transform)
        #Choice of reference does not mattter, the point is to modify the offset
        paramEstsTmp = nleqslv(x = numeric(ncol(confMat)*(numVar-1)), fn = scoreConfoundersComp, jac =
                                jacConfoundersComp, confMat = confMat, data = data,
                            control = control, meanVarTrend = meanVarTrend,
                            marginModel = marginModel, biasReduction = biasReduction,
                            allowMissingness = allowMissingness)
        if(paramEstsTmp$termcd!=1){
            # iter = 0
            # while(paramEstsTmp$termcd!=1 && ((iter = iter+1) < maxItFilt)){
            #     paramEstsTmp = nleqslv(x = rnorm(ncol(confMat)*(numVar-1), sd = 1.1^iter), fn = scoreConfoundersComp,
            #                            jac = jacConfoundersComp, confMat = confMat, data = data,
            #                            control = control, meanVarTrend = meanVarTrend,
            #                            marginModel = marginModel, biasReduction = biasReduction,
            #                            allowMissingness = allowMissingness)
            # }
            paramEstsTmp2 = BBsolve(par = paramEstsTmp$x, fn = scoreConfoundersComp,
                                   confMat = confMat, data = data,
                                   meanVarTrend = meanVarTrend,
                                   marginModel = marginModel, biasReduction = biasReduction,
                                   allowMissingness = allowMissingness)
            if(paramEstsTmp2$convergence!=0){
                warning("Estimation of conditioning parameters did not converge!
                    Try stricter filtering using the prevCutOff argument")
            }
        }
        paramEsts = cbind(0, matrix(if(paramEstsTmp$termcd==1) paramEstsTmp$x else paramEstsTmp2$par,
                                    nrow = ncol(confMat)))

    } else {
        paramEsts = vapply(FUN.VALUE = numeric(ncol(confMat)), seq_len(numVar),
                       function(i){
            tmp = try(nleqslv(x = numeric(ncol(confMat)),
            fn = scoreConfounders, jac = jacConfounders, offSet = offSet[,i],
            distribution = distribution, meanVarTrend = meanVarTrend,
            data = data[,i], confMat = confMat, control = control)$x, silent = TRUE)
            #There may be numerical inaccuracy, try with numerical jacobian
            if(inherits(tmp, "try-error")){
            tmp = nleqslv(x = paramEsts[,i],
                              fn = scoreConfounders, jac = NULL, offSet = offSet[,i],
                              distribution = distribution, meanVarTrend = meanVarTrend,
                              data = data[,i], confMat = confMat, control = control)$x
            }
            return(tmp)})
    }
   paramEsts
}