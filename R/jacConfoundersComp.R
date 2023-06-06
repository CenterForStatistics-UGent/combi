#' Jacobian for conditioning under compositionality
#' @inheritParams scoreConfounders
#' @importFrom tensor tensor
#' @return the jacobian matrix
#' @param marginModel,biasReduction,subtractMax The marginal mode, and booleans
#' indicating bias reduction and maximum subtraction
#' @param confMat,data,meanVarTrend arguments belonging to views
jacConfoundersComp = function(x, confMat, data, meanVarTrend, marginModel,
                              allowMissingness, biasReduction, subtractMax = TRUE){
    parMat = matrix(x, ncol(confMat), ncol(data)-1)
    parMat = cbind(0, parMat)#rowMeans(paramMatrix)
    logCompMat = matrix(marginModel$colOff, nrow(data), ncol(data), byrow = TRUE) +
        + confMat %*% parMat
    if(subtractMax) logCompMat = logCompMat - apply(logCompMat, 1, max)
    CompMat0 = exp(logCompMat)
    Sum = rowSums(CompMat0)
    CompMat = CompMat0/Sum
    libSizes = exp(marginModel$rowOff)
    mu = CompMat*libSizes
    if(allowMissingness){
        isNA = is.na(data)
        data[isNA] = mu[isNA]
    }
    meanEval = meanVarTrend(CompMat, outerProd = FALSE)*libSizes
    S = 1 - CompMat0/Sum
    muMat = mu*(data-mu)/meanEval
    meanEvalDeriv = meanVarTrend(CompMat, deriv = 1L)
    muMatDev = ((data-2*mu) - meanEvalDeriv*mu*(data-mu)/meanEval)/meanEval
    #Other taxon

    Aml = arrayMult(CompMat0, confMat)
    offDiag = tensor(Aml, arrayMult((muMat/Sum*CompMat0 - muMatDev*mu*S)/Sum, confMat, ncol(confMat)), 1,1)
    tmp4 = aperm(offDiag, c(2,1,4,3))
    tmp5  = matrix(c(tmp4), ncol(confMat)*ncol(data), byrow = TRUE)

    onDiag = -tensor(Aml, arrayMult(muMat/(Sum^2)*(Sum-CompMat0) -
                                        muMatDev*libSizes*(Sum-CompMat0)/Sum^2*S, confMat, ncol(confMat)), 1,1)
    onDiagMat = matrix(aperm(onDiag, c(2,1,4,3)), ncol(confMat)*ncol(data), byrow = TRUE)
    blockDiag = as.logical(kronecker(matrix(TRUE, ncol(confMat), ncol(confMat)), X = diag(TRUE, ncol(data))))
    tmp5[blockDiag] = onDiagMat[blockDiag]
    return(tmp5[-seq_len(ncol(confMat)), -seq_len(ncol(confMat))])
}