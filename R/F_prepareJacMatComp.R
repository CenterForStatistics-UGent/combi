#' prepare the jacobian for the latent variabels compostional
prepareJacMatComp = function(mu, paramEsts, CompMat0, meanVarTrend, data, libSizes){
Sum = rowSums(CompMat0)
CompMat = CompMat0/Sum
CompMat[CompMat==0] = .Machine$double.eps
meanEval = meanVarTrend(CompMat, libSizes, outerProd = FALSE)*libSizes
CC = c(CompMat0 %*% paramEsts)
foo = mu*CC/Sum
# meanEvalSq = meanEval^2
# meanEvalSq[meanEvalSq < .Machine$double.eps] = .Machine$double.eps
# CCsq = CC^2
# CCsq[is.infinite(CCsq)] = 1/.Machine$double.eps
tmp = -((rowMultiply(mu, paramEsts) - foo)^2*
            ((data - mu)*meanVarTrend(CompMat, deriv = 1L)/meanEval + 1) -
            (rowMultiply(mu, paramEsts^2) - 2*rowMultiply(foo, paramEsts) +
                 2*mu*(CC/Sum)^2 - mu*c(CompMat0 %*% paramEsts^2)/Sum)*(data - mu))/meanEval
if(anyNA(tmp)){
    foo = 2
}
tmp
}