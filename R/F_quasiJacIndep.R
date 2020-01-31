#' The jacobian for column offset estimation
#' @inheritParams quasiScoreIndep
#' @return the jacobian matrix
quasiJacIndep = function(x, data, otherMargin, meanVarTrend, col, libSizes,...){
    mu = exp(buildMuMargins(x = x, otherMargin = otherMargin, col = col))
    if(anyNA(mu)){
        return(rep(1e16, length(x)))
    }
    tmp = diag((if(col) colSums else rowSums)(prepareJacMat(data = data, mu = mu,
                               meanVarTrend = meanVarTrend,
                               CompMat = exp(if(col) x else otherMargin),
                               libSizes = libSizes, ...), na.rm  = TRUE))
    tmp
  }