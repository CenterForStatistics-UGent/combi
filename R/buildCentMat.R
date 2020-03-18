#' A function to build a centering matrix based on a dataframe
#'
#' @param object an modelDI object or dataframe
#'
#' @return a centering matrix consisting of ones and zeroes,
#'  or a list with components
#' \item{centMat}{a centering matrix consisting of ones and zeroes}
#' \item{datFrame}{The dataframe with factors with one level removed}
buildCentMat = function(object) {
        nFactorLevels = vapply(object, FUN.VALUE = double(1),
                               function(x) {
                                   if (is.factor(x))
                                       nlevels(x) else 1
                               })
        # Number of levels per factor
        oneLevelID = vapply(object, FUN.VALUE = TRUE, function(x) {
            length(unique(x)) == 1
        })
        object[, oneLevelID] = NULL  #Drop factors with one level
        if (any(oneLevelID)) {
            warning("The following variables were not included in the analyses because they have only one value: \n",
                paste(object[oneLevelID], sep = " \n"), immediate. = TRUE)
        }
    # Already prepare the matrix that defines the equations for
    # centering the coefficients of the dummy variables
    centMat = t(vapply(FUN.VALUE = numeric(sum(nFactorLevels)),
        seq_along(nFactorLevels), function(i) {
            c(rep.int(0L, sum(nFactorLevels[seq(0, i - 1)])),
                rep.int(if (nFactorLevels[i] == 1) 0L else 1L,
                nFactorLevels[i]), rep.int(0L, sum(nFactorLevels[-seq(1,
                    i)])))
        }))
    centMat = if (all(rowSums(centMat) == 0)) {
        matrix(0L, 1, sum(nFactorLevels))
    } else {
        centMat[rowSums(centMat) > 0, , drop = FALSE]
    }
    list(centMat = centMat, datFrame = object)
}