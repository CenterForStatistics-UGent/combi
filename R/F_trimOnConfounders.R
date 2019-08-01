#' Trim based on confounders to avoid taxa with only zero counts
#'
#' @param confounders a nxt confounder matrix
#' @param data the data matrix
#' @param prevCutOff a scalar between 0 and 1, the prevalence cut off
#' @param minFraction a scalar between 0 and 1,
#'  each taxon's total abundance should equal at least the number of samples n
#'   times minFraction, otherwise it is trimmed
#' @param n the number of samples
#'
#' Should be called prior to fitting the independence model
#'
#' @return A trimmed data matrix nxp'
trimOnConfounders = function(confounders,
                             data, prevCutOff, minFraction, n){
    # Over taxa
    trimmingID = apply(data, 2, function(x) {
        # Over confounding variables
        any(apply(confounders, 2, function(conf) {
            tapply(X = x, INDEX = conf, FUN = function(y) {
                mean(y == 0, na.rm = TRUE) >= prevCutOff |
                    sum(y, na.rm = TRUE) < (n * minFraction)
            })  #Any all-zero subgroup?
        }))
    })
    if (sum(!trimmingID) <= 1) {
        stop("All taxa would be trimmed,
        please provide a covariate with less levels,
        or reduce the prevalence cut-off! \n")
    }
    data[, !trimmingID]  #Return trimmed data
}