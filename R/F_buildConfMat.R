#' Build confounder design matrices with and without intercepts
#' @param confounders A dataframe of confounding variables
#' #' For the preliminary trimming, we do not include an intercept,
#' but we do include all the levels of the factors using contrasts=FALSE:
#'  we want to do the trimming in every subgroup, so no hidden reference levels
#'   For the filtering we just use a model with an intercept and
#'    treatment coding, here the interest is only in adjusting the offset
#'
#' @return a list with components
#' \item{confModelMatTrim}{A confounder matrix without intercept, with all
#'  levels of factors present. This will be used to trim out taxa that have
#'   zero abundances in any subgroup defined by confounders}
#' \item{confModelMat}{A confounder matrix with intercept,
#' and with reference levels for factors absent.
#' This will be used to fit the model to modify the independence model,
#' and may include continuous variables}

buildConfMat = function(confounders){
    if(is.null(confounders)){
        return(NULL)
    }
    if (anyNA(confounders)) {
        stop("Confounders contain missing values!\n")
    }
    # No intercept or continuous variables for preliminary trimming
    confModelMatTrim = model.matrix(object = as.formula(paste("~",
            paste(names(confounders)[vapply(FUN.VALUE = TRUE, confounders,
              is.factor)], collapse = "+"), "-1")), data = confounders,
                contrasts.arg = lapply(confounders[vapply(FUN.VALUE = TRUE,
                  confounders, is.factor)], contrasts, contrasts = FALSE))
    # With intercept for filtering
    confModelMat = model.matrix(object = as.formula(paste("~",
              paste(names(confounders), collapse = "+"))), data = confounders,
            contrasts.arg = lapply(confounders[vapply(FUN.VALUE = TRUE,
                      confounders, is.factor)], contrasts, contrasts = TRUE))
    list(confModelMatTrim = confModelMatTrim, confModelMat = confModelMat)
}