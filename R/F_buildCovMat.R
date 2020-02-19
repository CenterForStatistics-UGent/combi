#' A function to build the covariate matrix of the constraints
#'
#' @param datFrame the dataframe with which the covariate matrix is to be built
#'
#' In this case we will 1) Include dummy's for every level of the
#'  categorical variable, and force them to sum to zero.
#'  This is needed for plotting and required for
#'   reference level indepenent normalization.
#'   2) Exclude an intercept. The density function f()
#'   will provide this already.
#'
#' @return a list with components
#' \item{covModelMat}{The model matrix}
#' \item{datFrame}{The dataframe used to construct the model matrix}
#'
#' @import stats
buildCovMat = function(datFrame) {

    logVec = vapply(FUN.VALUE = TRUE, datFrame, is.logical)
    intVec = vapply(FUN.VALUE = TRUE, datFrame, is.integer)
    charVec = vapply(FUN.VALUE = TRUE, datFrame, is.character)

    if (any(logVec)) {
        datFrame[, logVec] = lapply(datFrame[logVec], as.factor)
        # Convert logicals to factors warning('Logicals converted to
        # factors! \n', immediate. = TRUE). No warning needed
    }
    if (any(intVec)) {
        datFrame[, intVec] = lapply(datFrame[intVec], as.numeric)
        # Convert integers to numeric
        warning("Integer values treated as numeric! \n", immediate. = TRUE)
    }
    if (any(charVec)) {
        datFrame[, charVec] = lapply(datFrame[charVec], factor)
        # Convert characters to factor
        warning("Character vectors treated as factors! \n", immediate. = TRUE)
    }
    nFactorLevels = vapply(FUN.VALUE = integer(1), datFrame,
                           function(x) {
                               if (is.factor(x))
                                   nlevels(x) else 0L
                           })  #Number of levels per factor
    covariatesNames = names(datFrame)

    covariatesNames = covariatesNames[!(vapply(FUN.VALUE = TRUE,
                                               datFrame, is.factor) & (nFactorLevels == 1L))]
    # Drop factors with one level
    if (any(nFactorLevels==1L)){
        warning(immediate. = TRUE, "Factors with only one level dropped!")
    }
    nFactorLevels = nFactorLevels[covariatesNames]
    datFrame = datFrame[, covariatesNames, drop = FALSE]
    #Check for alias structures
    checkAlias(datFrame, covariatesNames)

    if (any(vapply(FUN.VALUE = TRUE, datFrame, is.factor) & (nFactorLevels <
                                                             2))) {
        warning("The following variables were not included in the analyses
            because they are factors with only one level: \n",
                paste(covariatesNames[vapply(FUN.VALUE = TRUE, datFrame,
                                        is.factor) & (nFactorLevels < 2)], sep = " \n"),
                immediate. = TRUE, call. = FALSE)
    }
    # Center and scale the continuous covariates
    datFrame[vapply(FUN.VALUE = TRUE, datFrame, is.numeric)] =
        scale(datFrame[vapply(FUN.VALUE = TRUE,datFrame, is.numeric)])

    covModelMat = model.matrix(
        object = formula(paste("~", paste(covariatesNames,
                                          collapse = "+"), "-1")), data = datFrame,
        contrasts.arg = lapply(datFrame[vapply(datFrame,
                                               is.factor, FUN.VALUE = TRUE)], contrasts, contrasts = FALSE))
    if (NCOL(covModelMat) == 1)
        stop("A constrained ordination with only one variable is meaningless.
             Please provide more covariates or perform an unconstrained analysis.",
             call. = FALSE)
    if((NCOL(covModelMat) - sum((vapply(FUN.VALUE = TRUE,
                                datFrame, is.factor))))> NROW(covModelMat))
        stop("More parameters than samples in constrained analysis!\n ")
    list(covModelMat = covModelMat, datFrame = datFrame)
    }