context("compIntegrate input")

data(hmp2)
test_that("compIntegrate throws warnings when ordination model is fitted", {
 expect_warning(compInt(data = list("microbiome" = microPruneVir), distributions = "quasi",
                            compositional = TRUE, verbose = TRUE, nCores = 1, M = 2))
})
test_that("compIntegrate throws error for wrong input type", {
    expect_error(compInt(data = microPruneVir, distributions = "quasi",
                           compositional = TRUE, verbose = TRUE, nCores = 1, M = 2))
})
#
# test_that("RCM throws errors when only one covariate with one level supplied", {
#   expect_error(RCM(Zeller, covariates = "Age", k = 1))
# })
#
n  = 50;p = 100
tmpMat  = matrix(rnbinom(n*p, mu = 5, size = 1),n,p)
test_that("compIntegrate throws errors when no row names provided", {
    expect_error(
        compInt(list(tmpMat), distributions = "quasi",
                         compositional = TRUE)
        )
})
rownames(tmpMat) = paste("sample", seq_len(n))
tmpMat[1,1] = NA

test_that("compIntegrate throws errors when NAs present in data matrix", {
  expect_error(
      compInt(list(tmpMat), distributions = "quasi",
                       compositional = TRUE)
      )
})

test_that("compIntegrate runs when NAs present in data matrix
          and allowMissingess = TRUE", {
    expect_invisible(
        compInt(list(tmpMat), distributions = "quasi",
                         compositional = TRUE, allowMissingness = TRUE)
        )
})
