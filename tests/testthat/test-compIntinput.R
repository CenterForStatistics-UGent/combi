context("compIntegrate input")

data(Zhang)
test_that("compIntegrate throws warnings when ordination model is fitted", {
 expect_warning(compInt(data = list("microbiome" = zhangMicrobio), distributions = "quasi",
                            compositional = TRUE, nCores = 1, M = 2))
})
test_that("compIntegrate throws error for wrong input type", {
    expect_error(compInt(data = zhangMicrobio, distributions = "quasi",
                           compositional = TRUE, nCores = 1, M = 2))
})
n  = 50; p = 100
tmpMat  = matrix(rnbinom(n*p, mu = 5, size = 1),n,p)
test_that("compIntegrate throws errors when no row names provided", {
    expect_error(
        compInt(list(tmpMat), distributions = "quasi",
                         compositional = TRUE)
        )
})
tmpMat2 = tmpMat
rownames(tmpMat2) = paste("sample", seq_len(n))
tmpMat2[1,1] = NA

test_that("compIntegrate throws errors when NAs present in data matrix", {
  expect_error(
      compInt(list(tmpMat2), distributions = "quasi",
                       compositional = TRUE)
      )
})

test_that("compIntegrate runs when NAs present in data matrix
          and allowMissingess = TRUE", {
    expect_s3_class(
        compInt(list(tmpMat2), distributions = "quasi",
                         compositional = TRUE, allowMissingness = TRUE),
        "compInt"
        )
})
