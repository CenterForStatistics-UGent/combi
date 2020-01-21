context("combi input")

data(Zhang)
test_that("combi throws warnings when ordination model is fitted", {
 expect_warning(combi(data = list("microbiome" = zhangMicrobio), distributions = "quasi",
                            compositional = TRUE, nCores = 1, maxIt = 3))
})
test_that("combi throws error for wrong input type", {
    expect_error(combi(data = zhangMicrobio, distributions = "quasi",
                           compositional = TRUE, nCores = 1, M = 2))
})
n  = 50; p = 100
tmpMat  = matrix(rnbinom(n*p, mu = 5, size = 1),n,p)
test_that("combi throws errors when no row names provided", {
    expect_error(
        combi(list(tmpMat), distributions = "quasi",
                         compositional = TRUE)
        )
})
tmpMat2 = tmpMat
rownames(tmpMat2) = paste("sample", seq_len(n))
tmpMat2[1,1] = NA

test_that("combi throws errors when NAs present in data matrix", {
  expect_error(
      combi(list(tmpMat2), distributions = "quasi",
                       compositional = TRUE)
      )
})

test_that("combi runs when NAs present in data matrix
          and allowMissingess = TRUE", {
    expect_s3_class(
        combi(list(tmpMat2), distributions = "quasi",
                         compositional = TRUE, allowMissingness = TRUE),
        "combi"
        )
})
