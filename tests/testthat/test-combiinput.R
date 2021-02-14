context("combi input")
data(Zhang)
test_that("combi throws warnings when ordination model is fitted", {
 expect_warning(combi(data = list("microbiome" = zhangMicrobio), distributions = "quasi",
                            compositional = TRUE, maxIt = 2L))
})
test_that("combi throws error for wrong input type", {
    expect_error(combi(data = zhangMicrobio, distributions = "quasi",
                           compositional = TRUE, M = 2))
})
n  = 50; p = 100
tmpMat  = matrix(rnbinom(n*p, mu = 5, size = 1),n,p)
test_that("combi throws errors when no row names provided", {
    expect_error(
        combi(list(tmpMat), distributions = "quasi", compositional = TRUE)
        )
})
tmpMat2 = tmpMat
rownames(tmpMat2) = paste("sample", seq_len(n))
tmpMat2[1,1] = NA

test_that("combi throws errors when NAs present in data matrix", {
  expect_error(
      combi(list(tmpMat2), distributions = "quasi", compositional = TRUE)
      )
})

test_that("combi runs when NAs present in data matrix
          and allowMissingess = TRUE", {
    expect_s3_class(
        combi(list(tmpMat2), distributions = "quasi",
                    compositional = TRUE, allowMissingness = TRUE, maxIt = 2L),
        "combi"
        )
      expect_s3_class(
              combi(list(tmpMat2), distributions = "gaussian",
                    compositional = FALSE, allowMissingness = TRUE, maxIt = 2L),
              "combi"
      )
})

test_that("combi throws errors when gaussian distributions are mixed with
          compositional data, or quasi likelihood with non-compositional data", {
            expect_error(
              combi(list(tmpMat2), distributions = "quasi",
                    compositional = FALSE, allowMissingness = TRUE, maxIt = 2L)
            )
            expect_error(
              combi(list(tmpMat2), distributions = "gaussian",
                    compositional = TRUE, allowMissingness = TRUE, maxIt = 2L)
            )
          })

test_that("Polynomial mean-variance model works", {
    expect_s3_class(combi(data = list("microbiome" = zhangMicrobio,
                                      "metabo" = zhangMetabo),
                          distributions = c("quasi", "gaussian"),
                         compositional = c(TRUE, FALSE), maxIt = 2L,
                         meanVarFit = "cubic"), "combi")
})
