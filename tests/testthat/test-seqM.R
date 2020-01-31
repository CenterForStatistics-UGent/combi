context("Index vector of lagrange multipliers of dimension k")

test_that("Number of lagrange multipliers are correct", {
  expect_equal(seqM(1), seq_len(2))
  expect_equal(seqM(2), 3:5)
  expect_equal(seqM(3), 6:9)
  expect_equal(seqM(1, normal = FALSE), seq_len(1))
  expect_equal(seqM(2, normal = FALSE), 2:3)
  expect_equal(seqM(3, normal = FALSE), 4:6)
  expect_equal(seqM(1, nLambda1s = 2), seq_len(3))
  expect_equal(seqM(3, nLambda1s = 3), 10:15)
})
