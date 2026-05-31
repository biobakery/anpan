
test_that("arma::solve matches base", {
  X = matrix(rnorm(100), ncol = 10)
  A = t(X) %*% X

  expect_equal(arma_solve(A),
               solve(A))
})
