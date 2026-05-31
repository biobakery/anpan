
test_that("arma::solve matches base", {
  X = matrix(rnorm(100), ncol = 10)
  A = t(X) %*% X

  expect_equal(arma_solve(A),
               solve(A))
})

test_that("arma woodbury implementation matches old version", {
  n = 30
  tr = ape::rtree(n)
  cor_mat = ape::vcv.phylo(tr, corr = TRUE)

  cor_mat_inv = arma_solve(cor_mat)

  j = 4
  o2 = (1:n)[-j] - 1

  arma_res = woodbury_s22_inv(cor_mat_inv, cor_mat, j-1, o2, c(j-1, o2))

  expect_equal(arma_res,
               unname(woodbury_s22_inv_old(cor_mat_inv, cor_mat, j)))
})
