woodbury_s22_inv_old = function(cor_mat_inv, cor_mat, j) {
  # Redone with Rcpp, see binom_ll.cpp

  # https://en.wikipedia.org/wiki/Woodbury_matrix_identity
  # >20x faster, but less numerically precise by a tiny amount.

  n = nrow(cor_mat)
  ord = c(j, (1:n)[-j])
  rcor_mat = cor_mat[ord,ord]
  rcor_mat_inv = cor_mat_inv[ord,ord]

  U = cbind(c(1, rep(0, n-1)),
            c(0, rcor_mat[-1, 1]))

  V = t(U)[2:1,]

  update_factor = solve(diag(2) - V %*% (rcor_mat_inv %*% U))

  woodbury_ans = rcor_mat_inv + (rcor_mat_inv %*% U) %*% update_factor %*% (V %*% rcor_mat_inv)

  woodbury_ans[-1,-1]
}

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
