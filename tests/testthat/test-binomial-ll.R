
test_that("Binomial loglik close to original implementation", {

  load(test_path("test_data/binom_pieces.RData"))
  load(test_path("test_data/binom_ll.RData"))

  test_ll = get_ll_mat(pieces$nested_df,
                       pieces$em,
                       pieces$cor_mat,
                       pieces$Lcov,
                       pieces$Xc,
                       pieces$offset_val,
                       pieces$metadata$outcome,
                       pieces$family, verbose = FALSE)

  diffs = test_ll - binom_ll

  max(abs(diffs))

  # observations 3 and 9 seem problematic
  which(diffs == max(abs(diffs)), arr.ind = TRUE)

  comps = dplyr::near(diffs, 0, tol = .002)
  # Have to increase the tolerance a bit because turns out the newer
  # implementation is more accurate for low sigma_phylo values.

  expect_true(all(comps))

})

test_that("Binomial loglik accurate", {

  load(test_path("test_data/binom_pieces2.RData"))
  load(test_path("test_data/binom_ll2.RData"))

  test_ll = get_ll_mat(pieces$nested_df,
                       pieces$em,
                       pieces$cor_mat,
                       pieces$Lcov,
                       pieces$Xc,
                       pieces$offset_val,
                       pieces$metadata$outcome,
                       pieces$family, verbose = FALSE)

  diffs = test_ll - binom_ll

  max(abs(diffs))

  # observations 3 and 9 seem problematic
  which(diffs == max(abs(diffs)), arr.ind = TRUE)

  comps = dplyr::near(diffs, 0)
  # Have to increase the tolerance a bit because turns out the newer
  # implementation is more accurate for low sigma_phylo values.

  expect_true(all(comps))

})

test_that("C++ inv_logit works", {
  x = rnorm(100)

  expect_equal(inv_logit(x),
               plogis(x))
})

test_that("C++ log integrand accurate", {
  tv = c(0.0256493935170817, 0.440822642358155, 0.0624806907205936,
         1, 1.05692176403115, 0, 2.506628274631, 1)
  effect_mean_j = tv[1]
  mu_bar_j = tv[2]
  sigma_bar_j = tv[3]
  yj = tv[4]
  lm_mean = tv[5]
  sqrt2pi = tv[7]


  expect_equal(integrand_logistic(effect_mean_j, mu_bar_j, sigma_bar_j, yj, lm_mean, 0, sqrt2pi , 1),
               -20.5146062669803)
})

test_that("C++ log integrand matches old version", {
  tv = c(0.0256493935170817, 0.440822642358155, 0.0624806907205936,
         1, 1.05692176403115, 0, 2.506628274631, 1)
  effect_mean_j = tv[1]
  mu_bar_j = tv[2]
  sigma_bar_j = tv[3]
  yj = tv[4]
  lm_mean = tv[5]
  sqrt2pi = tv[7]


  expect_equal(integrand_logistic(effect_mean_j, mu_bar_j, sigma_bar_j, yj, lm_mean, 0, sqrt2pi , 1),
               vec_integrand_logistic_old(effect_mean_j, mu_bar_j, sigma_bar_j, yj, lm_mean, 0, sqrt2pi , 1))
})

# test_that("C++ observation-wise integration at iteration i correct", {
#   load(test_path("test_data/lli_data.RData"))
#
#   res = llij_binom(lm_means,
#                    sigma12x22_inv_mat,
#                    cor21_mat,
#                    sigma_phylo,
#                    i_df$phylo_effects[[1]],
#                    Y,
#                    effect_means,
#                    sqrt(2*pi))
# })
