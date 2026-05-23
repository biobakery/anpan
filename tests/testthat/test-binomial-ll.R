
test_that("Binomial loglik calculation accurate", {

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

  expect_true(all(dplyr::near(test_ll - binom_ll,
                                  0)))

})
