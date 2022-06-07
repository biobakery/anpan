test_that("pglmm runs", {
  meta = data.frame(x = rnorm(100), sample_id = paste0("t", 1:100))
  tr = ape::rtree(100)
  res = anpan_pglmm(meta, tr,
                    outcome = "x",
                    iter_sampling = 10,
                    iter_warmup = 10,
                    show_plot_cor_mat = FALSE,
                    show_plot_tree = FALSE)

  expect_is(res, "list")
})
