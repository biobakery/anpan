test_that("multiplication works", {
  bug_pwy_dat = data.frame(pwy = factor(rep(paste0("pwy", 1:5),
                                     each = 100),
                                  levels = paste0("pwy", 1:5)),
                           log10_species_abd = rnorm(500),
                           log10_pwy_abd = rnorm(500),
                           group01 = sample(c(0,1),
                                            size = 500,
                                            replace = TRUE))
  res = anpan_pwy_ranef(bug_pwy_dat,
                        group_ind = "group01",
                        iter_sampling = 10,
                        iter_warmup = 10)

  expect_s3_class(res, "data.frame")
})
