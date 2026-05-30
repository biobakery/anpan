test_that("pglmm runs", {
  meta = data.frame(x = rnorm(100), sample_id = paste0("t", 1:100))
  tr = ape::rtree(100)

  # There will be warnings because the chains are so short.

  # This is mainly to force the models to compile during tests. This is done so
  # root can compile the models on the Guacamole VMs.
  opt = options("cmdstanr_warn_inits")[[1]]

  options(cmdstanr_warn_inits = FALSE)
  suppressMessages(suppressWarnings({res = anpan_pglmm(meta, tr,
                    outcome = "x",
                    iter_sampling = 10,
                    iter_warmup = 10,
                    show_plot_cor_mat = FALSE,
                    show_plot_tree = FALSE,
                    show_post = FALSE,
                    show_messages = FALSE,
                    show_exceptions = FALSE,
                    verbose = FALSE,
                    run_diagnostics = FALSE,
                    )}))

  options("cmdstanr_warn_inits" = opt)
  expect_type(res, "list")
})

test_that("logistic pglmm runs", {
  meta = data.frame(x = as.logical(rbinom(100, 1, prob = .5)), sample_id = paste0("t", 1:100))
  tr = ape::rtree(100)

  opt = options("cmdstanr_warn_inits")[[1]]
  options(cmdstanr_warn_inits = FALSE)

  # There will be warnings because the chains are so short
  suppressMessages(suppressWarnings({res = anpan_pglmm(meta, tr,
                    family = "binomial",
                    outcome = "x",
                    iter_sampling = 10,
                    iter_warmup = 10,
                    show_plot_cor_mat = FALSE,
                    show_plot_tree = FALSE,
                    show_post = FALSE,
                    show_messages = FALSE,
                    show_exceptions = FALSE,
                    verbose = FALSE,
                    run_diagnostics = FALSE,
                    )}))

  options("cmdstanr_warn_inits" = opt)
  expect_type(res, "list")
})
