## code to prepare `DATASET` dataset goes here

library(anpan)
library(tibble)
library(dplyr)
devtools::load_all("~/projects/anpan")

set.seed(123)

n = 20
sigma_phylo = 1
sigma_resid = 1

tr = ape::rtree(n)

cor_mat = ape::vcv.phylo(tr, corr = TRUE)

set.seed(42)

covariate = rnorm(n)

linear_term = 1 * covariate + rnorm(n, mean = 0, sd = sigma_resid)

true_phylo_effects = sigma_phylo * draw_mvnorm(Sigma = cor_mat)

metadata = tibble(sample_id = colnames(cor_mat),
                  covariate = covariate,
                  outcome   = linear_term + true_phylo_effects)

m2 = metadata |>
  mutate(outcome = rbinom(n, 1, plogis(outcome)))

result = anpan_pglmm(meta_file       = metadata,
                     tree_file       = tr,
                     outcome         = "outcome",
                     covariates      = "covariate",
                     family          = "gaussian",
                     bug_name        = "sim_bug",
                     reg_noise       = TRUE,
                     loo_comparison  = TRUE,
                     run_diagnostics = FALSE,
                     refresh         = 500,
                     show_plot_tree  = FALSE,
                     show_post       = FALSE,
                     iter_sampling  = 100)

gauss_ll = result$loo$pglmm_ll_df |> as.matrix()
gauss_ddf = result$pglmm_fit$draws(format = "data.frame")
gauss_loo = get_pglmm_loo(gauss_ll, gauss_ddf)

# Other arguments to get_ll_mat() ----
em = result$pglmm_fit$summary(variables = "phylo_effect", mean = mean)$mean

cor_mat = result$cor_mat

chol_res = safely_chol(cor_mat)
Lcov = t(chol_res$result)

Xc = matrix(nrow = nrow(metadata),
            ncol = 1)

mx = mean(metadata$covariate)

if (ncol(Xc) > 0) {
  # can't remember how to use seq_along here
  for (i in 2) {
    Xc[,i-1] = metadata$covariate - mx[i-1]
  }
}

offset_val = rep(0, n)

family = "gaussian"

draw_dt = result$pglmm_fit$draws(format = "data.frame") |>
  tibble::as_tibble() |>
  select(-tidyselect::matches("std_phylo|yrep|log_lik|z_|lin_pred")) |>
  as.data.table()

phylo_eff_cols = grep("phylo_effect", names(draw_dt), value = TRUE)
beta_cols = grep("^beta", names(draw_dt), value = TRUE)
# other_cols = grep("phylo_eff|^beta", x = names(draw_dt), value = TRUE, invert = TRUE)

draw_dt[,phylo_effects := list(list(unlist(.SD))), by = `.draw`, .SDcols = phylo_eff_cols]
draw_dt[,beta          := list(list(.SD)), by = `.draw`, .SDcols = beta_cols]

nested_df = draw_dt[,!..phylo_eff_cols][,!..beta_cols] |>
  tibble::as_tibble()

if (ncol(nested_df$beta[[1]]) != 0) {
  nested_df$beta = purrr::map(nested_df$beta,
                            ~matrix(unlist(.x), ncol = 1))
}

get_ll_mat(nested_df, em, cor_mat, Lcov, Xc, offset_val, metadata$outcome, family)

pieces = list(nested_df = nested_df,
              em = em,
              cor_mat = cor_mat,
              Lcov = Lcov,
              Xc = Xc,
              offset_val = offset_val,
              metadata = metadata,
              family = family)

save(pieces, file = test_path("test_data/gauss_pieces.RData"))
save(gauss_ll, file = test_path("test_data/gauss_ll.RData"))
# save(gauss_ddf, file = test_path("gauss_ddf.RData"))
save(gauss_loo, file = test_path("test_data/gauss_loo.RData"))

# Save pieces needed for gaussian ll tests ----

set.seed(123)

r2 = anpan_pglmm(meta_file       = m2,
                     tree_file       = tr,
                     outcome         = "outcome",
                     covariates      = "covariate",
                     family          = "binomial",
                     bug_name        = "sim_bug",
                     reg_noise       = TRUE,
                     loo_comparison  = TRUE,
                     run_diagnostics = FALSE,
                     refresh         = 500,
                     show_plot_tree  = FALSE,
                     show_post       = FALSE,
                     iter_sampling  = 100)

binom_ll = r2$loo$pglmm_ll_df |> as.matrix()
binom_ddf = r2$pglmm_fit$draws(format = "data.frame")
binom_loo = get_pglmm_loo(binom_ll, binom_ddf)

usethis::use_data(binom_ll, overwrite = TRUE, compress = "xz", internal = TRUE)
usethis::use_data(binom_ddf, overwrite = TRUE, compress = "xz", internal = TRUE)
usethis::use_data(binom_loo, overwrite = TRUE, compress = "xz", internal = TRUE)


