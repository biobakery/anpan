## code to prepare integrated log-likelihood matrices for the binomial
## integrated log-likelihood. Uses the updated Cpp version since it turns out to
## be more numerically accurate on top of being faster.

# pak::pak("biobakery/anpan@abd7c94")
library(anpan)
library(tibble)
library(dplyr)
library(data.table)
library(testthat)

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

# Save pieces needed for binomial ll tests ----

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

binom_ddf = r2$pglmm_fit$draws(format = "data.frame")

em = r2$pglmm_fit$summary(variables = "phylo_effect",
                              mean = mean)$mean
family = "binomial"

draw_dt = r2$pglmm_fit$draws(format = "data.frame") |>
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
cor_mat = r2$cor_mat

chol_res = anpan:::safely_chol(cor_mat)
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

offset_val = rep(0, nrow(m2))

binom_ll = anpan:::get_ll_mat(nested_df, em, cor_mat, Lcov, Xc, offset_val, m2$outcome, family,
                              verbose = FALSE)
binom_loo = anpan:::get_pglmm_loo(binom_ll, binom_ddf)

pieces = list(nested_df = nested_df,
              em = em,
              cor_mat = cor_mat,
              Lcov = Lcov,
              Xc = Xc,
              offset_val = offset_val,
              metadata = m2,
              family = family)

save(pieces, file = test_path("test_data/binom_pieces2.RData"),
     compress = "zstd")

save(binom_ll, file = test_path("test_data/binom_ll2.RData"),
     compress = "zstd")

save(binom_loo, file = test_path("test_data/binom_loo2.RData"),
     compress = "zstd")
