#' Estimate the latent outcome model
#' @description Run a latent outcome model on a distance matrix and a model
#'   matrix.
#' @param model_mat a model matrix. Should not have an intercept column.
#' @param dist_mat a distance matrix.
#' @param iter_sampling the number of post-warmup draws to take
#' @param iter_warmup the number of warmup draws
#' @param sigma_factor higher --> lower prior spread on lm residuals
#' @param algorithm either "mcmc" or "variational"
#' @param verbose verbose setting
#' @param ... other parameters to pass to cmdstanr::sample()
#' @details The variational fit is faster but biased.
#' @return a named list. draw_df is set of MCMC draws of each parameter.
#'   summary_df is the posterior summary. sample_order is a tibble containing
#'   the sample order used by the model. y_scale is the empirically estimated
#'   prior scale of the latent y values and residiuals. term calls is a data
#'   frame of ANOVA-like "calls" for each term in the model formula, indicating
#'   whether any of its associated coefficients are non-zero by 90% interval
#'   excluding 0.
#' @export
anpan_latent_outcome = function(model_mat, dist_mat,
                                iter_sampling = 1000,
                                iter_warmup = 1000,
                                sigma_factor = 5,
                                verbose = TRUE,
                                algorithm = 'mcmc',
                                ...) {

  message("\nThis function is experimental and might be removed!\n")

  Sys.sleep(1)

  if (all(model_mat[,1] == 1)) {
    message("Removing intercept column from model matrix")
    model_mat = model_mat[,-1]
  }

  if (all(colnames(dist_mat) == rownames(model_mat))) {
    message("Samples are in the correct order!")
  } else {
    message("Column names of dist_mat don't match the rownames of the model matrix. Reordering the model matrix.")
    model_mat = model_mat[colnames(dist_mat),]
  }

  Sys.sleep(1)

  A = -1/2 * dist_mat^2

  clust = model_mat |> dist() |> hclust()

  n = nrow(model_mat)

  C = diag(n) - (1/n * matrix(1, nrow = n, ncol = n))
  G = C %*% A %*% C
  y_scale = sqrt(sd(G[upper.tri(G)]))

  y_init = seq(-2*y_scale, 2*y_scale, length.out = n)[sort(clust$order, index.return = TRUE)$ix]

  data_list = list(N = n,
                   d = ncol(model_mat),
                   X = model_mat,
                   dist_mat = dist_mat,
                   y_scale = y_scale,
                   sigma_factor = sigma_factor)

  # Initialize the latent y values at near the best ordering of x values (which
  # presumes a y = 1*x relationship). This is crucial to getting this toy
  # example to work.
  init_list = replicate(4,
                        list(y = y_init + rnorm(n, sd = .1)),
                        simplify = FALSE)

  latent_y_model = cmdstanr::cmdstan_model('inst/stan/latent_y_glm.stan',
                                           stanc_options = list("O1"))

  refresh_setting = if (verbose) 100 else 0

  if (algorithm == "mcmc") {
    latent_y_fit = latent_y_model$sample(data = data_list,
                                         init = init_list,
                                         parallel_chains = 4,
                                         max_treedepth = 12,
                                         adapt_delta = .80,
                                         iter_warmup = iter_warmup,
                                         iter_sampling = iter_sampling,
                                         refresh = refresh_setting,
                                         ...)
  } else {
    init_list = replicate(1,
                          list(y = y_init + rnorm(n, sd = .1)),
                          simplify = FALSE)
    latent_y_fit = latent_y_model$variational(data = data_list,
                                              init = init_list,
                                              refresh = refresh_setting,
                                              ...)
  }

  draw_df = latent_y_fit$draws(format = 'data.frame')
  summ_df = latent_y_fit$summary()

  if (algorithm == "mcmc") {
    sign_df = draw_df |>
      select(`y[1]`, `.chain`) |> # TODO: make this select a y near one end of the spectrum
      group_by(`.chain`) |>
      summarise(orig_sign_y = mean(sign(`y[1]`)))

    if (n_distinct(sign_df$orig_sign_y) > 1) {

      if (!all(sign_df$orig_sign_y %in% c(-1,1))) {
        warning("Detected inconsistent sign among chains AND within chain. The model likely doesn't fit your data well.")
      }

      message("Detected inconsistent sign among chains. Manually adjusting sign of affected chains.")

      Sys.sleep(2)
      joined = draw_df |>
        left_join(sign_df, by = '.chain')

      draw_df = joined |> mutate(across(matches("beta|y"), ~.x * sign(joined$orig_sign_y)))
    }

    summ_df = posterior::summarise_draws(draw_df)
  }

  # Needs better handling of distinct contrast groups rather than treating every
  # column of model_mat as separate

  if (!is.null(attr(model_mat, "contrasts"))) {
    contrast_df = tibble(contrast_i = attr(model_mat, 'assign'),
                         term = names(attr(model_mat, "contrasts"))[contrast_i]) |>
      unique()

    term_calls = summ_df |>
      filter(grepl("beta", variable)) |>
      mutate(contrast_i = attr(model_mat, 'assign')) |>
      right_join(contrast_df, by = 'contrast_i') |>
      group_by(term) |>
      summarise(influential = any(!(q5 < 0 & q95 > 0)))

  } else {

    term_calls = summ_df |>
      filter(grepl("beta", variable)) |>
      mutate(term = colnames(model_mat)) |>
      group_by(term) |>
      summarise(influential = !(q5 < 0 & q95 > 0))
  }

  res = list(draw_df = draw_df,
             sample_order = tibble(sample_id = rownames(model_mat)),
             y_scale = y_scale,
             summary_df = summ_df,
             term_calls = term_calls)

  return(res)
}
