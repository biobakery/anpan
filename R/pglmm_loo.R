# Standard loo doesn't work with PGLMMs because of the high number of random
# effects + high number of constraints. You have to use the integrated loo
# method from section 3.6.1 from here: Vehtari, Aki, Tommi Mononen, Ville
# Tolvanen, Tuomas Sivula, and Ole Winther. “Bayesian Leave-One-Out
# Cross-Validation Approximations for Gaussian Latent Variable Models,” n.d.,
# 38.

# See also the thread where I had to ask for help lol:
# https://discourse.mc-stan.org/t/integrated-loo-with-a-pglmm/27271


get_pglmm_loo = function(ll_mat, draw_df) {
  # run loo on the log-likelihood matrix
  loo::loo(x = ll_mat,
           r_eff = loo::relative_eff(exp(ll_mat),
                                     chain_id = draw_df$`.chain`))
}

# For each posterior iteration, compute the log-likelihood of the
# observations.
get_ll_mat = function(draw_df, max_i, effect_means, cor_mat, Xc, Y) {

  n_obs = length(effect_means)

  # precompute some arrays that are re-usable on each iteration.
  # See this link for notation. sigma12x22_inv = sigma12 %*% solve(sigma22) for each observation j
  # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
  sigma12x22_inv_arr = array(data = NA,
                             dim = c(1, n_obs - 1, n_obs))

  cor21_arr = array(dim = c(n_obs - 1, 1, n_obs))

  arr_list = furrr::future_map(1:n_obs, precompute_arrays,
                        cor_mat = cor_mat, n_obs = n_obs,
                        .progress = TRUE)

  for (j in 1:n_obs) {
    sigma12x22_inv_arr[,,j] = arr_list[[j]][[1]]
    cor21_arr[,,j] = arr_list[[j]][[2]]
  }

  # set up the progressr and list over posterior iterations

  p = progressr::progressor(along = 1:max_i)

  draw_split = draw_df[1:max_i,] %>%
    dplyr::group_split(`.draw`)

  # future map over posterior iterations

  ll_list = furrr::future_map(draw_split,
                       function(.x) {
                         p()
                         log_lik_terms_i(i_df = .x,
                                         effect_means = effect_means,
                                         cor_mat = cor_mat,
                                         Xc = Xc, Y = Y,
                                         sigma12x22_inv_arr = sigma12x22_inv_arr,
                                         cor21_arr = cor21_arr)
                       })

  # put the result into a matrix
  res = ll_list |>
    unlist() |>
    matrix(nrow = length(ll_list),
           byrow = TRUE)

  colnames(res) = paste("log_lik[", 1:n_obs, "]", sep = "")

  return(res)
}

# compute multivariate normal conditional components
precompute_arrays = function(j, cor_mat, n_obs) {


  ord = c(j, (1:n_obs)[-j])

  r12 = cor_mat[ord, ord][1,-1, drop = FALSE]
  r22_inv = solve(cor_mat[-j,-j])

  sigma12x22_inv_arr_j = r12 %*% r22_inv
  cor21_arr_j = t(r12)

  return(list(sigma12x22_inv_arr_j,
              cor21_arr_j))
}

# compute the log likelihood for all observations in a single posterior iteration i
log_lik_terms_i = function(i_df,
                           effect_means,
                           cor_mat,
                           Xc, Y,
                           sigma12x22_inv_arr,
                           cor21_arr) {

  p = length(effect_means)
  cov_mat = i_df$sigma_phylo^2 * cor_mat

  lm_means = Xc %*% matrix(i_df$beta[[1]], ncol = 1) + i_df$Intercept

  # j = index over observations
  j_df = tibble(j              = 1:p,
                lm_mean        = c(lm_means),
                sigma12x22_inv = purrr::map(j, ~matrix(sigma12x22_inv_arr[,,.x], nrow = 1)),
                sigma21        = purrr::map(j, ~matrix(i_df$sigma_phylo^2 * cor21_arr[,,.x], ncol = 1)),
                effects_mj     = purrr::map(j, ~matrix(i_df$phylo_effects[[1]][-.x], ncol = 1)),
                sigma_resid    = i_df$sigma_resid,
                yj             = Y,
                effect_mean_j  = effect_means,
                cov_mat_jj     = purrr::map_dbl(j, ~cov_mat[.x,.x]))

  # Map over all p leaves
  purrr::pmap_dbl(j_df,
                  log_lik_i_j)
}

# Compute the log-likelihood of a single observation at a single iteration
# i = posterior iteration
# j = index over observations / phylogenetic effects
log_lik_i_j = function(j, lm_mean, sigma12x22_inv, sigma21,
                       effects_mj, # effects minus j
                       sigma_resid, yj,
                       effect_mean_j, cov_mat_jj) {

  p = length(effects_mj) + 1
  # ord = c(j, (1:p)[-j])

  # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
  sigma_bar_j = cov_mat_jj - c(sigma12x22_inv %*% sigma21)

  mu_bar_j = c(sigma12x22_inv %*% (effects_mj))
  #^ a = the other phylo effects from iteration i, mu2 = 0

  # Evaluate the log likelihood at 3 points then fit a parabola. This method
  # likely won't work for non-gaussian outcomes.
  parab_x_points = effect_mean_j + c(-.25,0,.25)
  offset_j_terms = vec_integrand(phylo_effect_vec = parab_x_points,
                                 mu_bar_j         = mu_bar_j,
                                 sigma_bar_j      = sigma_bar_j,
                                 sigma_resid      = sigma_resid,
                                 yj               = yj,
                                 lm_term          = lm_mean,
                                 offset_term      = 0,
                                 log              = TRUE)

  A = parab_x_points |>
    sapply(\(x) c(x^2, x, 1)) |>
    t()

  abc = solve(A, offset_j_terms) # max at -b/2a

  ll_max = -abc[2] / (2*abc[1])

  # Compute the integrated log-likelihood value directly from the parabola.
  # https://www.wolframalpha.com/input?i=integrate+exp%28a*x%5E2+%2B+b*x+%2B+c%29+from+-inf+to+inf
  lik = sqrt(pi) * exp(abc[3] - (abc[2]^2)/(4*abc[1])) / sqrt(-abc[1])

  return(log(lik))
}

vec_integrand = function(phylo_effect_vec,
                         mu_bar_j, sigma_bar_j,      # phylo term components
                         sigma_resid, yj, lm_term,   # LM term components
                         offset_term,
                         log = FALSE) {
  # This function for the integrand is vectorized because that's what
  # integrate() needs. The log argument keeps the result on the log scale, which
  # can be helpful for finding an offset for numerical precision.

  # This is for a gaussian outcome variable. It will likely need to be changed
  # for non-gaussian outcomes.

  phylo_term = dnorm(phylo_effect_vec,
                     mean = mu_bar_j,
                     sd = sqrt(sigma_bar_j),
                     log = TRUE)

  model_mean = c(lm_term) + phylo_effect_vec

  fit_term = dnorm(x = yj,
                   mean = model_mean,
                   sd = sigma_resid,
                   log = TRUE)

  res = phylo_term + fit_term - offset_term

  if (!log) res = exp(res)

  return(res)
}
