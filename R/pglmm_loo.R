# Standard loo doesn't work with PGLMMs because of the high number of random effects + high number
# of constraints. You have to integrate out the phylogenetic effect for the likelihood of each
# observation to get the "integrated importance weights" from section 3.6.1 from here: Vehtari, Aki,
# Tommi Mononen, Ville Tolvanen, Tuomas Sivula, and Ole Winther. “Bayesian Leave-One-Out
# Cross-Validation Approximations for Gaussian Latent Variable Models,” n.d., 38.

# The integral is relatively simple, but needs to be evaluated numerically for non-Gaussian
# outcomes. It has to integrate the likelihood on the identity scale (not the log scale), so you
# have to use an offset trick like LogSumExp(). LogIntExp() I guess.

# You have to figure out the distribution of each phylogenetic effect conditional on the others.
# This involves a bit of linear algebra that can be mostly pre-computed. See the wikipedia link
# below.

# See also the thread where I had to ask for help lol:
# https://discourse.mc-stan.org/t/integrated-loo-with-a-pglmm/27271

# run loo on the log-likelihood matrix
get_pglmm_loo = function(ll_mat, draw_df) {
  loo::loo(x = ll_mat,
           r_eff = loo::relative_eff(exp(ll_mat),
                                     chain_id = draw_df$`.chain`))
}


# For each posterior iteration, compute the log-likelihood of each observation.
get_ll_mat = function(draw_df, max_i, effect_means, cor_mat, Lcov, Xc, Y, family, verbose = TRUE) {

  n_obs = length(effect_means)

  # precompute some arrays that are re-usable on each iteration.
  # See this link for notation. sigma12x22_inv = sigma12 %*% solve(sigma22) for each observation j
  # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
  sigma12x22_inv_arr = array(data = NA,
                             dim = c(1, n_obs - 1, n_obs))

  cor21_arr = array(dim = c(n_obs - 1, 1, n_obs))

  if (verbose) message("- 1/2 precomputing conditional covariance arrays")
  # Using future_map is safe here because even if anpan_pglmm_batch is run in a
  # future, nested futures run sequentially.
  p = progressr::progressor(steps = n_obs)

  cor_mat_inv = chol2inv(t(Lcov))
  # ^ Stan and R have different conventions on where the triangle should be for a cholesky factor...

  arr_list = furrr::future_map(1:n_obs,
                               function(.x) {
                                 res = precompute_arrays(j = .x, cor_mat = cor_mat, cor_mat_inv = cor_mat_inv)
                                 p()
                                 return(res)},
                               .options = furrr::furrr_options(globals = c("cor_mat", "cor_mat_inv")))

  for (j in 1:n_obs) {
    sigma12x22_inv_arr[,,j] = arr_list[[j]][[1]]
    cor21_arr[,,j] = arr_list[[j]][[2]]
  }

  # set up the progressr and list over posterior iterations
  if (verbose) message("- 2/2 computing integrated importance weights for loo CV")

  p = progressr::progressor(along = 1:max_i)

  draw_split = draw_df[1:max_i,] |>
    dplyr::group_split(`.draw`)

  # future map over posterior iterations

  if (family == "gaussian") {
    global_list = c("effect_means",
                    "cor_mat",
                    "Xc", "Y",
                    "sigma12x22_inv_arr",
                    "cor21_arr",
                    "family",
                    "log_lik_terms_i",
                    "log_lik_i_j_gaussian")
  } else {
    global_list = c("effect_means",
                    "cor_mat",
                    "Xc", "Y",
                    "sigma12x22_inv_arr",
                    "cor21_arr",
                    "family",
                    "log_lik_terms_i",
                    "inv_logit",
                    "log_lik_i_j_logistic",
                    "vec_integrand_logistic",
                    "p",
                    'safely_integrate')

  }

  ll_list = furrr::future_map(draw_split,
                       function(.x) {
                         res = log_lik_terms_i(i_df = .x,
                                               effect_means = effect_means,
                                               cor_mat = cor_mat,
                                               Xc = Xc, Y = Y,
                                               sigma12x22_inv_arr = sigma12x22_inv_arr,
                                               cor21_arr = cor21_arr,
                                               family = family)
                         p()
                         return(res)
                       },
                       .options = furrr::furrr_options(globals = global_list))

  # put the result into a matrix
  res = ll_list |>
    unlist() |>
    matrix(nrow = length(ll_list),
           byrow = TRUE)

  colnames(res) = paste("log_lik[", 1:n_obs, "]", sep = "")

  return(res)
}

woodbury_s22_inv = function(cor_mat_inv, cor_mat, j) {
  # https://en.wikipedia.org/wiki/Woodbury_matrix_identity

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

# compute multivariate normal conditional components
precompute_arrays = function(j, cor_mat, cor_mat_inv) {

  n = nrow(cor_mat)

  ord = c(j, (1:n)[-j])

  r12 = cor_mat[ord, ord][1,-1, drop = FALSE]
  r22_inv = woodbury_s22_inv(cor_mat_inv, cor_mat, j)

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
                           cor21_arr,
                           family = 'gaussian') {

  p = length(effect_means)
  cov_mat = i_df$sigma_phylo^2 * cor_mat

  if (ncol(Xc) > 0) {
    covariate_term = Xc %*% matrix(i_df$beta[[1]], ncol = 1)
  } else {
    covariate_term = 0
  }

  lm_means = covariate_term + i_df$centered_cov_intercept

  # j = index over observations
  if (family == 'gaussian') {
    # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
    j_df = data.table(j = seq_len(p),
                      l = c(lm_means),
                      sigma12x22_inv = lapply(seq_len(p), function(.x) matrix(sigma12x22_inv_arr[,,.x], nrow = 1)),
                      sigma21        = lapply(seq_len(p), function(.x) matrix(i_df$sigma_phylo^2 * cor21_arr[,,.x], ncol = 1)),
                      effects_mj     = lapply(seq_len(p), function(.x) matrix(i_df$phylo_effects[[1]][-.x], ncol = 1)),
                      sigma_resid    = i_df$sigma_resid,
                      yj             = Y,
                      effect_mean_j  = effect_means,
                      cov_mat_jj     = diag(cov_mat)) |>
      dplyr::transmute(m1 = mapply(function(.x, .y) c(.x %*% .y),
                                   sigma12x22_inv, effects_mj),
                       s1 = mapply(function(.x, .y, .z) sqrt(.x - c(.y %*% .z)),
                                   cov_mat_jj, sigma12x22_inv, sigma21),
                       l  = l,
                       yj = yj,
                       s2 = sigma_resid)

    ll_ij_fun = log_lik_i_j_gaussian
  } else {
    j_df = data.table(j              = seq_len(p),
                      lm_mean        = c(lm_means),
                      sigma12x22_inv = lapply(seq_len(p), function(.x) matrix(sigma12x22_inv_arr[,,.x], nrow = 1)),
                      sigma21        = lapply(seq_len(p), function(.x) matrix(i_df$sigma_phylo^2 * cor21_arr[,,.x], ncol = 1)),
                      effects_mj     = lapply(seq_len(p), function(.x) matrix(i_df$phylo_effects[[1]][-.x], ncol = 1)),
                      yj             = Y,
                      effect_mean_j  = effect_means,
                      cov_mat_jj     = diag(cov_mat),
                      sqrt2pi        = sqrt(2*pi))
    ll_ij_fun = log_lik_i_j_logistic
  }

  # Map over all p leaves
  purrr::pmap_dbl(j_df,
                  ll_ij_fun)
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

log_lik_i_j_gaussian = function(m1, s1, l, yj, s2) {

  # Get the coefficients of the quadratic log-likelihood of the two components of a gaussian PGLMM
  # m1, s1 = mean, sd of the PGLMM correlation term
  # yj, l, s2 = outcome, linear fit term, sigma_resid for the normally distributed outcome
  # Sum of two parabolas = a new parabola. Get the coefficients of that.

  abc = c(-1/2 * (1/s1^2 + 1/s2^2),
          m1 / s1^2 + (yj - l) / s2^2,
          -1/2*(m1^2/s1^2)
            - 1/2 * (l^2 - 2*l*yj + yj^2)/s2^2
            - log(2*pi*s1^2)/2
            - log(2*pi*s2^2)/2)


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

  # This is for a gaussian outcome variable.

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

safely_integrate = purrr::safely(integrate)
safely_optim = purrr::safely(optim)

log_lik_i_j_logistic = function(j, lm_mean, sigma12x22_inv, sigma21,
                                effects_mj, # effects minus j
                                yj,
                                effect_mean_j, cov_mat_jj, sqrt2pi) {

  p = length(effects_mj) + 1

  # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
  sigma_bar_j = sqrt(cov_mat_jj - c(sigma12x22_inv %*% sigma21))

  mu_bar_j = c(sigma12x22_inv %*% (effects_mj))
  #^ a = the other phylo effects from iteration i, mu2 = 0

  offset_j = vec_integrand_logistic(phylo_effect_vec = effect_mean_j,
                                    mu_bar_j         = mu_bar_j,
                                    sigma_bar_j      = sigma_bar_j,
                                    yj               = yj,
                                    lm_term          = lm_mean,
                                    offset_term      = 0,
                                    sqrt2pi          = sqrt2pi,
                                    log              = TRUE)

  offset_check = abs(offset_j) > 25

  if (!offset_check) {
    int_res = safely_integrate(vec_integrand_logistic,
                               lower = -Inf, upper = Inf,
                               mu_bar_j         = mu_bar_j,
                               sigma_bar_j      = sigma_bar_j,
                               yj               = yj,
                               lm_term          = lm_mean,
                               offset_term      = offset_j,
                               sqrt2pi          = sqrt2pi,
                               log              = FALSE)
  }

  if (offset_check ||
      !is.null(int_res$error) ||
      !is.finite(int_res$result$value) ||
      int_res$result$value == 0 ||
      int_res$result$value < 1e-4 ||
      int_res$result$value > 100) {
    opt_res = safely_optim(par = mu_bar_j,
                           fn = vec_integrand_logistic,
                           method = "L-BFGS-B",
                           control = list(fnscale = -1),
                           mu_bar_j         = mu_bar_j,
                           sigma_bar_j      = sigma_bar_j,
                           yj               = yj,
                           lm_term          = lm_mean,
                           offset_term      = 0,
                           sqrt2pi          = sqrt2pi,
                           log              = TRUE)

    if (!is.null(opt_res$error)) {
      opt_res = safely_optim(par = effect_mean_j,
                             fn = vec_integrand_logistic,
                             method = "L-BFGS-B",
                             control = list(fnscale = -1),
                             mu_bar_j         = mu_bar_j,
                             sigma_bar_j      = sigma_bar_j,
                             yj               = yj,
                             lm_term          = lm_mean,
                             offset_term      = 0,
                             sqrt2pi          = sqrt2pi,
                             log              = TRUE)
    }

    opt_res = opt_res$result

    offset_j = opt_res$value

    int_res = safely_integrate(vec_integrand_logistic,
                               lower = -Inf, upper = Inf,
                               mu_bar_j         = mu_bar_j,
                               sigma_bar_j      = sigma_bar_j,
                               yj               = yj,
                               lm_term          = lm_mean,
                               offset_term      = offset_j,
                               sqrt2pi          = sqrt2pi,
                               log              = FALSE)

    if (!is.null(int_res$error) ||
        int_res$result$value == 0 ||
        int_res$result$value < 1e-4 ||
        int_res$result$value > 100) {
      # If it still fails, it's a really sharp integral. Take a parabolic
      # approximation about the optimum to find integration limits.

      integrand_pts = mu_bar_j + c(-.001, 0, .001)
      offset_j_terms = vec_integrand_logistic(phylo_effect_vec = integrand_pts,
                                              mu_bar_j         = mu_bar_j,
                                              sigma_bar_j      = sigma_bar_j,
                                              yj               = yj,
                                              lm_term          = lm_mean,
                                              offset_term      = 0,
                                              sqrt2pi          = sqrt2pi,
                                              log              = TRUE)

      # Set up the linear system
      A = c(integrand_pts) |> sapply(\(x) c(x^2, x, 1)) |> t()

      abc = solve(A, offset_j_terms)
      ll_max = opt_res$par

      # https://www.wolframalpha.com/input?i=solve+a*%28x%2Bd%29%5E2+%2B+b*%28x%2Bd%29+%2B+c+%3D+k-10+for+d
      # d is the +/- delta value over which the log-likelihood of the phylogenetic effect[ij] will decrease by 10
      d = -(2*(abc[1]*ll_max) + abc[2] + sqrt(abc[2]^2 - 4*abc[1]*(abc[3]-opt_res$value+10))) / (2*abc[1])

      int_res = safely_integrate(vec_integrand_logistic,
                                 lower = ll_max - d,
                                 upper = ll_max + d,
                                 mu_bar_j         = mu_bar_j,
                                 sigma_bar_j      = sigma_bar_j,
                                 yj               = yj,
                                 lm_term          = lm_mean,
                                 offset_term      = offset_j,
                                 sqrt2pi          = sqrt2pi,
                                 log              = FALSE)
    }
  }

  ll_ij = log(int_res$result$value) + offset_j

  ll_ij

}

inv_logit = function(x) 1 / (1 + exp(-x))

vec_integrand_logistic = function(phylo_effect_vec,
                                  mu_bar_j, sigma_bar_j,      # phylo term components
                                  yj, lm_term,   # LM term components, NO SIGMA_RESID NOW!
                                  offset_term,
                                  sqrt2pi,
                                  log = FALSE) {
  # This function for the integrand is vectorized because that's what
  # integrate() needs. The log argument keeps the result on the log scale, which
  # can be helpful for finding an offset for numerical precision.

  # This is for a binomial outcome variable.

  # phylo_term = dnorm(phylo_effect_vec,
  #                    mean = mu_bar_j,
  #                    sd = sigma_bar_j,
  #                    log = TRUE)

  # Doing the log normal density manually is faster than dnorm(log=TRUE). This is NOT true for dbinom().
  phylo_term = -1/2 * ((phylo_effect_vec - mu_bar_j) / sigma_bar_j)^2 - log(sigma_bar_j * sqrt2pi)

  model_mean = c(lm_term) + phylo_effect_vec

  fit_term = dbinom(x = yj,
                    size = 1,
                    prob = inv_logit(model_mean),
                    log = TRUE)

  res = phylo_term + fit_term - offset_term

  if (!log) res = exp(res)

  return(res)
}

