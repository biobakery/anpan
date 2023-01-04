olap_cor_mat_and_meta = function(cor_mat,
                                 meta_file,
                                 covariates,
                                 outcome,
                                 omit_na = FALSE,
                                 trim_pattern = NULL,
                                 verbose = TRUE) {

  if (is.data.frame(meta_file)) {

    meta = meta_file |>
      dplyr::select(dplyr::all_of(c("sample_id", covariates, outcome)))

    if (omit_na) meta = na.omit(meta)

  } else {
    meta = read_meta(meta_file,
                     select_cols = c("sample_id", covariates, outcome),
                     omit_na = omit_na)
  }

  if (!is.null(trim_pattern)) {
    rownames(cor_mat) = gsub(trim_pattern, "", rownames(cor_mat))
    colnames(cor_mat) = gsub(trim_pattern, "", colnames(cor_mat))
  }

  overlapping_samples = intersect(rownames(cor_mat),
                                  meta$sample_id)

  if (!all(rownames(cor_mat) %in% overlapping_samples)) {

    n_tree_olap = sum(rownames(cor_mat) %in% overlapping_samples)

    if (verbose) message(paste0("Dropping ",
                                sum(!(rownames(cor_mat) %in% overlapping_samples)),
                                " tips from the correlation matrix (out of ",
                                length(rownames(cor_mat)),
                                ") not present in the metadata."))
    cor_mat = cor_mat[overlapping_samples,
                      overlapping_samples]

  }

  model_input = meta |>
    dplyr::filter(sample_id %in% overlapping_samples)

  if (nrow(model_input) < nrow(meta) & verbose) {
    message(paste0("Dropping ",
                   nrow(meta) - nrow(model_input),
                   ' samples from the metadata (out of ',
                   nrow(meta),
                   ' in total) not present in the tree.' ))
  }

  res = list(cor_mat,
             model_input)

  return(res)

}

#' Overlap a tree and metadata
#' @description Overlap a tree's leaves with the observations in metadata,
#'   returning an intersected tree and metadata data.table.
#' @return A list with two elements: the tree and metadata, both cut down to
#'   samples that are present in the other.
#' @inheritParams anpan_pglmm
#' @export
olap_tree_and_meta = function(tree_file,
                              meta_file,
                              covariates,
                              outcome,
                              omit_na = FALSE,
                              trim_pattern = NULL,
                              verbose = TRUE) {


  if (class(tree_file) == "phylo") {
    bug_tree = tree_file
  } else{
    bug_tree = ape::read.tree(tree_file)
  }

  if (is.data.frame(meta_file)) {

    meta = meta_file |>
      dplyr::select(dplyr::all_of(c("sample_id", covariates, outcome)))

    if (omit_na) meta = na.omit(meta)

  } else {
    meta = read_meta(meta_file,
                     select_cols = c("sample_id", covariates, outcome),
                     omit_na = omit_na)
  }

  if (!is.null(trim_pattern)) bug_tree$tip.label = gsub(trim_pattern, "", bug_tree$tip.label)

  overlapping_samples = intersect(bug_tree$tip.label,
                                  meta$sample_id)

  if (!all(bug_tree$tip.label %in% overlapping_samples)) {

    n_tree_olap = sum(bug_tree$tip.label %in% overlapping_samples)

    if (verbose) message(paste0("Dropping ",
                                sum(!(bug_tree$tip.label %in% overlapping_samples)),
                                " tips from the tree (out of ",
                                length(bug_tree$tip.label),
                                ") not present in the metadata."))

    bug_tree = ape::drop.tip(bug_tree,
                             tip = bug_tree$tip.label[!(bug_tree$tip.label %in% overlapping_samples)])

  }

  model_input = meta |>
    dplyr::filter(sample_id %in% overlapping_samples)

  if (nrow(model_input) < nrow(meta) & verbose) {
    message(paste0("Dropping ",
                   nrow(meta) - nrow(model_input),
                   ' samples from the metadata (out of ',
                   nrow(meta),
                   ' in total) not present in the tree.' ))
  }

  res = list(bug_tree,
             model_input)

  return(res)
}

#' Get the correlation matrix of a tree
#' @param bug_tree a tree object of class "phylo"
#' @seealso [ape::read.tree()]
#' @export
get_cor_mat = function(bug_tree) {
  cor_mat = ape::vcv.phylo(bug_tree, corr = TRUE)

  # d = sqrt(diag(cov_mat))
  # cor_mat = diag(1/d) %*% cov_mat %*% diag(1/d)
  # dimnames(cor_mat) = dimnames(cov_mat)
  return(cor_mat)
}

safely_chol = purrr::safely(chol)

#' @title Run a phylogenetic generalized linear mixed model
#' @md
#' @param tree_file either a path to a tree file readable by [ape::read.tree()] or an object of
#'   class "phylo" that is already read into R. Ignored if \code{cor_mat} is supplied.
#' @param meta_file either a data frame of metadata or a path to file containing the metadata
#' @param cor_mat a correlation matrix provided as an alternative to a tree.
#' @param trim_pattern optional pattern to trim from tip labels of the tree
#' @param family string giving the name of the distribution of the outcome variable (usually
#'   "gaussian" or "binomial")
#' @param beta_sd prior standard deviation parameters on the normal distribution for each covariate
#'   in the GLM component
#' @param save_object logical indicating whether to save the model fit object
#' @param out_dir if saving, directory where to save
#' @param reg_noise logical indicating whether to regularize the ratio of sigma_phylo to sigma_resid
#'   with a Gamma prior
#' @param reg_gamma_params the shape and rate parameters of the Gamma prior on the noise term ratio.
#'   Default: c(1,2)
#' @param plot_ext extension to use when saving plots
#' @param show_plot_cor_mat show a plot of the correlation matrix derived from the tree
#' @param show_plot_tree show a plot of the tree overlaid with the outcome.
#' @param show_post show a plot of the tree overlaid with the outcome and posterior distribution on
#'   phylogenetic effects.
#' @param show_yrep show a plot of the tree overlaid with the outcome and the posterior predictive
#'   distribution for each observation if plotting the tree
#' @param ... other arguments to pass to [cmdstanr::sample()]
#' @param loo_comparison logical indicating whether to compare the phylogenetic model against a base
#'   model (without the phylogenetic term) using [loo::loo_compare()]
#' @param run_diagnostics logical indicating whether to run [cmdstanr::cmdstan_diagnose()] and
#'   [loo::pareto_k_table()] to check the MCMC and loo diagnostics respectively.
#' @param sigma_phylo_scale standard deviation of half-normal prior on \code{sigma_phylo} for
#'   logistic PGLMMs when \code{family = 'binomial'}. Increasing this value can easily lead to
#'   overfitting.
#' @param int_prior_scale standard deviation of the 0-centered normal prior for the intercept of the
#'   model (with centered covariates)
#' @returns A list containing the model input (in the order passed to the model), estimated
#'   correlation matrix, the pglmm fit object, and (if \code{loo_comparison} is on) the base fit
#'   object and the associated loo objects.
#' @details the tip labels of the tree must be the sample ids from the metadata. You can use the
#'   \code{trim_pattern} argument to automatically trim off any consistent pattern from the tip
#'   labels if necessary.
#'
#'   The dots can be used to pass e.g. parallel_chains=4 to make the chains run in parallel.
#'
#'   The prior for the intercept is a normal distribution centered on the mean of the outcome
#'   variable with a standard deviation of 3*sd(outcome variable)
#'
#'   The default error distribution for the outcome is "gaussian". You could change this to a
#'   phylogenetic logistic regression by changing \code{family} to "binomial" for example.
#'
#'   It's normal to see some warnings during warmup, particularly about "Scale vector is inf".
#'
#'   This function tries to get the \code{bug_name} argument from the tree file, but if it's not
#'   provided you may need to set it yourself.
#'
#'   If using \code{beta_sd} with a categorical predictor with >2 levels, only specify a single
#'   element in beta_sd. This appropriate element will get repeated as necessary.
#'
#'   Common cmdstanr options that one might want to pass in via ... include \code{refresh = 500}
#'   (show fewer MCMC progress updates), \code{adapt_delta = .98} (avoid divergences at the cost of
#'   possibly needing more iterations to get convergence) and \code{parallel_chains = 4} (run the
#'   MCMC chains in parallel).
#'
#'   If you want to use the PGLMM log-likelihood data frame with functions from the \code{loo}
#'   package, you'll need to convert it to a matrix with \code{as.matrix()}. It's converted to a
#'   tibble internally so that it prints nicely.
#'
#'   If \code{int_prior_scale} isn't specified, it defaults to 1 for binary outcomes and 1 standard
#'   deviation of the outcome for gaussian outcomes.
#' @examples
#' meta = data.frame(x = rnorm(100), sample_id = paste0("t", 1:100))
#' tr = ape::rtree(100)
#' anpan_pglmm(meta, tr,
#' outcome = "x",
#' iter_sampling = 10, # Use more in practice!
#' iter_warmup = 10,
#' show_plot_cor_mat = FALSE,
#' show_plot_tree = FALSE)
#' @inheritParams anpan
#' @seealso [anpan_pglmm_batch()], [loo::loo()], [cmdstanr::sample()]
#' @export
anpan_pglmm = function(meta_file,
                       tree_file = NULL,
                       cor_mat = NULL,
                       outcome,
                       covariates = NULL,
                       out_dir = NULL,
                       trim_pattern = NULL,
                       bug_name = NULL,
                       omit_na = FALSE,
                       family = "gaussian",
                       show_plot_cor_mat = TRUE,
                       show_plot_tree = TRUE,
                       show_post = TRUE,
                       show_yrep = FALSE,
                       save_object = FALSE,
                       verbose = TRUE,
                       loo_comparison = TRUE,
                       run_diagnostics = TRUE,
                       reg_noise = TRUE,
                       reg_gamma_params = c(1,2),
                       plot_ext = "pdf",
                       beta_sd = NULL,
                       int_prior_scale = 1,
                       sigma_phylo_scale = 0.333,
                       ...) {

  n_steps = 2 + loo_comparison + run_diagnostics

  if (verbose) message(paste0("(1/", n_steps, ") Checking inputs."))

  if (save_object && is.null(out_dir)) stop("To save the fit you must specify an output directory")

  if (save_object && !dir.exists(out_dir)) {
    message("Creating specified output directory...")
    dir.create(out_dir)
  }

  if (reg_noise && family != "gaussian") {
    warning("You can't regularize the noise ratio with non-gaussian outcomes. Setting reg_noise = FALSE.")
    reg_noise = FALSE
  }

  if (is.null(cor_mat)) {
    cor_mat_provided = FALSE
  } else {
    cor_mat_provided = TRUE
  }

  if (is.null(cor_mat)) {
    olap_list = olap_tree_and_meta(tree_file,
                                   meta_file,
                                   covariates,
                                   outcome,
                                   omit_na,
                                   trim_pattern = trim_pattern,
                                   verbose)

    if (is.null(olap_list[[1]])) {
      stop("Couldn't find any overlapping samples between the sample_id column in the metadata and the tips of the first tree. Maybe you need to set the trim_pattern argument?")
    }

    bug_tree = olap_list[[1]]
    model_input = olap_list[[2]]
  } else {
    olap_list = olap_cor_mat_and_meta(cor_mat,
                                      meta_file,
                                      covariates,
                                      outcome,
                                      omit_na,
                                      trim_pattern,
                                      verbose)
    cor_mat = olap_list[[1]]
    model_input = olap_list[[2]]
  }

  if (family == "binomial" && dplyr::n_distinct(model_input[[outcome]]) != 2) {
    stop('family == "binomial" but couldn\'t find 2 distinct outcomes in the outcome variable.')
  }

  if (family == "binomial" && is.character(model_input[[outcome]])) {
    stop('family == "binomial" but the type of the outcome variable is character. Please change to either 0/1, FALSE/TRUE, or an appropriately ordered factor.')
  }

  if (family == "binomial" && is.factor(model_input[[outcome]])) {
    message("family == \"binomial\" and outcome variable is a factor. Converting to 0/1 with the following mapping:")
    mapping_df = data.frame(input_levels = levels(model_input[[outcome]]),
                        mapped_values = as.numeric(unique(sort(model_input[[outcome]]))) - 1)

    message(paste0(capture.output(mapping_df),
                   sep = "\n"))

    model_input[[outcome]] = as.numeric(model_input[[outcome]]) - 1

  }

  if (!omit_na && nrow(na.omit(model_input)) < nrow(model_input)) {
    stop("omit_na == FALSE but NAs present in metadata. Either set omit_na = TRUE or fix the metadata.")
  }

  if (cor_mat_provided) {

    chol_res = safely_chol(cor_mat)

    if (!is.null(chol_res$error)) stop("Could not compute the Cholesky factorization of the correlation matrix. It's probably not positive definite up to numerical precision. Try olap_tree_and_meta() and get_cor_mat() to examine the correlation matrix directly.")

    Lcov = t(chol_res$result)

  } else {
    cor_mat = get_cor_mat(bug_tree)

    chol_res = safely_chol(cor_mat)

    if (!is.null(chol_res$error)) stop("Could not compute the Cholesky factorization of the correlation matrix. It's probably not positive definite up to numerical precision. Try olap_tree_and_meta() and get_cor_mat() to examine the correlation matrix directly.")

    Lcov = t(chol_res$result)
  }

  if (!(class(tree_file) == "phylo") && is.null(bug_name) && !is.null(tree_file)) {
    bug_name = get_bug_name(tree_file,
                            remove_pattern = ".tre$|.tree$")
  } else {
    bug_name = basename(tempfile()) #
  }

  # Prepare model formulae
  if (!is.null(covariates)) {
    cov_str = paste(covariates, collapse = " + ")
    base_formula = as.formula(paste0(outcome, " ~ ", cov_str))
  } else {
    cov_str = NULL
    base_formula = as.formula(paste0(outcome, " ~ 1"))
  }

  if (!exists("chains")) {
    chains = 4
  }

  # TODO figure out how to set priors as a function of family and # and structure of covariates
  if (family == "gaussian") {
    outcome_mean = mean(model_input[[outcome]])
    outcome_sd = sd(model_input[[outcome]])

    base_path = system.file("stan", "cont_no_phylo_term.stan",
                            package = 'anpan',
                            mustWork = TRUE)

    init_list = replicate(chains,
                          list(sigma_resid = 1, sigma_phylo = 1),
                          simplify = FALSE)

    base_init = replicate(chains,
                          list(sigma_resid = 1),
                          simplify = FALSE)

    if (reg_noise == FALSE) {
      model_path = system.file("stan", "cont_pglmm.stan",
                               package = 'anpan',
                               mustWork = TRUE)
    } else {
      model_path = system.file("stan", "cont_pglmm_noise_reg.stan",
                               package = 'anpan',
                               mustWork = TRUE)

      if (any(reg_gamma_params < 0)) {
        stop("Parameters for the regularizing gamma distribution (reg_gamma_params) must be positive.")
      }
    }
  } else if (family == "binomial") {

    init_list = replicate(chains,
                          list(sigma_phylo = 1),
                          simplify = FALSE)

    base_init = NULL

    model_path = system.file("stan", "bin_pglmm.stan",
                             package = 'anpan',
                             mustWork = TRUE)
    base_path = system.file("stan", "bin_no_phylo_term.stan",
                             package = 'anpan',
                             mustWork = TRUE)
  }

  if (!is.null(out_dir)) {
    pglmm_dir = file.path(out_dir, bug_name, 'pglmm_fit');
    base_dir = file.path(out_dir, bug_name, 'base_fit');
    dir.create(pglmm_dir, recursive = TRUE)
    dir.create(base_dir, recursive = TRUE)
  } else {
    pglmm_dir = NULL
    base_dir = NULL
  }

  pglmm_model = cmdstanr::cmdstan_model(stan_file = model_path,
                                        quiet = TRUE)

  if (loo_comparison) {
    base_model = cmdstanr::cmdstan_model(stan_file = base_path,
                                         quiet = TRUE)
  }

  model_input$sample_id = factor(model_input$sample_id,
                                 levels = rownames(Lcov))
  model_input = model_input |>
    arrange(sample_id)

  if (show_plot_cor_mat) {
    mat_for_plot = cor_mat
    mat_for_plot = mat_for_plot[model_input$sample_id, model_input$sample_id]

    if (verbose) message("Plotting correlation matrix...")
    p = plot_cor_mat(mat_for_plot,
                     bug_name)

    if (verbose) print(p)
    if (!is.null(out_dir)) {
      ggsave(p,
             filename = file.path(out_dir, paste0(bug_name, "_cor_mat.", plot_ext)),
             width = 6, height = 5)
    }
  }

  # data_list preparation section ----

  model_mat = model.matrix(base_formula, data = model_input)

  if (family == "gaussian") {

    data_list = list(N = nrow(model_input),
                     Y = model_input[[outcome]],
                     K = ncol(model_mat),
                     X = model_mat,
                     Lcov = Lcov,
                     int_mean = outcome_mean,
                     int_prior_scale = ifelse(!is.null(int_prior_scale), int_prior_scale, outcome_sd),
                     resid_scale = outcome_sd)

    if (reg_noise) data_list$reg_gamma_params = reg_gamma_params

  } else {

    data_list = list(N = nrow(model_input),
                     Y = model_input[[outcome]],
                     K = ncol(model_mat),
                     X = model_mat,
                     Lcov = Lcov,
                     sigma_phylo_scale = sigma_phylo_scale,
                     int_prior_scale = ifelse(!is.null(int_prior_scale), int_prior_scale, 1))

  }

  Xc = matrix(nrow = data_list$N,
              ncol = data_list$K - 1)

  mx = colMeans(data_list$X)[-1]

  if (ncol(Xc) > 0) {
    # can't remember how to use seq_along here
    for (i in 2:data_list$K) {
      Xc[,i-1] = data_list$X[,i] - mx[i-1]
    }
  }

  if (!is.null(beta_sd) &&
      length(beta_sd[as.vector(attr(model_mat, "assign")[-1])]) != (data_list$K - 1)) stop("The number of covariates provided doesn't match the length of the provided beta_sd argument.")

  if (is.null(beta_sd)) {
    # TODO write this to a log file...
    if (length(covariates) > 0) {
      prior_df = data.table(linear_term = colnames(data_list$X)[-1],
                            prior_sd = round(apply(Xc, 2, sd), digits = 3))
      message("Prior scale on covariate effects aren't specified. Setting to 1 standard deviation for each centered covariate. These values are:\n")
      message(paste0(capture.output(prior_df),
                     sep = "\n"))
      message("\n\nIt would be better to set the beta_sd argument based on scientific background knowledge.")
    }
    if (ncol(Xc) > 0) {
      data_list$beta_sd = apply(Xc, 2, sd)
    } else {
      data_list$beta_sd = numeric()
    }
  } else {
    data_list$beta_sd = beta_sd[as.vector(attr(model_mat, "assign")[-1])]
  }

  # model fitting section ----
  if (verbose) message(paste0("(2/", n_steps, ") Fitting model(s)."))

  pglmm_fit = pglmm_model$sample(data = data_list,
                                 chains = chains,
                                 init = init_list,
                                 output_dir = pglmm_dir,
                                 ...)

  if (loo_comparison) {
    base_fit = base_model$sample(data = data_list,
                                 chains = chains,
                                 init = base_init,
                                 output_dir = base_dir,
                                 ...)
  }

  # plotting section ----
  if (cor_mat_provided && is.null(tree_file) && show_plot_tree) {
    message("show_plot_tree = TRUE but no tree_file was supplied. Setting show_plot_tree = FALSE")
    show_plot_tree = FALSE
  }

  if (show_plot_tree) {

    p = plot_outcome_tree(tree_file,
                          meta_file,
                          covariates = covariates,
                          outcome = outcome,
                          omit_na = omit_na,
                          trim_pattern = trim_pattern,
                          verbose = FALSE)

    if (!is.null(out_dir)) {
      ggsave(p, filename = file.path(out_dir, paste0(bug_name, "_outcome_tree.", plot_ext)),
             width = 12, height = 8)
    }

    if (verbose) print(p)
  }

  if (show_post) {
    p_post = plot_tree_with_post(tree_file,
                                 meta_file,
                                 covariates = covariates,
                                 outcome = outcome,
                                 omit_na = omit_na,
                                 verbose = FALSE,
                                 fit = pglmm_fit,
                                 trim_pattern = trim_pattern,
                                 labels = levels(model_input$sample_id))

    if (!is.null(out_dir)) {
      ggsave(p_post, filename = file.path(out_dir, paste0(bug_name, "_posterior_tree.", plot_ext)),
             width = 12, height = 8)
    }

    if (verbose) print(p_post)
  }

  if (show_yrep) {
    p_post_pred = plot_tree_with_post_pred(tree_file,
                                           meta_file,
                                           covariates = covariates,
                                           outcome = outcome,
                                           omit_na = omit_na,
                                           verbose = FALSE,
                                           fit = pglmm_fit,
                                           trim_pattern = trim_pattern,
                                           labels = levels(model_input$sample_id))

    if (!is.null(out_dir)) {
      ggsave(p_post_pred, filename = file.path(out_dir, paste0(bug_name, "_posterior_predictive_tree.", plot_ext)),
             width = 12, height = 8)
    }

    if (verbose) print(p_post_pred)
  }

  # loo section ----
  if (loo_comparison) {
    if (verbose) message(paste0("(3/", n_steps, ") Evaluating loo comparison."))
    if (verbose) message("- 0/2 preparing loo inputs")

    if (!reg_noise && family == "gaussian") {
      warning("loo-based model comparison is unstable and/or anti-conservative without regularized noise. Do not trust the results if there are bad Pareto k diagnostic values.")
    }

    draw_dt = pglmm_fit$draws(format = "data.frame") |>
      tibble::as_tibble() |>
      select(-tidyselect::matches("std_phylo|yrep|log_lik|z_|lin_pred")) |>
      as.data.table()

    phylo_eff_cols = grep("phylo_effect", names(draw_dt), value = TRUE)
    beta_cols = grep("^beta", names(draw_dt), value = TRUE)
    # other_cols = grep("phylo_eff|^beta", x = names(draw_dt), value = TRUE, invert = TRUE)

    draw_dt[,phylo_effects := list(list(unlist(.SD))), by = `.draw`, .SDcols = phylo_eff_cols]
    draw_dt[,beta          := list(list(.SD)), by = `.draw`, .SDcols = beta_cols]

    draw_df = draw_dt[,!..phylo_eff_cols][,!..beta_cols] |>
      tibble::as_tibble()

    if (ncol(draw_df$beta[[1]]) != 0) {
      draw_df$beta = purrr::map(draw_df$beta,
                                ~matrix(unlist(.x), ncol = 1))
    }

    fit_summary = pglmm_fit$summary() |>
      tibble::as_tibble()

    # These are decent starting places for offset terms
    effect_means = fit_summary |>
      filter(grepl("^phylo_effect", variable)) |>
      pull(mean)

    ll_mat = get_ll_mat(draw_df,
                        max_i = nrow(draw_df),
                        effect_means = effect_means,
                        cor_mat = cor_mat,
                        Lcov = Lcov,
                        Xc = Xc,
                        Y = data_list$Y,
                        family = family,
                        verbose = verbose)

    pglmm_loo = get_pglmm_loo(ll_mat, draw_df)

    base_loo = base_fit$loo()

    message("loo comparison: ")
    comparison = loo::loo_compare(list(pglmm_fit = pglmm_loo,
                                       base_fit  = base_loo))
    print(comparison)

    if (rownames(comparison)[1] == 'pglmm_fit') {
      p1 = "The phylogenetic model seems to fit better, "
    } else {
      p1 = "The phylogenetic model seems to fit worse, "
    }

    if (abs(comparison[2,1] / comparison[2,2]) > 2) {
      p2 = "and the difference seems clear (more than 2 standard errors difference in ELPD)."
    } else {
      p2 = "but the difference doesn't seem clear (less than 2 standard errors difference in ELPD)."
    }

    if (abs(comparison[2,1] / comparison[2,2]) > 2 && abs(comparison[2,1]) < 4) {
      p3 = " However the ELPD difference is less than 4, so the difference is small."
    } else {
      p3 = NULL
    }

    message(paste0(p1, p2, p3))

  } else {

    # outcome_signal = NULL
    # hyp = NULL
    base_fit = NULL
    pglmm_loo = NULL
    base_loo = NULL
    comparison = NULL
    ll_mat = NULL
  }

  if (run_diagnostics) {
    if (verbose) message(paste0("(", n_steps, "/", n_steps, ") Running diagnostics:"))

    pglmm_fit$cmdstan_diagnose() # No need for print(), it already does itself

    if (loo_comparison) {
      print(loo::pareto_k_table(pglmm_loo))
    }
  }

  if (!is.null(out_dir)) {
    save(model_input, cor_mat,
         file = file.path(out_dir, paste0(bug_name, "_inputs.RData")))
  }

  if (save_object) {
    # V This is what to use once the pglmm_fit is done with cmdstanr
    pglmm_fit$save_object(file = file.path(out_dir, paste0(bug_name, "_pglmm_fit.RDS")))
    base_fit$save_object(file = file.path(out_dir, paste0(bug_name, "_base_fit.RDS")))
  }

  if (!save_object && !is.null(pglmm_dir)) {
    unlink(pglmm_dir, recursive = TRUE)
    unlink(base_dir, recursive = TRUE)
  }

  list(model_input = model_input,
       cor_mat     = cor_mat,
       pglmm_fit   = pglmm_fit,
       base_fit    = base_fit,
       loo         = list(pglmm_loo    = pglmm_loo,
                          pglmm_ll_df  = tibble::as_tibble(ll_mat),
                          base_loo     = base_loo,
                          comparison   = comparison))

}

safely_anpan_pglmm = purrr::safely(anpan_pglmm)

#' Run PGLMMs on a batch of tree files
#' @description This function fits phylogenetic generalized linear mixed models on a batch of tree
#'   files, using the same outcome and covariate arguments.
#' @param tree_dir string giving the path to a directory of tree files
#' @param seed random seed to pass to furrr_options()
#' @details \code{tree_dir} must contain ONLY tree files readable by ape::read.tree()
#'
#'   If any trees cause an error while fitting, these are saved into a data frame in a file
#'   \code{pglmm_errors.RData} in the output directory.
#'
#'   The Stan model fitting can't be parallelized via futures, so the most effective way to
#'   parallelize the model fitting AND the importance weight calculations is a nested future
#'   topology (e.g. \code{plan(list(sequential, tweak(multisession, workers = 4)))} ) and set
#'   parallel_chains = 4 . This will run sequentially over the trees, running the model fits with
#'   parallel chains for each tree, then compute the importance weights in the future multisession
#'   for each tree.
#' @returns a tibble listing results for each tree file in input directory that fit successfully.
#'   Columns give the number of leaves on the tree, diagnostic values, loo comparison values,
#'   formatted input data, correlation matrices, PGLMM and "base" model fits, and loo objects (in
#'   list columns where appropriate).
#' @inheritParams anpan_pglmm
#' @seealso \code{\link[ape:read.tree]{ape::read.tree}},
#'   \code{\link[ape:write.tree]{ape::write.tree}}, \code{\link[=anpan_pglmm]{anpan_pglmm()}}
#' @export
anpan_pglmm_batch = function(meta_file,
                             tree_dir,
                             outcome,
                             covariates = NULL,
                             out_dir = NULL,
                             trim_pattern = NULL,
                             omit_na = FALSE,
                             family = "gaussian",
                             show_plot_cor_mat = TRUE,
                             show_plot_tree = TRUE,
                             save_object = FALSE,
                             verbose = TRUE,
                             loo_comparison = TRUE,
                             run_diagnostics = FALSE,
                             reg_noise = TRUE,
                             plot_ext = "pdf",
                             show_yrep = FALSE,
                             show_post = TRUE,
                             reg_gamma_params = c(1,2),
                             beta_sd = NULL,
                             sigma_phylo_scale = 0.333,
                             seed = 123,
                             ...) {

  n_steps = 3

  # Checking inputs ---------------------------------------------------------

  if (verbose) message(paste0("(1/", n_steps, ") Checking inputs."))

  dot_list = list(...)
  if ("parallel_chains" %in% names(dot_list) && attr(future::plan(), "class")[2] != "sequential") {
    stop("Don't specify parallel_chains with a non-sequential plan. Use future::plan() to parallelize over trees instead. A two-level future topology may be used to fit the model with parallel chains parallelize loo evaluation within each tree e.g. plan(list(sequential, multisession)).")
  }

  call = match.call()

  fn_call_string = paste0(gsub(', (?!")',
                               ",\n            ",
                               as.character(enquote(call))[2],
                               perl = TRUE),
                          "\n")

  if (verbose & !interactive()) message(paste0("Now running:\n\n", fn_call_string))

  if (!is.null(out_dir) && !dir.exists(out_dir)) {
    if (verbose) message("* Creating output directory.")
    dir.create(out_dir)

    cat(fn_call_string,
        file = file.path(out_dir, "anpan_pglmm_batch_call.txt"),
        sep = "\n")
  }

  tree_files = list.files(tree_dir,
                          full.names = TRUE)

  # Get the model paths
  if (family == "gaussian") {
    base_path = system.file("stan", "cont_no_phylo_term.stan",
                            package = 'anpan',
                            mustWork = TRUE)
    if (reg_noise == FALSE) {
      model_path = system.file("stan", "cont_pglmm.stan",
                               package = 'anpan',
                               mustWork = TRUE)
    } else {
      model_path = system.file("stan", "cont_pglmm_noise_reg.stan",
                               package = 'anpan',
                               mustWork = TRUE)
    }
  } else if (family == "binomial") {
    model_path = system.file("stan", "bin_pglmm.stan",
                             package = 'anpan',
                             mustWork = TRUE)
    base_path = system.file("stan", "bin_no_phylo_term.stan",
                            package = 'anpan',
                            mustWork = TRUE)
  }

  # Check for overlaps in the first file, just to make sure the entire batch doesn't suffer from a
  # missing trim pattern.

  olap_list = olap_tree_and_meta(tree_files[1],
                                 meta_file,
                                 covariates,
                                 outcome,
                                 omit_na,
                                 trim_pattern = trim_pattern,
                                 verbose = FALSE)

  if (is.null(olap_list[[1]])) {
    stop("Couldn't find any overlapping samples between the sample_id column in the metadata and the tips of the first tree. Maybe you need to set the trim_pattern argument?")
  }


  # initial compilation -----------------------------------------------------

  if (verbose) message(paste0("(2/", n_steps, ") Performing initial model compilation(s)."))

  # Compile them once here so that they don't get compiled inside future_map()
  pglmm_model = cmdstanr::cmdstan_model(stan_file = model_path,
                                        quiet = TRUE)

  if (loo_comparison) {
    base_model = cmdstanr::cmdstan_model(stan_file = base_path,
                                         quiet = TRUE)
  }


  # Fitting PGLMMs ----------------------------------------------------------

  if (verbose) message(paste0("(3/", n_steps, ") Fitting PGLMMs."))

  p = progressr::progressor(along = tree_files)

  # global_list = c('meta_file', # TODO figure out if passing this is necessary. If so, find out what else needs to be added...
  #                 'outcome'  ,
  #                 'out_dir' ,
  #                 'trim_pattern' ,
  #                 'covariates' ,
  #                 'omit_na' ,
  #                 'family' ,
  #                 'show_plot_cor_mat' ,
  #                 'show_plot_tree' ,
  #                 'save_object' ,
  #                 'loo_comparison',
  #                 'reg_noise',
  #                 'plot_ext' ,
  #                 'show_yrep')

  safe_results = furrr::future_map(tree_files,
                                   function(.x) {
                                     res = safely_anpan_pglmm(tree_file = .x,
                                                              meta_file = meta_file,
                                                              outcome = outcome,
                                                              out_dir = out_dir,
                                                              trim_pattern = trim_pattern,
                                                              covariates = covariates,
                                                              omit_na = omit_na,
                                                              family = family,
                                                              show_plot_cor_mat = show_plot_cor_mat,
                                                              show_plot_tree = show_plot_tree,
                                                              save_object = save_object,
                                                              verbose = FALSE,
                                                              loo_comparison = loo_comparison,
                                                              run_diagnostics = run_diagnostics,
                                                              reg_noise = reg_noise,
                                                              beta_sd = beta_sd,
                                                              reg_gamma_params = reg_gamma_params,
                                                              plot_ext = plot_ext,
                                                              show_yrep = show_yrep,
                                                              show_post = show_post,
                                                              sigma_phylo_scale = sigma_phylo_scale,
                                                              ...)
                                     p()
                                     return(res)},
                                   .options = furrr::furrr_options(seed = seed))

  safe_res_df = purrr::transpose(safe_results) |>
    as_tibble() |>
    mutate(input_file = basename(tree_files)) |>
    relocate(input_file)

  errors = safe_res_df |>
    dplyr::filter(purrr::map_lgl(error, ~!is.null(.x)))

  if (nrow(errors) > 0) {
    save(errors,
         file = file.path(out_dir, 'pglmm_errors.RData'))

    if (verbose) message(paste0("* There were ", nrow(errors), " tree files that failed to fit. You can see the error messages they produced in pglmm_errors.RData in the output directory."))
  }

  worked = safe_res_df |>
    dplyr::filter(purrr::map_lgl(error, ~is.null(.x)))

  res_df = dplyr::bind_cols(worked["input_file"],
                            as_tibble(purrr::transpose(worked$result))) |>
    mutate(n                      = purrr::map_int(model_input, nrow),
           prop_bad_pareto_k      = purrr::map_dbl(loo, ~mean(.x$pglmm_loo$diagnostics$pareto_k > .7)),
           n_divergences          = purrr::map_dbl(pglmm_fit, ~sum(as.matrix(.x$sampler_diagnostics()[,,"divergent__"]))),
           best_model             = gsub("_fit", "", purrr::map_chr(loo, ~rownames(.x$comparison)[1])),
           elpd_diff              = purrr::map_dbl(loo, ~.x$comparison[2,1]),
           elpd_se                = purrr::map_dbl(loo, ~.x$comparison[2,2])) |>
    dplyr::select(input_file, n:elpd_se, model_input:loo)

  return(res_df)

}


#' Fit a subject-wise PGLMM
#' @description If you have a dataset where multiple samples come from the same individual, running
#'   a PGLMM on all the samples without accounting for this will return a strong spurious signal
#'   because samples from the same individual will (usually) be right next to each other on the tree
#'   and all show the same outcome value. This function uses a subject-sample map input to aggregate
#'   samples by subject, derive a subject-wise correlation matrix, then run a PGLMM on that.
#' @param subject_sample_map a data frame giving a mapping between sample_id (which must match the
#'   leaves of the tree up to the trim_pattern) and subject_id
#' @param meta_file either a data frame or a file. Can provide covariate and outcome variables for
#'   either subjects or samples.
#' @details The \code{meta_file} must contain at least one column named "sample_id" or "subject_id".
#'   If the metadata is inferred to be provided by sample, representative covariate and outcome
#'   variable values are selected for each subject in the manner described in
#'   \code{?anpan_repeated_measures()}
#'
#'
#' @inheritParams anpan_pglmm
#' @seealso \code{\link[=anpan_pglmm]{anpan_pglmm()}} \code{\link[=anpan_repeated_measures]{anpan_repeated_measures()}}
#' @export
anpan_subjectwise_pglmm = function(tree_file,
                                   meta_file,
                                   subject_sample_map,
                                   outcome,
                                   covariates = NULL,
                                   out_dir = NULL,
                                   trim_pattern = NULL,
                                   omit_na = FALSE,
                                   family = "gaussian",
                                   show_plot_cor_mat = TRUE,
                                   show_plot_tree = TRUE,
                                   show_post = TRUE,
                                   show_yrep = FALSE,
                                   save_object = FALSE,
                                   verbose = TRUE,
                                   loo_comparison = TRUE,
                                   reg_noise = TRUE,
                                   reg_gamma_params = c(1,2),
                                   plot_ext = "pdf",
                                   beta_sd = NULL,
                                   sigma_phylo_scale = 0.333,
                                   ...) {
  # tree tip labels must match sample_id values already

  if (!is.data.table(subject_sample_map)) subject_sample_map = as.data.table(subject_sample_map)

  olap_tree_map = olap_tree_and_meta(tree_file, subject_sample_map, covariates = "subject_id",
                                     trim_pattern = trim_pattern, outcome = NULL, verbose = FALSE)

  tree = olap_tree_map[[1]]
  subject_sample_map = olap_tree_map[[2]]

  if (!is.data.table(subject_sample_map)) subject_sample_map = as.data.table(subject_sample_map)

  n_subj = dplyr::n_distinct(subject_sample_map$subject_id)

  cor_mat = tree |>
    ape::vcv.phylo(corr = TRUE)

  subj_cor_mat = diag(n_subj)
  dimnames(subj_cor_mat) = replicate(2,
                                     unique(subject_sample_map$subject_id),
                                     simplify = FALSE)

  subj_df = combn(unique(subject_sample_map$subject_id),
                  2) |>
    t() |>
    as.data.table()

  # subj_df[subject_sample_map, on = c("V1" = "subject_id"), allow.cartesian = TRUE][subject_sample_map, on = c("V2" = "subject_id"), allow.cartesian = TRUE, ]

  cor_df = subj_df |>
    left_join(subject_sample_map, by = c("V1" = "subject_id")) |>
    left_join(subject_sample_map, by = c("V2" = "subject_id"), suffix = c("1", "2")) |>
    mutate(cor_val = cor_mat[cbind(sample_id1, sample_id2)])

  avg_cor_df = cor_df[,.(avg_cor = mean(cor_val)), by = .(V1, V2)]

  subj_cor_mat[cbind(avg_cor_df$V1,
                     avg_cor_df$V2)] = avg_cor_df$avg_cor
  subj_cor_mat[cbind(avg_cor_df$V2,
                     avg_cor_df$V1)] = avg_cor_df$avg_cor

  eig_subj_cor_mat = eigen(subj_cor_mat)

  if (any(eig_subj_cor_mat$values < 0)) {
    stop("The derived subject-wise correlation matrix was not positive definite.")
  }

  # https://stats.stackexchange.com/questions/165194/using-correlation-as-distance-metric-for-hierarchical-clustering
  subj_dist = sqrt(2 * (1-subj_cor_mat))

  subj_tree = ape::nj(subj_dist) |>
    ape::ladderize()

  tree_dend = ggdendro::dendro_data(subj_tree |> phylogram::as.dendrogram())

  ordered_subj_cor_mat = subj_cor_mat[tree_dend$labels$label, tree_dend$labels$label]

  # Get the subjectwise metadata, either from the meta_file or aggregating the meta_file

  if (is.character(meta_file) && file.exists(meta_file)) {
    metadata = fread(meta_file,
                     header = TRUE,
                     showProgress = FALSE)
  } else if (is.data.frame(meta_file)) {
    metadata = as.data.table(meta_file)
  } else {
    stop("The provided meta_file doesn't seems to be neither a file nor a data frame.")
  }

  if ("sample_id" %in% names(metadata)) {
    output_cols = c(covariates, outcome, "subject_id")

    if ("subject_id" %in% names(metadata)) metadata$subject_id = NULL

    metadata = dplyr::left_join(metadata, subject_sample_map, by = "sample_id")[,..output_cols] |>
      unique() |>
      dplyr::group_by(subject_id) |>
      dplyr::summarise(dplyr::across(.cols = dplyr::all_of(c(covariates, outcome)),
                                     summarise_metadata_variable)) |>
      dplyr::rename(sample_id = subject_id) |>
      as.data.table()
  } else if ("subject_id" %in% names(metadata)) {
    metadata = metadata |>
      dplyr::select(dplyr::all_of(c("subject_id", outcome, covariates))) |>
      dplyr::rename(sample_id = subject_id)
  } else {
    stop("Couldn't find subject_id or sample_id in the metadata")
  }

  # plot_half_cor_mat(ordered_subj_cor_mat)
  anpan_pglmm(meta_file         = metadata,
              tree_file         = subj_tree,
              cor_mat           = ordered_subj_cor_mat,
              outcome           = outcome,
              covariates        = covariates,
              out_dir           = out_dir,
              trim_pattern      = NULL,
              bug_name          = bug_name,
              omit_na           = omit_na,
              family            = family,
              show_plot_cor_mat = show_plot_cor_mat,
              show_plot_tree    = show_plot_tree,
              show_post         = show_post,
              show_yrep         = show_yrep,
              save_object       = save_object,
              verbose           = verbose,
              loo_comparison    = loo_comparison,
              reg_noise         = reg_noise,
              reg_gamma_params  = reg_gamma_params,
              plot_ext          = plot_ext,
              beta_sd           = beta_sd,
              sigma_phylo_scale = sigma_phylo_scale,
              ...)


}
