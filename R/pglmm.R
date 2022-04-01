olap_tree_and_meta = function(tree_file,
                              meta_file,
                              covariates,
                              outcome,
                              omit_na = FALSE,
                              trim_pattern,
                              verbose = TRUE) {

  if (class(tree_file) == "phylo") {
    bug_tree = tree_file
  } else{
    bug_tree = ape::read.tree(tree_file)
  }

  if (is.data.frame(meta_file)) {

    meta = meta_file %>%
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

  model_input = meta %>%
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

#' @title Run a phylogenetic generalized linear mixed model
#' @md
#' @param tree_file either a path to a tree file readable by [ape::read.tree()]
#'   or an object of class "phylo" that is already read into R
#' @param meta_file either a data frame of metadata or a path to file containing
#'   the metadata
#' @param trim_pattern optional pattern to trim from tip labels of the tree
#' @param family string giving the name of the distribution of the outcome
#'   variable (usually "gaussian" or "binomial")
#' @param save_object logical indicating whether to save the model fit object
#' @param out_dir if saving, directory where to save
#' @param reg_noise logical indicating whether to regularize the ratio of
#'   sigma_phylo to sigma_resid with a Gamma(1.33,2) prior
#' @param ... other arguments to pass to [cmdstanr::sample()]
#' @param loo_comparison logical indicating whether to compare the phylogenetic
#'   model against a base model (without the phylogenetic term) using
#'   [loo::loo_compare()]
#' @details the tip labels of the tree must be the sample ids from the metadata.
#'   You can use the \code{trim_pattern} argument to automatically trim off any
#'   consistent pattern from the tip labels if necessary.
#'
#'   The dots can be used to pass e.g. parallel_chains=4 to make the chains run
#'   in parallel.
#'
#'   The prior for the intercept is a normal distribution centered on the mean
#'   of the outcome variable with a standard deviation of 3*sd(outcome variable)
#'
#'   The default error distribution for the outcome is "gaussian". You could
#'   change this to a phylogenetic logistic regression by changing \code{family}
#'   to "binomial" for example.
#'
#'   It's normal to see some warnings during warmup, particularly about "Scale
#'   vector is inf".
#' @inheritParams anpan
#' @seealso [loo::loo()], [cmdstanr::sample()]
#' @export
anpan_pglmm = function(meta_file,
                       tree_file,
                       outcome,
                       out_dir = NULL,
                       trim_pattern = NULL,
                       covariates = NULL,
                       omit_na = FALSE,
                       family = "gaussian",
                       plot_cor_mat = TRUE,
                       save_object = FALSE,
                       verbose = TRUE,
                       loo_comparison = TRUE,
                       reg_noise = TRUE,
                       ...) {

  n_steps = ifelse(loo_comparison,
                   3, 2)

  if (verbose) message(paste0("(1/", n_steps, ") Checking inputs."))

  if (save_object && is.null(out_dir)) stop("To save the fit you must specify an output directory")

  if (save_object && !dir.exists(out_dir)) {
    message("Creating specified output directory...")
    dir.create(out_dir)
  }

  if (reg_noise && family != "gaussian") {
    stop("You can't regularize the noise ratio with non-gaussian outcomes.")
  }

  olap_list = olap_tree_and_meta(tree_file,
                                 meta_file,
                                 covariates,
                                 outcome,
                                 omit_na,
                                 trim_pattern,
                                 verbose)

  bug_tree = olap_list[[1]]
  model_input = olap_list[[2]]

  cov_mat = ape::vcv.phylo(bug_tree)

  d = sqrt(diag(cov_mat)) # Not sure if needed
  cor_mat = diag(1/d) %*% cov_mat %*% diag(1/d)
  dimnames(cor_mat) = dimnames(cov_mat)

  if (plot_cor_mat) {

    if (!(class(tree_file) == "phylo")) {
      bug_name = get_bug_name(tree_file,
                              remove_pattern = ".tre$|.tree$")
    } else {
      bug_name = NULL # no way to get the bug name in this case
    }

    p = make_cor_mat_plot(cor_mat,
                          bug_name)
    if (verbose) message("Plotting correlation matrix...")

    print(p)
    if (!is.null(out_dir)) {
      if (is.null(bug_name)) bug_name = basename(tempfile())

      ggsave(p,
             filename = file.path(out_dir, paste0(bug_name, "_cor_mat.png")),
             width = 6, height = 5)
    }
  }

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

  pglmm_model = cmdstanr::cmdstan_model(stan_file = model_path,
                                        quiet = TRUE)

  if (loo_comparison) {
    base_model = cmdstanr::cmdstan_model(stan_file = base_path,
                                         quiet = TRUE)
  }

  Lcov = t(chol(cor_mat))

  model_input$sample_id = factor(model_input$sample_id,
                                 levels = rownames(Lcov))
  model_input = model_input %>%
    arrange(sample_id)

  if (family == "gaussian") {
    data_list = list(N = nrow(model_input),
                     Y = model_input[[outcome]],
                     K = length(covariates) + 1,
                     X = model.matrix(base_formula, data = model_input),
                     Lcov = Lcov,
                     int_mean = outcome_mean,
                     resid_scale = outcome_sd)
  } else {
    data_list = list(N = nrow(model_input),
                     Y = model_input[[outcome]],
                     K = length(covariates) + 1,
                     X = model.matrix(base_formula, data = model_input),
                     Lcov = Lcov)
  }

  if (verbose) message(paste0("(2/", n_steps, ") Fitting model(s)."))



  pglmm_fit = pglmm_model$sample(data = data_list,
                                 chains = chains,
                                 init = init_list,
                                 ...)

  if (loo_comparison) {
    base_fit = base_model$sample(data = data_list,
                                 chains = chains,
                                 init = base_init,
                                 ...)
  }

  # TODO throw out the "Intercept" parameter and rename "b_Intercept" as
  # appropriate (and move it to the top...)

  if (loo_comparison) {
    if (verbose) message(paste0("(3/", n_steps, ") Evaluating loo comparison."))
    pglmm_loo = pglmm_fit$loo()
    base_loo = base_fit$loo()

    message("loo comparison: ")
    comparison = loo::loo_compare(list(pglmm_fit = pglmm_loo,
                                       base_fit = base_loo))
    print(comparison)

    if (rownames(comparison)[1] == 'pglmm_fit') {
      p1 = "The phylogenetic model seems to fit better, "
    } else {
      p1 = "The phylogenetic model seems to fit worse, "
    }

    if (abs(comparison[2,1] / comparison[2,2]) > 2) {
      p2 = "and the difference seems substantial (more than 2 standard errors difference in ELPD)."
    } else {
      p2 = "but the difference doesn't seem substantial (less than 2 standard errors difference in ELPD)."
    }

    message(paste0(p1, p2))

  } else {
    # outcome_signal = NULL
    # hyp = NULL
    base_fit = NULL
    pglmm_loo = NULL
    base_loo = NULL
    comparison = NULL
  }

  if (save_object) {
    # V This is what to use once the pglmm_fit is done with cmdstanr
    pglmm_fit$save_object(file = file.path(out_dir, paste0(tree_name, "_pglmm_fit.RDS")))
    base_fit$save_object(file = file.path(out_dir, paste0(tree_name, "_base_fit.RDS")))
  }

  list(pglmm_fit = pglmm_fit,
       base_fit = base_fit,
       loo = list(pglmm_loo = pglmm_loo, base_loo = base_loo, comparison = comparison))

}

# TODO anpan_pglmm_batch()
