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
#'   sigma_phylo to sigma_resid with a Gamma(1,2) prior
#' @param plot_ext extension to use when saving plots
#' @param show_plot_cor_mat show a plot of the correlation matrix derived from
#'   the tree
#' @param show_plot_tree show a plot of the tree overlaid with the outcome
#' @param show_yrep logical indicating whether to show the posterior predictive
#'   distribution for each observation if plotting the tree
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
#'
#'   This function tries to get the \code{bug_name} argument from the tree file,
#'   but if it's not provided you may need to set it yourself.
#' @inheritParams anpan
#' @seealso [loo::loo()], [cmdstanr::sample()]
#' @export
anpan_pglmm = function(meta_file,
                       tree_file,
                       outcome,
                       covariates = NULL,
                       out_dir = NULL,
                       trim_pattern = NULL,
                       bug_name = NULL,
                       omit_na = FALSE,
                       family = "gaussian",
                       show_plot_cor_mat = TRUE,
                       show_plot_tree = TRUE,
                       save_object = FALSE,
                       verbose = TRUE,
                       loo_comparison = TRUE,
                       reg_noise = TRUE,
                       plot_ext = "pdf",
                       show_yrep = TRUE,
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
    stop("You can't regularize the noise ratio with non-gaussian outcomes. Set reg_noise = FALSE to use a model without regularized noise.")
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

  d = sqrt(diag(cov_mat))
  cor_mat = diag(1/d) %*% cov_mat %*% diag(1/d)
  dimnames(cor_mat) = dimnames(cov_mat)

  if (!(class(tree_file) == "phylo") && is.null(bug_name)) {
    bug_name = get_bug_name(tree_file,
                            remove_pattern = ".tre$|.tree$")
  } else {
    bug_name = basename(tempfile()) #
  }

  if (show_plot_cor_mat) {
    p = plot_cor_mat(cor_mat,
                          bug_name)
    if (verbose) message("Plotting correlation matrix...")

    if (verbose) print(p)
    if (!is.null(out_dir)) {
      ggsave(p,
             filename = file.path(out_dir, paste0(bug_name, "_cor_mat.", plot_ext)),
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
                                 output_dir = pglmm_dir,
                                 ...)
  # TODO throw out the "Intercept" parameter and rename "b_Intercept" as
  # appropriate (and move it to the top...)

  if (loo_comparison) {
    base_fit = base_model$sample(data = data_list,
                                 chains = chains,
                                 init = base_init,
                                 output_dir = base_dir,
                                 ...)
  }

  if (show_plot_tree) {
    if (!show_yrep) {
      p = plot_tree(tree_file,
                    meta_file,
                    covariates = covariates,
                    outcome = outcome,
                    omit_na = omit_na,
                    verbose = FALSE)
    } else {
      p = plot_tree_with_post_pred(tree_file,
                                   meta_file,
                                   covariates = covariates,
                                   outcome = outcome,
                                   omit_na = omit_na,
                                   verbose = FALSE,
                                   fit = pglmm_fit,
                                   labels = levels(model_input$sample_id))
    }

    if (!is.null(out_dir)) {
      ggsave(p, filename = file.path(out_dir, paste0(bug_name, "_tree.", plot_ext)),
             width = 12, height = 8)
    }

    if (verbose) print(p)
  }

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
    pglmm_fit$save_object(file = file.path(out_dir, paste0(bug_name, "_pglmm_fit.RDS")))
    base_fit$save_object(file = file.path(out_dir, paste0(bug_name, "_base_fit.RDS")))
  }

  if (!is.null(pglmm_dir)){
    unlink(pglmm_dir, recursive = TRUE)
    unlink(base_dir, recursive = TRUE)
  }

  list(pglmm_fit = pglmm_fit,
       base_fit = base_fit,
       loo = list(pglmm_loo = pglmm_loo, base_loo = base_loo, comparison = comparison))

}

safely_anpan_pglmm = purrr::safely(anpan_pglmm)

#' Run PGLMMs on a batch of tree files
#' @description This function fits phylogenetic generalized linear mixed models
#'   on a batch of tree files, using the same outcome and covariate arguments.
#' @param tree_dir string giving the path to a directory of tree files
#' @param seed random seed to pass to furrr_options()
#' @details \code{tree_dir} must contain ONLY tree files readable by ape::read.tree()
#' @inheritParams anpan_pglmm
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
                             reg_noise = TRUE,
                             plot_ext = "pdf",
                             show_yrep = TRUE,
                             seed = 123,
                             ...) {

  n_steps = 3

  # Checking inputs ---------------------------------------------------------

  if (verbose) message(paste0("(1/", n_steps, ") Checking inputs."))

  dot_list = list(...)
  if ("parallel_chains" %in% names(dot_list)) {
    stop("Don't specify parallel_chains. Use future::plan() to parallelize instead.")
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

    readr::write_lines(fn_call_string,
                       file = file.path(out_dir, "anpan_pglmm_batch_call.txt"))
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


  # initial compilation -----------------------------------------------------

  if (verbose) message(paste0("(2/", n_steps, ") Performing initial model compilation."))

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

  safe_results = furrr::future_map(tree_files,
                                   function(.x) {
                                     p()
                                     safely_anpan_pglmm(tree_file = .x,
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
                                                        reg_noise = reg_noise,
                                                        plot_ext = plot_ext,
                                                        show_yrep = show_yrep,
                                                        parallel_chains = 1,
                                                        refresh = 0,
                                                        ...)},
                                   .options = furrr::furrr_options(seed = seed))

  safe_res_df = purrr::transpose(safe_results) %>%
    as_tibble() %>%
    mutate(input_file = basename(tree_files)) %>%
    relocate(input_file)

  errors = safe_res_df %>%
    dplyr::filter(purrr::map_lgl(error, ~!is.null(.x)))

  if (nrow(errors) > 0) {
    save(errors,
         file = file.path(out_dir, 'pglmm_errors.RData'))

    if (verbose) message(paste0("* There were ", nrow(errors), " tree files that failed to fit. You can see the error messages they produced in pglmm_errors.RData in the output directory."))
  }

  worked = safe_res_df %>%
    dplyr::filter(purrr::map_lgl(error, ~is.null(.x)))

  res_df = bind_cols(worked["input_file"],
                     as_tibble(purrr::transpose(worked$result)))

  return(worked)

}
