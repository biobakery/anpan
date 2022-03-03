#' @title Run a phylogenetic generalized linear mixed model
#' @param tree_file path to a .tre file
#' @param trim_pattern optional pattern to trim from tip labels of the tree
#' @param family family object for the error distribution
#' @param save_object logical indicating whether to save the model fit object
#' @param out_dir if saving, directory where to save
#' @param ... other arguments to pass to brms::brm
#' @param test_signal compute the phylogenetic signal
#' @details the tip labels of the tree must be the sample ids from the metadata.
#'   You can use the \code{trim_pattern} argument to automatically trim off any
#'   consistent pattern from the tip labels if necessary. The dots can be used
#'   to pass cores=4 to brms to make the chains run in parallel.
#'
#'   The prior for the intercept is a normal distribution centered on the mean
#'   of the outcome variable with a standard deviation of 3*sd(outcome variable)
#'
#'   The default error distribution for the outcome is \code{stats::gaussian()}.
#'   You could change this to a phylogenetic logistic regression by changing
#'   \code{family} to \code{brms::bernoulli()} for example.
#' @inheritParams anpan
#' @export
anpan_pglmm = function(meta_file,
                       tree_file,
                       outcome,
                       out_dir = NULL,
                       trim_pattern = NULL,
                       covariates = NULL,
                       omit_na = FALSE,
                       family = stats::gaussian(),
                       plot_cor_mat = TRUE,
                       save_object = FALSE,
                       verbose = TRUE,
                       test_signal = TRUE,
                       ...) {

  if (save_object && is.null(out_dir)) stop("To save the fit you must specify an output directory")

  meta = read_meta(meta_file,
                   select_cols = c("sample_id", covariates, outcome),
                   omit_na = omit_na)

  tree_name = basename(tree_file) %>%
    stringr::str_replace_all(pattern = c("RAxML_bestTree\\.|\\.StrainPhlAn3|\\.tre$|\\.tree$"),
                             replacement = "")

  bug_tree = ape::read.tree(tree_file)

  overlapping_samples = intersect(bug_tree$tip.label,
                                  meta$sample_id)
  model_input = meta %>%
    dplyr::filter(sample_id %in% overlapping_samples)

  if (!is.null(trim_pattern)) bug_tree$tip.label = gsub(trim_pattern, "", bug_tree$tip.label)

  cov_mat = ape::vcv.phylo(bug_tree)

  d = sqrt(diag(cov_mat)) # Not sure if needed
  cor_mat = diag(1/d) %*% cov_mat %*% diag(1/d)
  dimnames(cor_mat) = dimnames(cov_mat)

  if (plot_cor_mat) {
    bug_name = get_bug_name(tree_file,
                            remove_pattern = ".tre$|.tree$")
    p = make_cor_mat_plot(cor_mat,
                          bug_name)
    if (verbose) message("Plotting correlation matrix...")

    print(p)
    if (!is.null(out_dir)) {
      ggsave(p,
             filename = file.path(out_dir, paste0(bug_name, "_cor_mat.png")))
    }
  }

  if (nrow(model_input) < nrow(meta) & verbose) {
    message(paste0("Using ", nrow(model_input), ' samples present in the tree out of ', nrow(meta), " present in the metadata." ))
  }

  if (!is.null(covariates)) {
    cov_str = paste0(paste(covariates, collapse = " + "), " + ")
  } else {
    cov_str = NULL
  }

  model_formula = as.formula(paste0(outcome, " ~ ", cov_str, "(1|gr(sample_id, cov = cor_mat))"))
  if (is.null(cov_str)){
    base_formula = as.formula(paste0(outcome, " ~ 1"))
  } else {
    base_formula = as.formula(paste0(outcome, " ~ ", cov_str))
  }


  # TODO figure out how to set priors as a function of family and # and structure of covariates
  if (family$family == "gaussian") {
    n_mean = mean(model_input[[outcome]])
    n_sd = 3*sd(model_input[[outcome]])
    prior_vec = c(
      brms::prior_string(paste0("normal(",n_mean, ",", n_sd, ")"), "Intercept"),
      brms::prior(student_t(3, 0, 20), "sd"),
      brms::prior(student_t(3, 0, 20), "sigma")
    )
  } else if (family$family == "bernoulli") {
    prior_vec = NULL
  }

  pglmm_fit = brms::brm(formula = model_formula,
                        data = model_input,
                        family = family,
                        data2 = list(cor_mat = cor_mat),
                        prior = prior_vec,
                        backend = "cmdstanr",
                        ...)

  base_fit = brms::brm(formula = base_formula,
                       data = model_input,
                       family = family,
                       prior = prior_vec[-2],
                       backend = "cmdstanr",
                       ...)


  if (test_signal) {
    hyp_str = "sd_sample_id__Intercept^2 / (sd_sample_id__Intercept^2 + sigma^2) = 0"
    # hyp = brms::hypothesis(pglmm_fit, hyp_str, class = NULL)
    #
    # outcome_signal = pglmm_fit$fit %>%
    #   as.data.frame() %>%
    #   tibble::as_tibble() %>%
    #   dplyr::mutate(i = 1:n()) %>%
    #   dplyr::mutate(phy_signal = sd_sample_id__Intercept^2 / (sd_sample_id__Intercept^2 + sigma^2)) %>%
    #   dplyr::select(phy_signal)
    #
    # hp = outcome_signal %>%
    #   ggplot(aes(phy_signal)) +
    #   geom_density() +
    #   labs(title = 'Phylogenetic signal posterior') +
    #   theme_light()
    #
    # print(hp)

    pglmm_loo = loo::loo(pglmm_fit)
    base_loo = loo::loo(base_fit)

    message("loo comparison: ")
    comparison = loo::loo_compare(pglmm_loo,
                             base_loo)
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
    pglmm_loo = NULL
    base_loo = NULL
  }

  if (save_object) {
    save(pglmm_fit,
         file = file.path(out_dir, paste0(tree_name, "_pglmm_fit.RData")))
    save(pglmm_fit,
         file = file.path(out_dir, paste0(tree_name, "_base_fit.RData")))
    # V This is what to use once the pglmm_fit is done with cmdstanr
    # pglmm_fit$save_object(file = file.path(out_dir, paste0(tree_name, "_pglmm_fit.RDS")))
  }

  list(pglmm_fit = pglmm_fit,
       base_fit = base_fit,
       loo = list(pglmm_loo = pglmm_loo, base_loo = base_loo, comparison = comparison))

}
