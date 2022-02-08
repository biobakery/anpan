get_bug_name = function(bug_file,
                        remove_pattern = ".genefamilies.tsv") {
  # TODO make this function smarter.
  gsub(remove_pattern, "", basename(bug_file))
}

fit_glms = function(model_input, out_dir, covariates, outcome, bug_name) {

  if (dplyr::n_distinct(model_input[[outcome]]) == 2) {
    # TODO let the user specify a family that overrides this logic
    mod_family = "binomial"
  } else {
    mod_family = "gaussian"
  }

  glm_fits = model_input[,.(data_subset = list(.SD)), by = gene]

  # Progress won't be that hard: https://furrr.futureverse.org/articles/articles/progress.html#package-developers-1
  # p <- progressr::progressor(along = glm_fits$data_subset)
  # ^ That goes right here, then activate the p() calls commented out in fit_glm()
  glm_fits$glm_res = furrr::future_map(.x = glm_fits$data_subset,
                                       .f = safely_fit_glm,
                                       covariates = covariates,
                                       outcome = outcome,
                                       mod_family = mod_family)
  # TODO progress bar with progressr
  # What I have here doesn't work for some reason.


  failed = glm_fits[sapply(glm_fits$glm_res,
                           function(.x) !is.null(.x$error))]
  if (nrow(failed) > 0) {
    # TODO Write out the failures to a warning file with a message
    error_dir = file.path(out_dir, "errors")
    if (!dir.exists(error_dir)) dir.create(error_dir)
    save(failed,
         file = file.path(error_dir, paste0("failures_", bug_name, ".RData")))
  }

  worked = glm_fits[sapply(glm_fits$glm_res,
                           function(.x) is.null(.x$error))]

  # V Find the data.table way to unnest. tidyr is slow.
  all_terms = tidyr::unnest(worked[,.(gene, glm_res = lapply(glm_res, function(.x) .x$result))],
                            cols = c(glm_res))

  readr::write_tsv(all_terms,
                   file = file.path(out_dir, paste0(bug_name, "_all_terms.tsv.gz")))
  # TODO write this and the one two lines down to a "bug_results/" directory rather than the top level output directory

  bug_terms = all_terms %>%
    dplyr::filter(term == "presentTRUE") %>%
    dplyr::arrange(p.value) %>%
    dplyr::mutate(q_bug_wise = p.adjust(p.value, method = 'fdr')) %>%
    dplyr::select(-term)

  readr::write_tsv(bug_terms,
                   file = file.path(out_dir, paste0(bug_name, "_gene_terms.tsv.gz")))

  return(bug_terms)
}

check_prevalence_okay = function(gene_dat, outcome, prevalence_filter) {
  # TODO deprecate this function.
  min_prev_by_outcome = gene_dat %>% dplyr::select(dplyr::all_of(c("present", outcome))) %>%
    table() %>%
    apply(2, function(.x) .x / sum(.x)) %>%
    apply(2, min) # It has to be variable above the prevalence filter in BOTH groups

  if (n_distinct(gene_dat[[outcome]]) != 2 |                             # If only one outcome shows up
      n_distinct(gene_dat$present) == 1 |                         # or if the gene is omni-present or omni-absent
      all(min_prev_by_outcome < prevalence_filter)) {             # or if it doesn't get through the prevalence filter
    return(FALSE) # Then it's not okay
  } else {
    return(TRUE)
  }
}

#' Fit a GLM to one gene
fit_glm = function(gene_dat, covariates, outcome, out_dir,
                   mod_family) {

  # if (!check_prevalence_okay(gene_dat, outcome = outcome, prevalence_filter)) {
  #   return(data.table(term = character(),
  #                     estimate = numeric(),
  #                     std.error = numeric(),
  #                     statistic = numeric(),
  #                     p.value = numeric()))
  # }

  glm_formula = as.formula(paste0(outcome, " ~ ", paste(covariates, collapse = " + "), " + present"))

  res = glm(glm_formula, # TODO adjustable covariates
            data = gene_dat,
            family = 'binomial') %>%
    broom::tidy() %>%
    as.data.table()
  return(res)
}

safely_fit_glm = purrr::safely(fit_glm)

fit_anpan = function(model_input,
                     out_dir,
                     bug_name,
                     covariates,
                     outcome,
                     ...,
                     save_summ = FALSE) {

  main_formula = as.formula(paste0(outcome, " ~ stdCovariates + genes"))
  std_formula = as.formula(paste0("stdCovariates ~ 1 + ", paste(covariates, collapse = " + ")))


  form = brms::bf(main_formula, nl = TRUE) +
    brms::lf(std_formula, center = TRUE) +
    brms::lf(as.formula(paste0("genes ~ 0 + ",
                               paste(names(model_input)[!(names(model_input) %in% c(covariates, outcome, "sample_id"))],
                                     collapse = " + "))), cmc = FALSE)

  p = brms::set_prior("horseshoe(par_ratio = .005)",
                      class = 'b', nlpar = "genes") +
    brms::prior(normal(0,2),
                class = "b", nlpar = 'stdCovariates')

  if (dplyr::n_distinct(model_input[[outcome]]) == 2){
    # TODO allow the user to specify a family that overrides this logic
    mod_family = brms::bernoulli()
  } else {
    mod_family = stats::gaussian()
  }

  model_fit = brms::brm(form,
                        prior = p,
                        data = model_input,
                        family = mod_family,
                        backend = 'cmdstanr',
                        ...)

  summ = summarise_fit(model_fit)

  if (save_summ) {
    save(summ,
         file = file.path(out_dir, paste0(bug_name, ".RData"))) # TODO make the parts of this work
  }

  model_fit

}

#' Run anpan
#'
#' @description Run the anpan model on a single bug
#' @param bug_file path to a gene family file (usually probably from HUMAnN)
#' @param meta_file path to a metadata tsv
#' @param out_dir path to the desired output directory
#' @param model_type either "anpan" or "glm"
#' @param covariates covariates to account for (as a vector of strings)
#' @param discard_absent_samples logical indicating whether to discard samples when a bug is labelled as completely absent
#' @param outcome the name of the outcome variable
#' @param save_filter_stats logical indicating whether to save filter statistics
#' @param ... arguments to pass to brms::brm()
#' @details The specified metadata file must contain columns matching "sample_id"
#'   and the specified covariates and outcome variables.
#' @export
anpan = function(bug_file,
                 meta_file,
                 out_dir,
                 model_type = "anpan",
                 covariates = c("age", "gender"),
                 outcome = "crc",
                 filtering_method = "kmeans",
                 discard_absent_samples = TRUE,
                 annotation_file = NULL,
                 plot_ext = "png",
                 save_filter_stats = TRUE,
                 verbose = TRUE,
                 ...) {

  n_steps = 3
# Checks ------------------------------------------------------------------
  # TODO separate the checks out to a distinct function.
  if (verbose) message(paste0("\n(1/", n_steps, ") Preparing the mise en place (checking inputs)..."))

  if (!(model_type %in% c("anpan", "glm"))) stop('model_type must be either "anpan" or "glm"')

  bug_name = get_bug_name(bug_file)

  if (!dir.exists(out_dir)) {
    if (verbose) message("* Creating output directory.")
    dir.create(out_dir)
  }

  if (save_filter_stats) {
    filter_stats_dir = file.path(out_dir, "filter_stats")
    fs_plot_dir = file.path(filter_stats_dir, 'plots')
    fs_labs_dir = file.path(filter_stats_dir, 'labels')
    if (!dir.exists(filter_stats_dir)) {
      if (verbose) message("* Creating the filter stats directory in the output directory.")
      dir.create(filter_stats_dir)
    }
    if (!dir.exists(fs_plot_dir)) dir.create(fs_plot_dir)
    if (!dir.exists(fs_labs_dir)) dir.create(fs_labs_dir)
  }

  plot_dir = file.path(out_dir, 'plots')
  if (!dir.exists(plot_dir)) {
    if (verbose) message("* Creating output plots directory.")
    dir.create(plot_dir)
  }

  warnings_file = file.path(out_dir, "warnings.txt")
  if (!file.exists(warnings_file)) {
    file.create(warnings_file)
  }

# Filtering ---------------------------------------------------------------

  if (verbose) message(paste0("(2/", n_steps, ") Reading and filtering ", bug_file))
  metadata = read_meta(meta_file,
                       select_cols = c("sample_id", outcome, covariates))

  model_input = read_and_filter(bug_file, metadata = metadata,
                                pivot_wide = model_type == "anpan",
                                covariates = covariates,
                                outcome = outcome,
                                filtering_method = filtering_method,
                                discard_absent_samples = discard_absent_samples,
                                save_filter_stats = save_filter_stats,
                                filter_stats_dir = filter_stats_dir,
                                plot_ext = plot_ext,
                                verbose = verbose)

  if (is.null(model_input)){
    readr::write_lines(paste0(bug_file, " was skipped because no samples passed the filter criteria."),
                       file = warnings_file,
                       append = TRUE)
    return(data.table::data.table())
  }

  if (nrow(model_input) == 0) {
    # ^ if nothing passed the prevalence filter:
    readr::write_lines(paste0(bug_file, " contained no genes that passed the prevalence filter."),
                       file = warnings_file,
                       append = TRUE)

    return(data.table::data.table())
  }

  if (save_filter_stats) {
    if (verbose) message("* Saving filtered data in wide format. ")
    spread_formula = paste(paste(covariates, collapse = " + "), " + sample_id + ", outcome,  " ~ gene",
                           sep = "") %>%
      as.formula()

    wide_dat = dcast(model_input,
                     formula = spread_formula,
                     value.var = 'present')
    readr::write_tsv(wide_dat,
                     file = file.path(filter_stats_dir, paste0("filtered_", bug_name, ".tsv.gz")))
  }


# Fitting -----------------------------------------------------------------

  if (verbose) message(paste0("(3/", n_steps, ") Fitting models to filtered data"))
  res = switch(model_type,
               anpan = fit_anpan(model_input = model_input,
                                 out_dir = out_dir,
                                 covariates = covariates,
                                 outcome = outcome,
                                 bug_name = bug_name,
                                 ...),
               glm   = fit_glms(model_input, out_dir,
                                covariates = covariates,
                                outcome = outcome,
                                bug_name = bug_name))


# Summarizing -------------------------------------------------------------
# Only needed for anpan

# Write output ------------------------------------------------------------
# Done inside fit_glms()

  res$bug_name = bug_name

  return(res)
}

#' Apply anpan to a many bugs
#'
#' @description This function calls anpan() on each gene family file in the
#'   \code{bug_dir} directory and makes a composite data + results plot for
#'   each.
#'
#' @param bug_dir a directory of gene family files
#' @param plot_results logical indicating whether or not to plot the results
#' @param covariates character vector of covariates to include in the model
#' @param discard_absent_samples logical indicating whether to discard samples when a bug is labelled as completely absent
#' @param ... arguments to pass to brms::brm()
#' @inheritParams make_results_plot
#' @inheritParams anpan
#' @export
anpan_batch = function(bug_dir,
                       meta_file,
                       out_dir,
                       model_type = "anpan",
                       covariates = c("age", "gender"),
                       outcome = "crc",
                       filtering_method = "kmeans",
                       discard_absent_samples = TRUE,
                       annotation_file = NULL,
                       save_filter_stats = TRUE,
                       verbose = TRUE,
                       plot_results = TRUE,
                       plot_ext = "png",
                       q_threshold = NULL,
                       n_top = 50,
                       ...) {

  bug_files = get_file_list(bug_dir)
  # anpan is parallelized internally, so just map here.
  all_bug_terms = purrr::map(.x = bug_files,
                             .f = anpan,
                             meta_file = meta_file,
                             out_dir = out_dir,
                             model_type = model_type,
                             filtering_method = filtering_method,
                             discard_absent_samples = discard_absent_samples,
                             covariates = covariates,
                             outcome = outcome,
                             save_filter_stats = save_filter_stats,
                             annotation_file = annotation_file,
                             plot_ext = plot_ext,
                             verbose = verbose,
                             ...) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(q_global = p.adjust(p.value, method = "fdr")) %>%
    dplyr::relocate(bug_name, gene) %>%
    data.table::as.data.table()

  if (model_type == "glm") {
    make_p_value_histogram(all_bug_terms,
                           out_dir = out_dir,
                           plot_ext = plot_ext)
  }

  readr::write_tsv(all_bug_terms,
                   file = file.path(out_dir, 'all_bug_gene_terms.tsv.gz'))

  filter_stats_dir = file.path(out_dir, "filter_stats")
  plot_dir = file.path(out_dir, 'plots')
  if (plot_results) {
    purrr::pmap(all_bug_terms[,.(s = list(.SD)), by = bug_name],
                function(bug_name, s){make_results_plot(res = s,
                                                        bug_name = bug_name,
                                                        covariates = covariates,
                                                        outcome = outcome,
                                                        model_input = fread(file.path(filter_stats_dir, paste0("filtered_", bug_name, ".tsv.gz"))),
                                                        plot_dir = plot_dir,
                                                        annotation_file = annotation_file,
                                                        plot_ext = plot_ext,
                                                        n_top = n_top,
                                                        q_threshold = q_threshold)})
  }

  all_bug_terms

}

#' @export
summarise_anpan = function(model_fit) {
  post_sum = brms::posterior_summary(model_fit)

  post_df = post_sum %>%
    as.data.frame %>%
    tibble::rownames_to_column('param') %>%
    as.data.table()

  post_df
}
