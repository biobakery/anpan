get_bug_name = function(bug_file,
                        remove_pattern = ".genefamilies.tsv|.genefamilies.tsv.gz") {
  # TODO make this function smarter.
  gsub(remove_pattern, "", basename(bug_file))
}

fit_glms = function(model_input, out_dir, covariates, outcome, bug_name,
                    discretized_inputs = TRUE,
                    glm_fun,
                    fastglm_method = 1) {

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
                                       .f = glm_fun,
                                       covariates = covariates,
                                       outcome = outcome,
                                       mod_family = mod_family,
                                       discretized_inputs = discretized_inputs,
                                       fastglm_method = fastglm_method,
                                       .options = furrr::furrr_options(seed = 123))
  # TODO progress bar with progressr
  # What I have here doesn't work for some reason.

  failed = glm_fits[sapply(glm_fits$glm_res,
                           function(.x) !is.null(.x$error))]
  #TODO ^ detect failures better. Sometimes glm or fastglm "fails" but just returns NAs instead of erroring.

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

  write_tsv_no_progress(all_terms,
                   file = file.path(out_dir, paste0(bug_name, "_all_terms.tsv.gz")))
  # TODO write this and the one two lines down to a "bug_results/" directory rather than the top level output directory

  if (!discretized_inputs) {
    bug_term_name = "abd"
  } else {
    bug_term_name = "presentTRUE"
  }

  bug_terms = all_terms %>%
    dplyr::filter(term == bug_term_name) %>%
    dplyr::arrange(p.value) %>%
    dplyr::mutate(q_bug_wise = p.adjust(p.value, method = 'fdr')) %>%
    dplyr::select(-term)

  write_tsv_no_progress(bug_terms,
                   file = file.path(out_dir, paste0(bug_name, "_gene_terms.tsv.gz")))

  return(bug_terms)
}

#' Fit a GLM to one gene
fit_glm = function(gene_dat, covariates, outcome, out_dir,
                   discretized_inputs = TRUE,
                   mod_family,
                   fastglm_method = NULL) {
  # fastglm_method is NOT used in this function. It's an argument here to make
  # calling a varying function easier inside fit_glms().

  if (!discretized_inputs) {
    bug_covariate = "abd"
  } else {
    bug_covariate = "present"
  }

  glm_formula = as.formula(paste0(outcome, " ~ ", paste(covariates, collapse = " + "), " + ", bug_covariate))

  res = glm(glm_formula,
            data = gene_dat,
            family = mod_family) %>%
    broom::tidy() %>% # Only place I use broom, maybe do it manually
    as.data.table()

  return(res)
}

fit_fastglm = function(gene_dat, covariates, outcome, out_dir,
                       discretized_inputs = TRUE,
                       mod_family, fastglm_method = 1) {

  y = gene_dat[[outcome]]

  if (dplyr::n_distinct(y) == 2) y = 1*y

  if (!discretized_inputs) {
    bug_covariate = "abd"
  } else {
    bug_covariate = "present"
  }

  glm_formula = as.formula(paste0(" ~ ", paste(covariates, collapse = " + "), " + ", bug_covariate))
  x = model.matrix(glm_formula,
                   data = gene_dat)

  res = fastglm::fastglm(x = x, y = y,
                         family = mod_family,
                         method = fastglm_method) %>% summary %>%
    .[['coefficients']] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    data.table::as.data.table()

  names(res) = c("term", "estimate", "std.error", "statistic", "p.value")

  return(res)
}

safely_fit_glm = purrr::safely(fit_glm)
safely_fit_fastglm = purrr::safely(fit_fastglm)

clean_ushoe_summary = function(ushoe_fit,
                               cov_names,
                               gene_names) {

  res = ushoe_fit$summary() %>%
    as.data.table()

  res$variable["b_covariates_Intercept" == res$variable] = "intercept"
  res$variable[grepl("b_covariate", res$variable)] = cov_names[-1]

  res = res[!grepl("zb_genes|Intercept_covariates|hs_local", variable)]
  res$index = as.numeric(gsub("[A-Za-z[:punct:]]+", "", res$variable))

  res$gene = gene_names[res$index]

  res = res %>%
    dplyr::select(param = variable, gene, mean:ess_tail)

  return(res)
}

fit_horseshoe = function(model_input,
                         out_dir,
                         bug_name,
                         covariates,
                         outcome,
                         save_fit = TRUE,
                         skip_large = TRUE,
                         ...) {

  if (dplyr::n_distinct(model_input[[outcome]]) == 2) {
    # TODO allow the user to specify a family that overrides this logic
    model_path = system.file("stan", "logistic_ushoe.stan",
                             package = 'anpan',
                             mustWork = TRUE)
    ushoe_model = cmdstanr::cmdstan_model(stan_file = model_path, quiet = TRUE)
  } else {
    model_path = system.file("stan", "gaussian_ushoe.stan",
                             package = 'anpan',
                             mustWork = TRUE)
    ushoe_model = cmdstanr::cmdstan_model(stan_file = model_path, quiet = TRUE)
  }

  if (skip_large && ncol(model_input) > (10002 + length(covariates))) {
    warnings_file = file.path(out_dir, "warnings.txt")
    readr::write_lines(paste0(bug_name, " was skipped because there are over ten thousand genes after filtering. Add skip_large = FALSE to disable this behavior."),
                       file = warnings_file,
                       append = TRUE)
    return(NULL)
  }

  cov_formula = as.formula(paste0("~ 1 + ", paste(covariates, collapse = " + ")))

  X_genes = model_input %>%
    dplyr::select(-dplyr::all_of(c('sample_id', outcome, covariates))) %>%
    as.matrix()
  X_genes = 0+X_genes # convert to 0/1
  X_covariates = model.matrix(cov_formula, data = model_input)

  data_list = list(N = nrow(model_input),
                   Y = as.numeric(model_input[[outcome]]),
                   K_covariates = 1 + length(covariates), # + 1 for intercept
                   X_covariates = X_covariates,
                   K_genes = ncol(model_input) - (2 + length(covariates)),  # sample_id + outcome + length(covariates)
                   X_genes = X_genes,
                   hs_df_genes = 1,
                   hs_df_global_genes = 1,
                   hs_df_slab_genes = 4,
                   hs_scale_global_genes = .01 / sqrt(nrow(model_input)),
                   hs_scale_slab_genes = 3,
                   prior_only = 0)

  ushoe_fit = ushoe_model$sample(data = data_list,
                                 iter_sampling = 1000, # TODO make these options user-accessible
                                 iter_warmup = 500,
                                 parallel_chains = 4,
                                 adapt_delta = .9,
                                 refresh = 0)

  res = clean_ushoe_summary(ushoe_fit,
                            colnames(X_covariates),
                            colnames(X_genes))

  if (save_fit) {
    fit_dir = file.path(out_dir, "fits")
    if (!dir.exists(fit_dir)) dir.create(fit_dir)
    ushoe_fit$save_object(file = file.path(fit_dir, paste0(bug_name, "_fit.RDS")))
  }

  return(res)

}

#' Run anpan
#'
#' @description Run the anpan model on a single bug
#' @param bug_file path to a gene family file (usually probably from HUMAnN)
#' @param meta_file path to a metadata tsv
#' @param out_dir path to the desired output directory
#' @param prefiltered_dir an optional directory to pre-filtered data from an
#'   earlier run to skip the filtering step
#' @param model_type either "horseshoe", "glm", or "fastglm"
#' @param outcome the name of the outcome variable
#' @param covariates covariates to account for (as a vector of strings)
#' @param skip_large logical indicating whether to skip bugs with over 5k genes.
#'   Only used when model_type = "horseshoe".
#' @param save_fit logical indicating whether to save horseshoe fit objects.
#'   Only used when model_type = "horseshoe".
#' @param discard_absent_samples logical indicating whether to discard samples
#'   when a bug is labelled as completely absent
#' @param omit_na logical indicating whether to omit incomplete cases
#' @param filtering_method method to use for filtering samples. Either "kmeans"
#'   or "none"
#' @param discretize_inputs logical indicating whether to discretize the input
#'   abundance measurements (0/nonzero --> FALSE/TRUE) before passing them to
#'   the modelling function
#' @param save_filter_stats logical indicating whether to save filter statistics
#' @param ... arguments to pass to [cmdstanr::sample()] if applicable
#' @details The specified metadata file must contain columns matching
#'   "sample_id" and the specified covariates and outcome variables.
#' @export
anpan = function(bug_file,
                 meta_file,
                 out_dir,
                 prefiltered_dir = NULL,
                 model_type = "fastglm",
                 covariates = c("age", "gender"),
                 outcome = "crc",
                 omit_na = FALSE,
                 filtering_method = "kmeans",
                 discretize_inputs = TRUE,
                 skip_large = TRUE,
                 save_fit = TRUE,
                 discard_absent_samples = TRUE,
                 plot_ext = "png",
                 save_filter_stats = TRUE,
                 verbose = TRUE,
                 ...) {

  n_steps = 3
# Checks ------------------------------------------------------------------
  # TODO separate the checks out to a distinct function.
  if (verbose) message(paste0("\n(1/", n_steps, ") Preparing the mise en place (checking inputs)..."))

  if (!(model_type %in% c("horseshoe", "glm", "fastglm"))) stop('model_type must be either "horseshoe", "glm", or "fastglm"')

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

  metadata = read_meta(meta_file,
                       select_cols = c("sample_id", outcome, covariates),
                       omit_na = omit_na)

# Filtering ---------------------------------------------------------------

  if (!discretize_inputs) {
    bug_covariate = "abd"
  } else {
    bug_covariate = "present"
  }

  if (!is.null(prefiltered_dir)) {

    if (verbose) message(paste0("(2/", n_steps, ") Reading ", bug_file, " from the provided directory of pre-filtered data."))

    prefiltered_bug = list.files(prefiltered_dir, full.names = TRUE, pattern = bug_name)

    if (length(prefiltered_bug) == 0) {
      readr::write_lines(paste0(bug_file, " was skipped because no matching file was found in the pre-filtered data."),
                         file = warnings_file,
                         append = TRUE)
      if (verbose) message(paste0("(3/", n_steps, ") No matching file found in in pre-filtered data directory - Model fitting skipped"))
      return(data.table::data.table())
    }

    model_input = fread(prefiltered_bug)

    if (model_type %in% c("glm", "fastglm")) {
      model_input = data.table::melt(model_input,
                                     id.vars = c(covariates, outcome, "sample_id"),
                                     variable.name = "gene",
                                     value.name = bug_covariate)
    }

  } else {
    if (verbose) message(paste0("(2/", n_steps, ") Reading and filtering ", bug_file))

    model_input = read_and_filter(bug_file, metadata = metadata,
                                  pivot_wide = model_type == "horseshoe",
                                  covariates = covariates,
                                  outcome = outcome,
                                  filtering_method = filtering_method,
                                  discretize_inputs = discretize_inputs,
                                  discard_absent_samples = discard_absent_samples,
                                  save_filter_stats = save_filter_stats,
                                  filter_stats_dir = filter_stats_dir,
                                  plot_ext = plot_ext,
                                  verbose = verbose)

    if (is.null(model_input)) {
      readr::write_lines(paste0(bug_file, " was skipped because no samples passed the filter criteria."),
                         file = warnings_file,
                         append = TRUE)
      if (verbose) message(paste0("(3/", n_steps, ") Nothing passed filters - Model fitting skipped"))
      return(data.table::data.table())
    }

    if (nrow(model_input) == 0) {
      # ^ if nothing passed the prevalence or kmeans filters:
      readr::write_lines(paste0(bug_file, " contained no genes that the prevalence filter."),
                         file = warnings_file,
                         append = TRUE)
      if (verbose) message(paste0("(3/", n_steps, ") Nothing passed filters - Model fitting skipped"))
      return(data.table::data.table())
    }

    if (save_filter_stats) {
      if (verbose) message("* Saving filtered data in wide format. ")

      spread_formula = paste(paste(covariates, collapse = " + "), " + sample_id + ", outcome,  " ~ gene",
                             sep = "") %>%
        as.formula()

      wide_dat = dcast(model_input,
                       formula = spread_formula,
                       value.var = bug_covariate)

      write_tsv_no_progress(wide_dat,
                            file = file.path(filter_stats_dir, paste0("filtered_", bug_name, ".tsv.gz")))
    }
  }


# Fitting -----------------------------------------------------------------

  if (verbose) message(paste0("(3/", n_steps, ") Fitting models to filtered data"))
  res = switch(model_type,
               horseshoe = fit_horseshoe(model_input = model_input,
                                         out_dir = out_dir,
                                         covariates = covariates,
                                         outcome = outcome,
                                         bug_name = bug_name,
                                         skip_large = skip_large,
                                         save_fit = save_fit,
                                         ...),
               glm   =   fit_glms(model_input, out_dir,
                                  covariates = covariates,
                                  outcome = outcome,
                                  bug_name = bug_name,
                                  discretized_inputs = discretize_inputs,
                                  glm_fun = safely_fit_glm),
               fastglm = fit_glms(model_input, out_dir,
                                  covariates = covariates,
                                  outcome = outcome,
                                  bug_name = bug_name,
                                  discretized_inputs = discretize_inputs,
                                  glm_fun = safely_fit_fastglm))


# Summarizing -------------------------------------------------------------
# Only needed for anpan

# Write output ------------------------------------------------------------
# Done inside fit_glms()

  res$bug_name = bug_name

  return(res)
}

safely_anpan = purrr::safely(anpan)

#' Apply anpan to a many bugs
#'
#' @description This function calls anpan() on each gene family file in the
#'   \code{bug_dir} directory and makes a composite data + results plot for
#'   each.
#'
#' @param bug_dir a directory of gene family files
#' @param plot_result logical indicating whether or not to plot the results
#' @param covariates character vector of covariates to include in the model
#' @param prefiltered_dir an optional directory to pre-filtered data from an
#'   earlier run to skip the filtering step
#' @param discard_absent_samples logical indicating whether to discard samples
#'   when a bug is labelled as completely absent
#' @param ... arguments to pass to [cmdstanr::sample()] if applicable
#' @details \code{bug_dir} should be a directory of gene (or SNVs or pathways)
#'   abundance files, one for each bug.
#' @inheritParams plot_results
#' @inheritParams anpan
#' @export
anpan_batch = function(bug_dir,
                       meta_file,
                       out_dir,
                       prefiltered_dir = NULL,
                       model_type = "fastglm",
                       covariates = c("age", "gender"),
                       outcome = "crc",
                       omit_na = FALSE,
                       filtering_method = "kmeans",
                       discretize_inputs = TRUE,
                       discard_absent_samples = TRUE,
                       skip_large = TRUE,
                       save_fit = TRUE,
                       annotation_file = NULL,
                       save_filter_stats = TRUE,
                       verbose = TRUE,
                       plot_result = TRUE,
                       plot_ext = "png",
                       q_threshold = NULL,
                       n_top = 50,
                       width = 10,
                       height = 8,
                       ...) {

  call = match.call()

  fn_call_string = paste0(gsub(', (?!")',
                               ",\n            ",
                               as.character(enquote(call))[2],
                               perl = TRUE),
                          "\n")

  if (verbose & !interactive()) message(paste0("Now running:\n\n", fn_call_string))

  if (!dir.exists(out_dir)) {
    if (verbose) message("* Creating output directory.")
    dir.create(out_dir)
  }

  readr::write_lines(fn_call_string,
                     file = file.path(out_dir, "anpan_batch_call.txt"))

  bug_files = get_file_list(bug_dir)

  p = progressr::progressor(along = bug_files)

  anpan_results = purrr::map(.x = bug_files,
                             .f = function(.x) {
                               anpan_res = safely_anpan(.x,
                                                 meta_file = meta_file,
                                                 out_dir = out_dir,
                                                 prefiltered_dir = prefiltered_dir,
                                                 model_type = model_type,
                                                 skip_large = skip_large,
                                                 save_fit = save_fit,
                                                 filtering_method = filtering_method,
                                                 discretize_inputs = discretize_inputs,
                                                 discard_absent_samples = discard_absent_samples,
                                                 covariates = covariates,
                                                 outcome = outcome,
                                                 omit_na = omit_na,
                                                 save_filter_stats = save_filter_stats,
                                                 annotation_file = annotation_file,
                                                 plot_ext = plot_ext,
                                                 verbose = verbose,
                                                 ...)
                               p()
                               return(anpan_res)
                             }) |>
    purrr::transpose() |>
    as_tibble()

  worked = anpan_results |>
    dplyr::filter(purrr::map_lgl(error, is.null))

  errors = anpan_results |>
    dplyr::mutate(bug_file = basename(bug_files)) |>
    dplyr::filter(purrr::map_lgl(result, is.null)) |>
    dplyr::relocate(bug_file) |>
    dplyr::select(-result)

  if (nrow(errors) > 0) {
    save(errors,
         file = file.path(out_dir, 'errors.RData'))
  }

  all_bug_terms = worked$result %>%
    dplyr::bind_rows() %>%
    dplyr::relocate(bug_name, gene) %>%
    data.table::as.data.table()

  if (model_type %in% c("glm", "fastglm")) {
    all_bug_terms$q_global = p.adjust(all_bug_terms$p.value, method = "fdr")

    plot_p_value_histogram(all_bug_terms,
                           out_dir = out_dir,
                           plot_ext = plot_ext)
  }

  if (!is.null(annotation_file)) {
    anno = fread(annotation_file)

    all_bug_terms = anno[all_bug_terms, on = 'gene'] %>%
      dplyr::relocate(bug_name, gene) %>%
      dplyr::relocate(annotation, .after = dplyr::last_col())
  }

  write_tsv_no_progress(all_bug_terms,
                        file = file.path(out_dir, 'all_bug_gene_terms.tsv.gz'))

  filter_stats_dir = file.path(out_dir, "filter_stats")
  plot_dir = file.path(out_dir, 'plots')

  if (verbose) message("Saving results plots to output directory...")
  if (plot_result) {
    plotting_input = all_bug_terms[,.(s = list(.SD)), by = bug_name]
    p = progressr::progressor(steps = nrow(plotting_input))
    plot_list = furrr::future_pmap(plotting_input,
                                   function(bug_name, s){plot_res = safely_plot_results(res = s,
                                                                                        bug_name = bug_name,
                                                                                        covariates = covariates,
                                                                                        outcome = outcome,
                                                                                        model_input = fread(file.path(filter_stats_dir, paste0("filtered_", bug_name, ".tsv.gz"))),
                                                                                        discretize_inputs = discretize_inputs,
                                                                                        plot_dir = plot_dir,
                                                                                        annotation_file = annotation_file,
                                                                                        plot_ext = plot_ext,
                                                                                        n_top = n_top,
                                                                                        q_threshold = q_threshold,
                                                                                        cluster = 'both',
                                                                                        show_trees = TRUE,
                                                                                        width = width,
                                                                                        height = height)
                                                         p()
                                                         return(plot_res)})
  }

  return(all_bug_terms)

}
