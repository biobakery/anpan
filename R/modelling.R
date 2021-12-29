get_bug_name = function(bug_file,
                        remove_pattern = ".genefamilies.tsv") {
  gsub(remove_pattern, "", basename(bug_file))
}


fit_glms = function(model_input, out_dir, bug_name) {
  glm_fits = model_input[,.(data_subset = list(.SD)), by = gene]

  # Progress won't be that hard: https://furrr.futureverse.org/articles/articles/progress.html#package-developers-1
  # p <- progressr::progressor(steps = dplyr::n_distinct(model_input$gene))
  # ^ That goes right here, then activate the p() calls commented out in fit_glm()
  glm_fits$glm_res = furrr::future_map(glm_fits$data_subset,
                                       function(.x){
                                         # p()
                                         safely_fit_glm(.x)
                                       }) # TODO progress bar with progressr


  failed = glm_fits[sapply(glm_fits$glm_res,
                           \(.x) !is.null(.x$error))]
  if (nrow(failed) > 0){
    # TODO Write out the failures to a warning file with a message
  }

  worked = glm_fits[sapply(glm_fits$glm_res,
                           \(.x) is.null(.x$error))]

  # V Find the data.table way to unnest. tidyr is slow.
  all_terms = tidyr::unnest(worked[,.(gene, glm_res = lapply(glm_res, \(.x) .x$result))],
                            cols = c(glm_res))

  readr::write_tsv(all_terms,
                   file = file.path(out_dir, paste0(bug_name, "_all_terms.tsv")))

  bug_terms = all_terms %>%
    dplyr::filter(term == "presentTRUE") %>%
    dplyr::arrange(p.value) %>%
    dplyr::mutate(q = p.adjust(p.value, method = 'fdr')) %>%
    dplyr::select(-term)

  readr::write_tsv(bug_terms,
                   file = file.path(out_dir, paste0(bug_name, "_gene_terms.tsv")))

  return(bug_terms)
}

check_prevalence_okay = function(gene_dat, prevalence_filter) {

  min_prev_by_outcome = table(gene_dat[,.(present, crc)]) %>%
    apply(2, \(.x) .x / sum(.x)) %>%
    apply(2, min) # It has to be variable above the prevalence filter in BOTH groups

  if (n_distinct(gene_dat$crc) != 2 |                             # If only one outcome shows up
      n_distinct(gene_dat$present) == 1 |                         # or if the gene is omni-present or omni-absent
      all(min_prev_by_outcome < prevalence_filter)) {             # or if it doesn't get through the prevalence filter
    return(FALSE) # Then it's not okay
  } else {
    return(TRUE)
  }
}

#' Fit a GLM to one gene
fit_glm = function(gene_dat, out_dir, prevalence_filter = .05) {

  if (!check_prevalence_okay(gene_dat, prevalence_filter)) {
    return(data.table(term = character(),
                      estimate = numeric(),
                      std.error = numeric(),
                      statistic = numeric(),
                      p.value = numeric()))
  }

  res = glm(crc ~ age + gender + present, # TODO adjustable covariates
            data = gene_dat,
            family = 'binomial') %>%
    broom::tidy() %>%
    as.data.table()
  return(res)
}

safely_fit_glm = purrr::safely(fit_glm)

fit_anpan = function(model_input, bug_name, tpc = 4, ncore = 4, out_path,
                    save_summ = FALSE) {
  form = brms::bf(n_crc ~ stdCovariates + unirefs, nl = TRUE) +
    brms::lf(stdCovariates  ~ 1 + age + gender , center = TRUE) +
    brms::lf(as.formula(paste0("unirefs ~ 0 + ",
                               paste(names(model_input)[-(1:5)],
                                     collapse = " + "))), cmc = FALSE)

  p = brms::set_prior("horseshoe(par_ratio = .005)",
                      class = 'b', nlpar = "unirefs") +
    brms::prior(normal(0,2),
                class = "b", nlpar = 'stdCovariates')

  model_fit = brms::brm(form,
            prior = p,
            data = model_input,
            family = brms::bernoulli(),
            refresh = 50,
            control = list(adapt_delta = .9),
            backend = 'cmdstanr',
            chains = 4,
            cores = ncore,
            threads = brms::threading(tpc))

  summ = summarise_fit(model_fit)

  if (save_summ) {
    save(summ,
         file = file.path(out_path, paste0(bug_name, ".RData"))) # TODO make the parts of this work
  }

  model_fit

}

#' Run anpan
#'
#' @description Run the anpan model on a single bug
#' @param bug_file path to a gene family file (usually probably from HUMAnN)
#' @param meta_file path to a metadata tsv. Must contain the specify covariates
#' @param model_type either "anpan" or "glm"
#' @param covariates covariates to account for CURRENTLY IGNORED
#' @param save_filter_stats logical indicating whether to save filter statistics
#' @export
anpan = function(bug_file,
                 meta_file,
                 out_dir,
                 model_type = "anpan",
                 covariates = c("age", "gender"),
                 filtering_method = "kmeans",
                 annotation_file = NULL,
                 plot_ext = "png",
                 save_filter_stats = TRUE) {


# Checks ------------------------------------------------------------------
  message("Preparing the mise en place (checking inputs)...")

  if (!(model_type %in% c("anpan", "glm"))) stop('model_type must be either "anpan" or "glm"')

  bug_name = get_bug_name(bug_file)
  n_lines = R.utils::countLines(bug_file)

  if (!dir.exists(out_dir)) {
    message("* Creating output directory.")
    dir.create(out_dir)
  }

  if (save_filter_stats) {
    filter_stats_dir = file.path(out_dir, "filter_stats")
    fs_plot_dir = file.path(filter_stats_dir, 'plots')
    if (!dir.exists(filter_stats_dir)) {
      message("* Creating the filter stats directory in the output directory.")
      dir.create(filter_stats_dir)
    }
    if (!dir.exists(fs_plot_dir)) dir.create(fs_plot_dir)
  }

  plot_dir = file.path(out_dir, 'plots')
  if (!dir.exists(plot_dir)) {
    message("* Creating output plots directory.")
    dir.create(plot_dir)
  }


# Filtering ---------------------------------------------------------------

  message("Reading and filtering the data.")
  model_input = read_and_filter(bug_file, read_meta(meta_file),
                                pivot_wide = model_type == "anpan",
                                filtering_method = filtering_method,
                                save_filter_stats = save_filter_stats,
                                filter_stats_dir = filter_stats_dir)


# Fitting -----------------------------------------------------------------

  res = switch(model_type,
               anpan = fit_anpan(model_input, out_dir),
               glm = fit_glms(model_input, out_dir, bug_name = bug_name))

  make_data_plot(res, covariates, model_input, plot_dir = plot_dir,
                 annotation_file = annotation_file,
                 bug_name = bug_name, plot_ext = plot_ext)


# Summarizing -------------------------------------------------------------
# Only needed for anpan

# Write output ------------------------------------------------------------
# Done inside fit_glms()

  return(res)
}

#' @export
anpan_batch = function(bug_dir,
                       meta_file,
                       out_dir,
                       model_type = "anpan",
                       covariates = c("age", "gender"),
                       filtering_method = "kmeans",
                       annotation_file = NULL,
                       plot_ext = "png",
                       save_filter_stats = TRUE) {

  bug_files = get_file_list(bug_dir)
  # anpan is parallelized internally, so just map here.
  all_bug_terms = purrr::map(bug_files,
                             anpan,
                             meta_file = meta_file,
                             out_dir = out_dir,
                             model_type = model_type,
                             filtering_method = filtering_method,
                             covariates = covariates,
                             save_filter_stats = save_filter_stats,
                             annotation_file = annotation_file,
                             plot_ext = plot_ext) %>%
    purrr::imap(\(.x, .y) mutate(.x, bug_name = basename(bug_files)[.y])) %>%
    bind_rows %>%
    mutate(q_global = p.adjust(p.value, method = "fdr"))

  readr::write_tsv(all_bug_terms,
                   file = file.path(out_dir, 'all_bug_gene_terms.tsv'))

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
