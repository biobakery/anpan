get_bug_name = function(bug_file,
                        remove_pattern = ".genefamilies.tsv$|.genefamilies.tsv.gz$|.tsv$|.tsv.gz$") {
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
  # ^ That goes right here, then activate the p() calls commented out in fit_fastglm()
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
  #TODO ^ detect failures better. Sometimes fastglm "fails" but just returns NAs instead of erroring.

  if (nrow(failed) > 0) {
    # TODO Write out the failures to a warning file with a message
    error_dir = file.path(out_dir, "errors")
    if (!dir.exists(error_dir)) dir.create(error_dir)
    save(failed,
         file = file.path(error_dir, paste0("failures_", bug_name, ".RData")))
  }

  worked = glm_fits[sapply(glm_fits$glm_res,
                           function(.x) is.null(.x$error))]

  unnest_input = worked[,.(gene, glm_res = lapply(glm_res, function(.x) .x$result))]
  all_terms = unnest_input[,data.table::rbindlist(glm_res), by = gene] |>
    tibble::as_tibble()

  # TODO add a check here if ^ is empty and exit gracefully. Necessary for anpan().

  write_tsv_no_progress(all_terms,
                   file = file.path(out_dir, paste0(bug_name, "_all_terms.tsv.gz")))
  # TODO write this and the one two lines down to a "bug_results/" directory rather than the top level output directory

  if (!discretized_inputs) {
    bug_term_name = "abd"
  } else {
    bug_term_name = "presentTRUE"
  }

  bug_terms = all_terms |>
    dplyr::filter(term == bug_term_name) |>
    dplyr::arrange(p.value) |>
    dplyr::mutate(q_bug_wise = p.adjust(p.value, method = 'fdr')) |>
    dplyr::select(-term) |>
    data.table::as.data.table()

  write_tsv_no_progress(bug_terms,
                   file = file.path(out_dir, paste0(bug_name, "_gene_terms.tsv.gz")))

  return(bug_terms)
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

  res = summary(fastglm::fastglm(x = x, y = y,
                                 family = mod_family,
                                 method = fastglm_method))[['coefficients']] |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    data.table::as.data.table()

  names(res) = c("term", "estimate", "std.error", "statistic", "p.value")

  return(res)
}

safely_fit_fastglm = purrr::safely(fit_fastglm)

clean_ushoe_summary = function(ushoe_fit,
                               cov_names,
                               gene_names) {

  res = ushoe_fit$summary() |>
    as.data.table()

  res$variable["b_covariates_Intercept" == res$variable] = "intercept"
  res$variable[grepl("b_covariate", res$variable)] = cov_names[-1]

  res = res[!grepl("zb_genes|Intercept_covariates|hs_local", variable)]
  res$index = as.numeric(gsub("[A-Za-z[:punct:]]+", "", res$variable))

  res$gene = gene_names[res$index]

  res = res |>
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

    cat(paste0(bug_name, " was skipped because there are over ten thousand genes after filtering. Add skip_large = FALSE to disable this behavior."),
        file = warnings_file,
        append = TRUE,
        sep = "\n")

    return(NULL)
  }

  cov_formula = as.formula(paste0("~ 1 + ", paste(covariates, collapse = " + ")))

  X_genes = model_input |>
    dplyr::select(-dplyr::all_of(c('sample_id', outcome, covariates))) |>
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
                   hs_scale_slab_genes = 1,
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
#' @description Run the anpan gene model on a single bug
#' @param bug_file path to a gene family file (usually probably from HUMAnN)
#' @param meta_file path to a metadata tsv
#' @param out_dir path to the desired output directory
#' @param genomes_file optional file giving gene presence/absence of representative isolate genomes
#' @param prefiltered_dir an optional directory to pre-filtered data from an earlier run to skip the
#'   filtering step
#' @param model_type either "horseshoe" or "fastglm"
#' @param outcome the name of the outcome variable
#' @param covariates covariates to account for (as a vector of strings)
#' @param skip_large logical indicating whether to skip bugs with over 5k genes. Only used when
#'   model_type = "horseshoe".
#' @param save_fit logical indicating whether to save horseshoe fit objects. Only used when
#'   model_type = "horseshoe".
#' @param discard_poorly_covered_samples logical indicating whether to discard samples where the genes of a bug are poorly covered
#' @param omit_na logical indicating whether to omit incomplete cases of the metadata
#' @param filtering_method method to use for filtering samples. Either "kmeans" or "none"
#' @param discretize_inputs logical indicating whether to discretize the input abundance
#'   measurements (0/nonzero --> FALSE/TRUE) before passing them to the modelling function
#' @param save_filter_stats logical indicating whether to save filter statistics
#' @param ... arguments to pass to [cmdstanr::sample()] if applicable
#' @details The specified metadata file must contain columns matching "sample_id" and the specified
#'   covariates and outcome variables.
#'
#'   If provided, \code{genomes_file} is used to refine the filtering process. The format must be
#'   genes as rows, with the first column giving the gene id (usually a UniRef90 identifier), and
#'   subsequent columns representing isolate genomes. The entries of the isolate genome columns
#'   should give 0/1 indicators of whether or not the gene is present in the isolate. The gene
#'   counts present in these isolates are used to establish the typical number of genes present in a
#'   strain of the species and a lower threshold on the number of acceptable gene observations. If
#'   >=5 isolate genomes are available, the lower threshold is 2 standard deviations below the mean,
#'   otherwise it is 2/3 of the mean.
#' @returns a data.table of model statistics for each gene
#' @seealso [anpan_batch()]
#' @export
anpan = function(bug_file,
                 meta_file,
                 out_dir,
                 genomes_file = NULL,
                 prefiltered_dir = NULL,
                 model_type = "fastglm",
                 covariates = c("age", "gender"),
                 outcome = "crc",
                 omit_na = FALSE,
                 filtering_method = "kmeans",
                 discretize_inputs = TRUE,
                 skip_large = TRUE,
                 save_fit = TRUE,
                 discard_poorly_covered_samples = TRUE,
                 plot_ext = "png",
                 save_filter_stats = TRUE,
                 verbose = TRUE,
                 ...) {

  n_steps = 3
# Checks ------------------------------------------------------------------
  # TODO separate the checks out to a distinct function.
  if (verbose) message(paste0("\n(1/", n_steps, ") Preparing the mise en place (checking inputs)..."))

  if (!(model_type %in% c("horseshoe", "fastglm"))) stop('model_type must be either "horseshoe" or "fastglm"')

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

  if (!(is.numeric(metadata[[outcome]]) || is.logical(metadata[[outcome]]))) {
    error_msg = paste0("The specified outcome variable in the metadata is neither numeric nor logical. The specified outcome variable is class: ",
                       class(metadata[[outcome]])[1])
    stop(error_msg)
  }

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
      cat(paste0(bug_file, " was skipped because no matching file was found in the pre-filtered data."),
          file = warnings_file,
          append = TRUE,
          sep = "\n")

      if (verbose) message(paste0("(3/", n_steps, ") No matching file found in in pre-filtered data directory - Model fitting skipped"))
      return(data.table::data.table())
    }

    model_input = fread(prefiltered_bug,
                        header = TRUE)

    if (model_type %in% c("fastglm")) {
      model_input = data.table::melt(model_input,
                                     id.vars = c(covariates, outcome, "sample_id"),
                                     variable.name = "gene",
                                     value.name = bug_covariate)
    }

  } else {
    if (verbose) message(paste0("(2/", n_steps, ") Reading and filtering ", bug_file))

    model_input = read_and_filter(bug_file,
                                  metadata               = metadata,
                                  pivot_wide             = model_type == "horseshoe",
                                  covariates             = covariates,
                                  outcome                = outcome,
                                  genomes_file           = genomes_file,
                                  filtering_method       = filtering_method,
                                  discretize_inputs      = discretize_inputs,
                                  discard_poorly_covered_samples = discard_poorly_covered_samples,
                                  save_filter_stats      = save_filter_stats,
                                  filter_stats_dir       = filter_stats_dir,
                                  plot_ext               = plot_ext,
                                  verbose                = verbose)

    if (is.null(model_input)) {
      cat(paste0(bug_file, " was skipped because no samples passed the filter criteria."),
          file = warnings_file,
          append = TRUE,
          sep = "\n")

      if (verbose) message(paste0("(3/", n_steps, ") Nothing passed filters - Model fitting skipped"))
      return(data.table::data.table())
    }

    if (nrow(model_input) == 0) {
      # ^ if nothing passed the prevalence or kmeans filters:
      cat(paste0(bug_file, " contained no genes that passed the prevalence filter."),
          file = warnings_file,
          append = TRUE,
          sep = "\n")

      if (verbose) message(paste0("(3/", n_steps, ") Nothing passed filters - Model fitting skipped"))
      return(data.table::data.table())
    }

    if (save_filter_stats) {
      if (verbose) message("* Saving filtered data in wide format. ")

      spread_formula = paste(paste(covariates, collapse = " + "), " + sample_id + ", outcome,  " ~ gene",
                             sep = "") |>
        as.formula()

      if (model_type == "horseshoe") {
        wide_dat = model_input
      } else {
        wide_dat = dcast(model_input,
                         formula = spread_formula,
                         value.var = bug_covariate)
      }

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
#' @param genomes_dir an optional directory of genome files
#' @param plot_result logical indicating whether or not to plot the results
#' @param covariates character vector of covariates to include in the model
#' @param prefiltered_dir an optional directory to pre-filtered data from an
#'   earlier run to skip the filtering step
#' @param discard_poorly_covered_samples logical indicating whether to discard samples where the genes of a bug are poorly covered
#' @param annotation_file a path to a file giving annotations for each gene
#' @param ... arguments to pass to [cmdstanr::sample()] if applicable
#' @details \code{bug_dir} should be a directory of gene (or SNV or pathway)
#'   abundance files, one for each bug.
#'
#'   \code{annotation} file must have two columns named "gene" and "annotation"
#'
#'   See \code{?anpan()} for the format / usage if providing genome files.
#' @inheritParams plot_results
#' @inheritParams anpan
#' @returns a data.table of model statistics for each bug:gene combination
#' @seealso [anpan()]
#' @export
anpan_batch = function(bug_dir,
                       meta_file,
                       out_dir,
                       genomes_dir = NULL,
                       prefiltered_dir = NULL,
                       model_type = "fastglm",
                       covariates = c("age", "gender"),
                       outcome = "crc",
                       omit_na = FALSE,
                       filtering_method = "kmeans",
                       discretize_inputs = TRUE,
                       discard_poorly_covered_samples = TRUE,
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

  if (!is.null(annotation_file)) {
    anno = fread(annotation_file, nrows = 3, header = TRUE)
    if (!all(c("gene", "annotation") %in% names(anno))) {
      stop("Couldn't find the gene and annotation columns in the supplied annotation file.")
    }
  }

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

  cat(fn_call_string,
      file = file.path(out_dir, "anpan_batch_call.txt"),
      sep = "\n")

  bug_files = get_file_list(bug_dir)

  if (!is.null(genomes_dir)) {
    genomes_files = get_file_list(genomes_dir)
  }

  metadata = read_meta(meta_file,
                       select_cols = c("sample_id", outcome, covariates),
                       omit_na = omit_na)

  if (!(is.numeric(metadata[[outcome]]) || is.logical(metadata[[outcome]]))) {
    error_msg = paste0("The specified outcome variable in the metadata is neither numeric nor logical. The specified outcome variable is class: ",
                       class(meta[[outcome]])[1])
    stop(error_msg)
  }

  p = progressr::progressor(along = bug_files)

  anpan_results = purrr::map(.x = bug_files,
                             .f = function(.x) {
                               bn = get_bug_name(.x)
                               if (!is.null(genomes_dir)) {
                                 genomes_file = grep(bn, genomes_files, value = TRUE)
                                 if (length(genomes_file) == 0) {
                                   genomes_file = NULL
                                   warning(paste0("No genome file found for ", bn, " in the provided genome directory."))
                                 }
                                 if (length(genomes_file) > 1) {
                                   genomes_file = genomes_file[1]
                                   warning(paste0("Multiple genome files found for ", bn, " , using the first to establish typical genome size"))
                                 }
                               } else {
                                 genomes_file = NULL
                               }
                               anpan_res = safely_anpan(.x,
                                                 meta_file = meta_file,
                                                 out_dir = out_dir,
                                                 genomes_file = genomes_file,
                                                 prefiltered_dir = prefiltered_dir,
                                                 model_type = model_type,
                                                 skip_large = skip_large,
                                                 save_fit = save_fit,
                                                 filtering_method = filtering_method,
                                                 discretize_inputs = discretize_inputs,
                                                 discard_poorly_covered_samples = discard_poorly_covered_samples,
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

  if (nrow(worked) > 0) {
    all_bug_terms = worked$result |>
      dplyr::bind_rows() |>
      dplyr::relocate(bug_name, gene) |>
      data.table::as.data.table()
  } else {
    stop("No models fit successfully, see errors.RData")
  }

  if (model_type %in% c("fastglm")) {
    all_bug_terms$q_global = p.adjust(all_bug_terms$p.value, method = "fdr")

    plot_p_value_histogram(all_bug_terms,
                           out_dir = out_dir,
                           plot_ext = plot_ext)
  }

  if (!is.null(annotation_file)) {
    anno = fread(annotation_file, header = TRUE)

    if (!("gene" %in% colnames(anno))) {
      warning('No "gene" column in annotation file. Annotations not joined onto result')

    } else if (any(duplicated(anno$gene))) {
      warning("Gene annotations are not unique. Annotations not joined onto result")
    } else {
      all_bug_terms = anno[all_bug_terms, on = 'gene'] |>
        dplyr::relocate(bug_name, gene) |>
        dplyr::relocate(annotation, .after = dplyr::last_col())
    }
  }

  write_tsv_no_progress(all_bug_terms,
                        file = file.path(out_dir, 'all_bug_gene_terms.tsv.gz'))

  filter_stats_dir = file.path(out_dir, "filter_stats")
  plot_dir = file.path(out_dir, 'plots')

  if (verbose) message("Saving results plots to output directory...")
  if (plot_result) {
    plotting_input = all_bug_terms[,.(s = list(.SD)), by = bug_name]
    p = progressr::progressor(steps = nrow(plotting_input))

    filtered_data_dir = if (!is.null(prefiltered_dir)) prefiltered_dir else filter_stats_dir

    filtered_file_list = list.files(filtered_data_dir, full.names = TRUE)

    if (length(covariates) > 2) {
      warning("Only using the first two covariates for annotation color bars on plots.")
      covariates = covariates[1:2]
    }

    plot_list = furrr::future_pmap(plotting_input,
                                   function(bug_name, s){plot_res = safely_plot_results(res = s,
                                                                                        bug_name = bug_name,
                                                                                        covariates = covariates,
                                                                                        outcome = outcome,
                                                                                        model_input = fread(grep(filtered_file_list,
                                                                                                                 pattern = bug_name,
                                                                                                                 value = TRUE),
                                                                                                            showProgress = FALSE,
                                                                                                            header = TRUE),
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


  # Check if there are any bugs with a lot of hits, if so, issue a warning.
  if (model_type == "fastglm") {
    qt = ifelse(is.null(q_threshold), .1, q_threshold)
    hit_prop_df = all_bug_terms[, .(prop_bug    = mean(q_bug_wise < qt),
                                    prop_global = mean(q_global < qt)),
                                by = bug_name][order(-prop_bug, -prop_global)]

    if (any(hit_prop_df$prop_global > .01) || any(hit_prop_df$prop_bug > .01)) {
      warning_msg = paste0("The bug(s) listed below (and in warnings.txt) had more than 1% of their genes exhibit an association with an FDR q-value below ",
                           qt,
                           " (either bug-wise or globally). This may indicate that there is some phylogenetic structure within the species confounded with the outcome. This might show up visually as large blocks on the results plot. You can try evaluating a PGLMM with anpan_pglmm() (or anpan_pglmm_batch()) to quantify the phylogenetic signal. If you don't have a phylogeny for the species, see the \"Estimating phylogenies from gene matrices\" section of the vignette.\n") |>
        strwrap(initial = '\n    ') |>
        paste(collapse = "\n") |>
        paste(paste(capture.output(hit_prop_df[prop_global > .01 | prop_bug > .01]),
                    collapse = "\n"),
              sep = "\n\n")

      warning(warning_msg)

      warnings_file = file.path(out_dir, "warnings.txt")

      cat(warning_msg,
          file = warnings_file,
          append = TRUE,
          sep = "\n")
    }
  }

  return(all_bug_terms)
}

summarise_metadata_variable = function(meta_var) {
  if (is.numeric(meta_var)) {
    return(mean(meta_var))
  } else if (is.character(meta_var) || is.factor(meta_var) || is.logical(meta_var)) {
    freq_table = sort(table(meta_var), decreasing = TRUE)
    res = names(freq_table[1])
    if (res %in% c("TRUE", "FALSE")) res = as.logical(res)
    if (is.factor(meta_var)) res = factor(res, levels = levels(meta_var))
    return(res)
  }
}

aggregate_by_subject = function(filtered_sample_file,
                                subject_dir,
                                subject_sample_map,
                                covariates,
                                outcome) {
  sample_df = filtered_sample_file |>
    fread(header = TRUE) |>
    melt(id.vars = c(covariates, outcome, "sample_id"),
         variable.name = "gene", value.name = "present")

  subject_sample_map[, bug_well_covered := sample_id %in% sample_df$sample_id]
  subject_sample_map[, n_samples := .N, by = subject_id]

  joined = merge(subject_sample_map,
        sample_df, by = "sample_id",
        all = TRUE)

  joined$present[!joined$bug_well_covered] = FALSE

  prop_df = joined[, .(present_prop = sum(bug_well_covered & present) / n_samples[1]), by = c("subject_id", "gene")][!is.na(gene)] |>
    dcast(subject_id ~ gene, value.var = "present_prop")

  select_cols = c(covariates, outcome, "sample_id")
  output_cols = c(covariates, outcome, "subject_id")

  subject_sample_map |>
    left_join(unique(sample_df[, ..select_cols]), by = "sample_id")

  other_cols = dplyr::inner_join(subject_sample_map, unique(sample_df[, ..select_cols]), by = "sample_id")[,..output_cols] |>
    unique() |>
    dplyr::group_by(subject_id) |>
    dplyr::summarise(dplyr::across(.cols = dplyr::all_of(c(covariates, outcome)),
                                   summarise_metadata_variable)) |>
    as.data.table()

  res = other_cols[prop_df, on = "subject_id"] |>
    dplyr::rename(sample_id = subject_id) # Needed to work with anpan()

  fwrite(res,
         file = file.path(subject_dir,
                          basename(filtered_sample_file)),
         sep = "\t")

  return(NULL)
}

read_filter_write = function(.x,
                             metadata,
                             covariates,
                             outcome,
                             filtering_method,
                             sample_wise_filter_stats_dir,
                             plot_ext = "pdf") {

  read_res = read_and_filter(.x,
                             metadata         = metadata,
                             covariates       = covariates,
                             outcome          = outcome,
                             filtering_method = "kmeans",
                             filter_stats_dir = sample_wise_filter_stats_dir,
                             plot_ext         = "pdf")

  if (is.null(read_res)) return(NULL)

  fwrite(read_res,
         file = file.path(sample_wise_filter_stats_dir,
                          paste0("filtered_", basename(.x))))
}

#' Use repeated measures to refine the gene model
#' @param subject_sample_map a data frame between the mapping between subject_id and sample_id
#' @details This function performs the standard anpan filtering on all samples, then uses the
#'   subject-sample map to compute the proportion of samples with the bug. This gives a gene
#'   _proportion_ matrix (instead of a presence/absence matrix) which is then passed to
#'   \code{anpan_batch(filtering_method = "none", discretize_inputs = FALSE)}.
#'
#'   In cases where subject metadata varies by sample, the mean is taken if the variable is numeric,
#'   otherwise it is tabulated and the most frequent category is selected as the subject-level
#'   metadata value. This tabulation will respect factor ordering if you'd like to alter the value
#'   selected in the event of ties.
#' @inheritParams anpan_batch
#' @export
anpan_repeated_measures = function(subject_sample_map,
                                   bug_dir,
                                   meta_file,
                                   out_dir,
                                   model_type = "fastglm",
                                   covariates = c("age", "gender"),
                                   outcome = "crc",
                                   omit_na = FALSE,
                                   filtering_method = "kmeans",
                                   discard_poorly_covered_samples = TRUE,
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

  cat(fn_call_string,
      file = file.path(out_dir, "anpan_batch_call.txt"),
      sep = "\n")

  if (is.character(subject_sample_map) && file.exists(subject_sample_map)) {
    subject_sample_map = fread(subject_sample_map, header = TRUE)
  } else if (!is.data.frame(subject_sample_map)) {
    stop("Couldn't read subject_sample_map from a file nor is it a data frame.")
  }

  if (!is.data.table(subject_sample_map)) subject_sample_map = as.data.table(subject_sample_map)

  if (!(all(c("subject_id", "sample_id") %in% names(subject_sample_map)))) {
    stop("Couldn't find the subject_id and sample_id columns in the subject_sample_map")
  }

  if (!is.null(annotation_file)) {
    anno = fread(annotation_file, nrows = 3, header = TRUE)
    if (!all(c("gene", "annotation") %in% names(anno))) {
      stop("Couldn't find the gene and annotation columns in the supplied annotation file.")
    }
  }

  bug_files = get_file_list(bug_dir)

  metadata = read_meta(meta_file,
                       select_cols = c("sample_id", outcome, covariates),
                       omit_na = omit_na)

  sample_wise_filter_stats_dir = file.path(out_dir, "sample_wise_filter_stats_dir")
  dir.create(sample_wise_filter_stats_dir)
  dir.create(file.path(sample_wise_filter_stats_dir, "plots"))
  dir.create(file.path(sample_wise_filter_stats_dir, "labels"))

  # TODO Add a progressor here
  bug_files |>
    purrr::walk(read_filter_write,
                metadata = metadata,
                covariates = covariates,
                outcome = outcome,
                filtering_method = "kmeans",
                sample_wise_filter_stats_dir = sample_wise_filter_stats_dir,
                plot_ext = "pdf")

  subject_dir = file.path(out_dir, "subject_dir")
  dir.create(subject_dir)

  # TODO Add a progressor here
  list.files(sample_wise_filter_stats_dir,
             full.names = TRUE,
             pattern = "filtered") |>
    lapply(aggregate_by_subject,
           subject_dir = subject_dir,
           subject_sample_map = subject_sample_map,
           covariates = covariates,
           outcome = outcome)

  # For each subject, multiply the proportion of samples that have the bug by the proportion of
  # samples that have the gene to get the final gene data. Pass that to anpan_batch with no
  # filtering and no discretizing.

  subj_metadata = read_meta(meta_file,
                                       select_cols = c("sample_id", outcome, covariates),
                                       omit_na = omit_na) |>
    dplyr::left_join(subject_sample_map, by = "sample_id") |>
    dplyr::select(-sample_id) |>
    dplyr::select(dplyr::all_of(c("subject_id", outcome, covariates))) |>
    unique() |>
    dplyr::rename(sample_id = subject_id)

  anpan_batch(bug_dir = subject_dir,
              meta_file = subj_metadata,
              model_type = model_type,
              omit_na = omit_na,
              skip_large = skip_large,
              discretize_inputs = FALSE,
              filtering_method = "none",
              prefiltered_dir = subject_dir,
              covariates = covariates,
              outcome = outcome,
              out_dir = file.path(out_dir, "model_output"))
}
