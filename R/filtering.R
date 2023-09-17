get_samp_stats_large = function(gf_file) {
  # If the gene family file is large (e.g. E coli), we have to get the quantiles in chunks
  # TODO Test this function specifically

  nc = strsplit(readLines(gf_file,
                          n = 1),
                '\t')[[1]] |>
    length()

  chunks = ceiling(10*((2:nc) / nc))

  for (i in 1:10) {
    gf = fread(gf_file,
               colClasses = list(character = 1, numeric = 2:nc),
               select = c(1, (2:nc)[chunks == i]),
               showProgress = FALSE,
               header = TRUE) |>
      dplyr::select_all(~gsub("_Abundance-RPKs", "", .))

    names(gf)[1] = "gene"
    col_classes = gf |> sapply(class)

    if (any(col_classes[-1] != "numeric")) {
      for (i in (which(col_classes[-1] != "numeric") + 1)) {
        gf[[i]] = as.numeric(gf[[i]])
      }
    }

    gf = gf[, c("u", "s") := tstrsplit(gene, split = "\\|")][,-"gene"] |>
      dplyr::relocate(u, s)

    s_ids = names(gf)[-c(1, 2)]

    gf = gf |>
      melt(id.vars = c("u", "s"),
           variable.name = 'sample_id',
           value.name = "abd")

    gf$sample_id = factor(gf$sample_id,
                         levels = s_ids)

    if (i == 1){
      dcast_input = gf[, labd := log10(abd)][, .(n_z = sum(abd == 0),
                                                 n_nz = sum(abd > 0),
                                                 qs = get_qs_and_n_hi_lo(labd[abd > 0]),
                                                 quantile = c('q10', 'q50', 'q90', 'n_hi', 'n_lo')), by = sample_id]

      chunk_set = data.table::dcast(data = dcast_input,
                                    formula = sample_id + n_z + n_nz ~ quantile, value.var = 'qs')[, n_diff := (n_lo - n_hi), by = sample_id][]

    } else {
      dcast_input = gf[, labd := log10(abd)][, .(n_z = sum(abd == 0),
                                                 n_nz = sum(abd > 0),
                                                 qs = get_qs_and_n_hi_lo(labd[abd > 0]),
                                                 quantile = c('q10', 'q50', 'q90', 'n_hi', 'n_lo')), by = sample_id]

      dcast_res = data.table::dcast(data = dcast_input,
                                    formula = sample_id + n_z + n_nz ~ quantile,
                                    value.var = 'qs')[, n_diff := (n_lo - n_hi), by = sample_id][]

      chunk_set = rbind(chunk_set,
                        dcast_res)
    }
  }

  chunk_set
}

get_qs_and_n_hi_lo = function(.x){
  qs = quantile(.x,
                probs = c(.1, .5, .9))

  n_hi = sum(.x > qs[3])
  n_lo = sum(.x > qs[1])
  c(qs, n_hi, n_lo)
}

get_samp_stats = function(gf){
  dcast_input = gf[, labd := log10(abd)][, .(n_z = sum(abd == 0),
                                             n_nz = sum(abd > 0),
                                             qs = get_qs_and_n_hi_lo(labd[abd > 0]),
                                             quantile = c('q10', 'q50', 'q90', 'n_hi', 'n_lo')), by = sample_id]

  res = data.table::dcast(data = dcast_input,
                          formula = sample_id + n_z + n_nz ~ quantile,
                          value.var = 'qs')[, n_diff := (n_lo - n_hi), by = sample_id][]
  return(res)
}

get_genomes_stats = function(genomes_file) {
  genomes_df = fread(genomes_file,
                     header = TRUE,
                     showProgress = FALSE)

  gene_counts = genomes_df[,-1][,lapply(.SD, function(.x) sum((.x + 0) != 0))] |>
    unlist()
  # ^ The "+ 0" is to convert logicals to numeric if necessary

  n_genomes = ncol(genomes_df) - 1

  if (n_genomes >= 5) {
    lt = mean(gene_counts) - 2*sd(gene_counts)
    lt_source = "2 SDs below mean genome size"
  } else {
    lt = mean(gene_counts)*2/3
    lt_source = "2/3 mean genome size"
  }

  data.table(n_genomes = n_genomes,
             mean_genes = mean(gene_counts),
             lower_threshold = lt,
             lt_source = lt_source)
}

filter_with_kmeans = function(gf,
                              samp_stats,
                              genomes_stats = NULL,
                              save_filter_stats,
                              filter_stats_dir,
                              bug_name,
                              plot_ext = plot_ext) {

  em_input = na.omit(samp_stats[,.(sample_id, n_nz, q50)])
  log_it = dplyr::n_distinct(gf$gene) > 50000
  if (log_it) em_input$n_nz = log1p(em_input$n_nz)

  if (dplyr::n_distinct(em_input$n_nz) == 1) {
    # This means that all samples either had ALL genes present or were ALL zero. Revert to 1D kmeans
    # in that case. This might happen if there are samples that are very deeply sequenced.

    message("* Samples passed to kmeans had either ALL genes present or NO genes present. Reverting to 1D clustering on median log abundance.")

    scrunch = 2 # Scrunch the y-axis of the plots to make sure k-means doesn't accidentally produce a horizontal decision boundary. (i.e cut off the top of the U shape)
    # TODO make sure this doesn't mess up anything downstream
    em_input$q50 = scale(em_input$q50) / scrunch

    km_res = kmeans(as.matrix(em_input$q50),
                    centers = 2)

    samp_stats$in_right = NA
    samp_stats$clust = NA
    present_clust = which.max(km_res$centers[,1])
    absent_clust = which.min(km_res$centers[,1])

    samp_stats$clust[!is.na(samp_stats$q50)] = km_res$cluster
    samp_stats$clust[is.na(samp_stats$q50)] = absent_clust

    samp_stats$in_right = samp_stats$clust == present_clust
    samp_stats$clust = NULL

  } else {
    em_input$n_nz = scale(em_input$n_nz)

    scrunch = 2 # Scrunch the y-axis of the plots to make sure k-means doesn't accidentally produce a horizontal decision boundary. (i.e cut off the top of the U shape)
    # TODO make sure this doesn't mess up anything downstream
    em_input$q50 = scale(em_input$q50) / scrunch

    km_res = kmeans(as.matrix(em_input[,-1]),
                    centers = matrix(c(-1,1,1,-1), nrow = 2))

    samp_stats$in_right = NA
    samp_stats$clust = NA

    present_clust = which.max(km_res$centers[,1])
    absent_clust = which.min(km_res$centers[,1])

    samp_stats$clust[!is.na(samp_stats$q50)] = km_res$cluster
    samp_stats$clust[is.na(samp_stats$q50)] = absent_clust

    if (!is.null(genomes_stats)) {
      n_flipped = samp_stats[n_nz < genomes_stats$lower_threshold & clust == present_clust] |> nrow()
      genomes_stats$flipped = n_flipped
      samp_stats[n_nz < genomes_stats$lower_threshold, clust := absent_clust]
    }

    samp_stats$in_right = samp_stats$clust == present_clust
    samp_stats$clust = NULL

    if (save_filter_stats) {
      plot_kmeans_dots(samp_stats = samp_stats,
                       plot_dir = file.path(filter_stats_dir, "plots"),
                       bug_name,
                       was_logged = log_it,
                       plot_ext = plot_ext,
                       genomes_stats = genomes_stats)
    }
  }

  filtered_gf = samp_stats[,.(sample_id, in_right)][gf, on = "sample_id"]
  filtered_gf$abd[!filtered_gf$in_right] = 0
  filtered_gf$present = filtered_gf$abd > 0

  if (!is.null(genomes_stats)) {
    # Set those below the threshold to "absent"
    # TODO
  }

  filtered_gf
}

#' Filter a gene family file
#' @description This is the function that anpan uses to filter an input
#'   genefamily data.table. This can be useful if you want to play around with
#'   different filtering options without repeatedly re-reading or checking the
#'   file.
#'
#' @param gf a gene family data.table
#' @param samp_stats a data.table of sample statistics
#' @inheritParams anpan
#' @inheritParams read_and_filter
#' @export
filter_gf = function(gf,
                     samp_stats,
                     filtering_method = "kmeans",
                     covariates = NULL,
                     outcome = NULL,
                     genomes_file = NULL,
                     discard_poorly_covered_samples = TRUE,
                     save_filter_stats = FALSE,
                     filter_stats_dir = NULL,
                     plot_ext = "pdf",
                     bug_name = NULL) {

  if (!(filtering_method %in% c("kmeans", "none"))) stop("Specified filtering method not implemented")

  if (filtering_method == 'kmeans') {

    if (!is.null(genomes_file)) {
      genomes_stats = get_genomes_stats(genomes_file)
    } else {
      genomes_stats = NULL
    }

    filtered_gf = filter_with_kmeans(gf = gf,
                                     samp_stats = samp_stats,
                                     genomes_stats = genomes_stats,
                                     save_filter_stats = save_filter_stats,
                                     filter_stats_dir = filter_stats_dir,
                                     bug_name = bug_name,
                                     plot_ext = plot_ext)
  } else if (filtering_method == "none") {
    filtered_gf = gf
    filtered_gf$present = filtered_gf$abd > 0
  }

  if (save_filter_stats & filtering_method != "none") {
    plot_lines(fgf = filtered_gf,
               bug_name = bug_name,
               covariates = covariates,
               outcome = outcome,
               omit_na = TRUE,
               plot_dir = file.path(filter_stats_dir, "plots"),
               plot_ext = plot_ext)
  }

  return(filtered_gf)
}

initial_prevalence_filter = function(gf,
                                     meta, outcome,
                                     bug_name,
                                     minmax_thresh,
                                     filter_stats_dir,
                                     verbose) {

  if (dplyr::n_distinct(meta[[outcome]]) == 2) {

    gf = gf[,.(sample_id, abd,
               varies_enough = sum(abd != 0) < (.N - minmax_thresh) & sum(abd != 0) > minmax_thresh),
            by = c("gene")]
  } else {
    gf = gf[,.(sample_id, abd,
               varies_enough = sum(abd != 0) < (.N - minmax_thresh) & sum(abd != 0) > minmax_thresh),
            by = gene]
  }

  if (!any(gf$varies_enough)) {
    return(data.table::data.table())
  }

  n_start = dplyr::n_distinct(gf$gene)
  gf = gf[(varies_enough)][,.(gene, sample_id, abd)]
  n_end = dplyr::n_distinct(gf$gene)

  if (verbose) message(paste0("* Initial prevalence filter dropped ", n_start - n_end, " genes out of ", n_start, " present in the input file."))
  initial_prev_filter = file.path(filter_stats_dir, "initial_prevalence_filter.tsv.gz")
  drop_df = data.table(bug = bug_name,
                       n_dropped_initial_prevalence_filter = n_start - n_end)

  if (!is.null(filter_stats_dir) && !file.exists(initial_prev_filter)) {
    write_tsv_no_progress(drop_df,
                     file = initial_prev_filter)
  } else if (!is.null(filter_stats_dir)) {
    write_tsv_no_progress(drop_df,
                     file = initial_prev_filter, append = TRUE)
  }

  gf
}

check_table = function(outcome_presence_table,
                       minmax_thresh) {

  if (!all(dim(outcome_presence_table) == c(2, 2))) return(FALSE)

  mins = apply(outcome_presence_table, 2, min)
  maxs = apply(outcome_presence_table, 2, max)
  n = sum(outcome_presence_table)

  ncol(outcome_presence_table) == 2 && all(mins > minmax_thresh) && all(maxs < (n - minmax_thresh))
}

final_prevalence_filter = function(filtered_gf,
                                   outcome,
                                   bn, # bug_name
                                   minmax_thresh,
                                   filter_stats_dir,
                                   verbose) {

  select_cols = c(outcome, "present", "gene")

  if (dplyr::n_distinct(filtered_gf[[outcome]]) <= 2) {
    to_check = filtered_gf[,..select_cols][,.(sd_table = list(table(.SD))), by = gene]
    to_check$varies_enough = sapply(to_check$sd_table, check_table,
                                    minmax_thresh = minmax_thresh)
  } else {
    to_check = filtered_gf[,..select_cols][,.(n_pres = sum(present),
                                              n_abs = sum(!present),
                                              .N), by = gene][, varies_enough := (n_pres > minmax_thresh) & (n_abs > minmax_thresh), by = gene][]
  }

  if (any(!to_check$varies_enough)) {
    n_drop = nrow(to_check[!(varies_enough)])

    if (verbose) message(paste0("* Final prevalence filter dropped ", n_drop, " genes."))

    final_filter_file = file.path(filter_stats_dir, "final_prevalence_filter.tsv.gz")

    if (!is.null(filter_stats_dir) && !file.exists(final_filter_file)) {
      write_tsv_no_progress(data.table(bug_name = bn,
                                       n_dropped = n_drop),
                            file = final_filter_file)

    } else if (!is.null(filter_stats_dir)){
      write_tsv_no_progress(data.table(bug_name = bn,
                                       n_dropped = n_drop),
                            file = final_filter_file,
                            append = TRUE)
    }
  }

  res = filtered_gf[gene %in% to_check[(varies_enough)]$gene]
  return(res)
}

#' Read and filter a gene family file
#'
#' @inheritParams anpan
#' @param minmax_thresh genes must have at least this many (or N - this many)
#'   non-zero observations or else be discarded. NULL defaults to \code{floor(0.005*nrow(metadata))}.
#' @param discard_poorly_covered_samples logical indicating whether to discard samples where the genes of a bug are poorly covered
#' @param pivot_wide logical indicating whether to return data in wide format
#' @param filtering_method either "kmeans" or "none"
#' @param filter_stats_dir directory to save filtering statistics to
#' @param plot_ext extension to use for plots
#' @export
read_and_filter = function(bug_file, metadata, # TODO make metadata optional for this step
                           pivot_wide = TRUE,
                           minmax_thresh,
                           covariates = NULL,
                           outcome = NULL,
                           genomes_file = NULL,
                           filtering_method = "kmeans",
                           discretize_inputs = TRUE,
                           discard_poorly_covered_samples = TRUE,
                           save_filter_stats = TRUE,
                           filter_stats_dir = NULL,
                           plot_ext = "pdf",
                           verbose = TRUE) {

  n_lines = R.utils::countLines(bug_file)
  is_large = n_lines > 200000

  if (is_large) {
    warning("This gene family file is huge. It would be better to use read_and_filter_large() here.")
    # TODO make the downstream functions have a chunked = TRUE argument
    return(NULL)
  }

  if (save_filter_stats & is.null(filter_stats_dir)) {
    stop("You said save the filter statistics but didn't provide filter_stats_dir")
  }

  bug_name = get_bug_name(bug_file)

  if (verbose) message(paste0("* Reading " , bug_file))

  gf = read_bug(bug_file, meta = metadata) # gf = gene families

  if (!("sample_id" %in% names(gf))) {
    warning("No samples matching metadata samples found in gene file")
    return(NULL)
  }

  # first filter: the initial prevalence filter -----------------------------

  gf = initial_prevalence_filter(gf,
                                 meta = metadata,
                                 outcome = outcome,
                                 bug_name = bug_name,
                                 minmax_thresh = minmax_thresh,
                                 filter_stats_dir = filter_stats_dir,
                                 verbose = verbose)
  if (nrow(gf) == 0) {
    return(NULL)
  }

  if (dplyr::n_distinct(gf$gene) == 1) filtering_method = "none" # TODO write out to warnings.txt if this happens

  if (filtering_method != "none" && !is_large) {
    if (verbose) message("* Computing sample statistics...")
    samp_stats = get_samp_stats(gf)
  } else {
    samp_stats = NA
  }

  # second filter: the sample filter ----------------------------------------
  # This filter uses sample statistics to drop samples where the bug is entirely
  # absent, not just lacking specific genes.

  filtered_gf = filter_gf(gf,
                          samp_stats        = samp_stats,
                          filtering_method  = filtering_method,
                          covariates        = covariates,
                          outcome           = outcome,
                          genomes_file      = genomes_file,
                          save_filter_stats = save_filter_stats,
                          filter_stats_dir  = filter_stats_dir,
                          plot_ext          = plot_ext,
                          bug_name          = bug_name)

  if (filtering_method != "none") { # TODO separate these four blocks out to a distinct function
    sample_labels = unique(filtered_gf[,.(sample_id, bug_well_covered = in_right)])
    n_poorly_covered = sum(!sample_labels$bug_well_covered)
  }

  if (verbose && filtering_method != "none" && n_poorly_covered > 0) {
    message(paste0("* ", n_poorly_covered, " samples out of ", nrow(sample_labels), " had poor coverage of the genes of ", bug_name))
  }

  if (save_filter_stats && filtering_method != "none") {
    write_tsv_no_progress(sample_labels,
                          file = file.path(filter_stats_dir, 'labels', paste0('sample_labels_', bug_name, '.tsv.gz')))
  }

  if (discard_poorly_covered_samples && filtering_method != "none") {
    filtered_gf = filtered_gf[(in_right)]
    if (n_poorly_covered != 0 & verbose) message("* Samples with ", bug_name, " poorly covered have been discarded.")
  }

  # third filter: final prevalence filter -----------------------------------

  select_cols = c(outcome, "sample_id")
  filtered_gf = final_prevalence_filter(metadata[, ..select_cols][filtered_gf, on = "sample_id"],
                                        outcome = outcome,
                                        minmax_thresh = minmax_thresh,
                                        filter_stats_dir = filter_stats_dir,
                                        bn = bug_name,
                                        verbose = verbose) |>
    dplyr::select(-dplyr::all_of(outcome))

  if (!discretize_inputs) {
    bug_covariate = "abd"
  } else {
    bug_covariate = "present"
  }

  select_cols = c("gene", bug_covariate, "sample_id", covariates, outcome)

  if (nrow(filtered_gf) == 0){
    return(NULL)
  }

  if (!is.null(metadata)) {
    joined = filtered_gf[metadata, on = 'sample_id', nomatch = 0] |>
      dplyr::select(dplyr::all_of(select_cols))
  } else {
    joined = filtered_gf |>
      dplyr::select(dplyr::all_of(select_cols))
  }

  if (pivot_wide) {
    if (!is.null(covariates)) {
      cov_str = paste(covariates, collapse = " + ")
      form_lhs = paste(cov_str, "sample_id", outcome, sep = " + ")
    } else {
      form_lhs = paste("sample_id", outcome, sep = " + ")
    }

    spread_formula = paste(form_lhs,
                           " ~ gene",
                           sep = "") |>
      as.formula()

    wide_dat = dcast(joined,
                     formula = spread_formula,
                     value.var = bug_covariate)

    return(wide_dat)
  } else {
    return(joined)
  }
}

check_filter_res = function(model_input,
                            bug_file,
                            warnings_file,
                            covariates, outcome,
                            already_wide,
                            bug_covariate,
                            filter_stats_dir) {

  if (is.null(model_input)) {
    cat(paste0(bug_file, " was skipped because no samples passed the filter criteria."),
        file = warnings_file,
        append = TRUE,
        sep = "\n")

    return(NULL)
  }

  if (nrow(model_input) == 0) {
    # ^ if nothing passed the prevalence or kmeans filters:
    cat(paste0(bug_file, " contained no genes that the prevalence filter."),
        file = warnings_file,
        append = TRUE,
        sep = "\n")

    return(NULL)
  }

  if (already_wide) {
    wide_dat = model_input
  } else {
    if (!is.null(covariates)) {
      cov_str = paste(covariates, collapse = " + ")
      form_lhs = paste(cov_str, "sample_id", outcome, sep = " + ")
    } else {
      form_lhs = paste("sample_id", outcome, sep = " + ")
    }

    spread_formula = paste(form_lhs,
                           " ~ gene",
                           sep = "") |>
      as.formula()

    wide_dat = dcast(model_input,
                     formula = spread_formula,
                     value.var = bug_covariate)
  }

  bug_name = get_bug_name(bug_file)

  write_tsv_no_progress(wide_dat,
                        file = file.path(filter_stats_dir, paste0("filtered_", bug_name, ".tsv.gz")))

  return(NULL)
}

safely_read_and_filter = purrr::safely(read_and_filter)

#' Filter a batch of files
#'
#' @description This function applies anpan::read_and_filter() to a set of files.
#' @return A list of filtered data frames
#' @inheritParams read_and_filter
#' @inheritParams anpan_batch
#' @export
filter_batch = function(bug_dir, meta_file,
                        filter_stats_dir,
                        pivot_wide = TRUE,
                        minmax_thresh = NULL,
                        covariates = NULL,
                        outcome = NULL,
                        filtering_method = "kmeans",
                        discretize_inputs = TRUE,
                        discard_poorly_covered_samples = TRUE,
                        omit_na = FALSE,
                        plot_ext = "pdf",
                        verbose = TRUE) {

  if (!dir.exists(filter_stats_dir)) {
    if (verbose) message("* Creating the filter stats directory.")
    dir.create(filter_stats_dir)
  }

  fs_plot_dir = file.path(filter_stats_dir, 'plots')
  fs_labs_dir = file.path(filter_stats_dir, 'labels')

  if (!dir.exists(fs_plot_dir)) dir.create(fs_plot_dir)
  if (!dir.exists(fs_labs_dir)) dir.create(fs_labs_dir)

  warnings_file = file.path(filter_stats_dir, "warnings.txt")
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

  bug_files = get_file_list(bug_dir)

  p = progressr::progressor(along = bug_files)

  filter_list = furrr::future_map(bug_files,
                                  ~{filter_res = safely_read_and_filter(.x,
                                                                 metadata = metadata,
                                                                 pivot_wide = pivot_wide,
                                                                 covariates = covariates,
                                                                 outcome = outcome,
                                                                 filtering_method = filtering_method,
                                                                 discretize_inputs = discretize_inputs,
                                                                 minmax_thresh = minmax_thresh,
                                                                 discard_poorly_covered_samples = discard_poorly_covered_samples,
                                                                 save_filter_stats = TRUE,
                                                                 filter_stats_dir = filter_stats_dir,
                                                                 plot_ext = plot_ext,
                                                                 verbose = verbose)
                                    if (!is.null(filter_res$error)) return(NULL)
                                    check_filter_res(filter_res$result,
                                                     bug_file = .x,
                                                     warnings_file = warnings_file,
                                                     covariates = covariates,
                                                     outcome = outcome,
                                                     already_wide = pivot_wide,
                                                     bug_covariate = bug_covariate,
                                                     filter_stats_dir = filter_stats_dir)
                                    p()
                                    return(filter_res$result)
                                  })

  return(filter_list)

}


