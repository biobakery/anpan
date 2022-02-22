get_samp_stats_large = function(gf_file) {
  # If the gene family file is large (e.g. E coli), we have to get the quantiles in chunks
  # TODO Test this function specifically

  nc = readLines(gf_file,
                 n = 1) %>%
    strsplit('\t') %>%
    .[[1]] %>%
    length

  chunks = ceiling(10*((2:nc) / nc))

  for (i in 1:10) {
    gf = fread(gf_file,
               colClasses = list(character = 1, numeric = 2:nc),
               select = c(1, (2:nc)[chunks == i]),
               showProgress = FALSE) %>%
      dplyr::select_all(~gsub("_Abundance-RPKs", "", .))

    names(gf)[1] = "gene"
    gf = gf %>%
      tidyr::separate(gene,
                      into = c("u", "s"),
                      sep = "\\|")

    s_ids = names(gf)[-c(1, 2)]

    gf = gf %>%
      melt(id.vars = c("u", "s"),
           variable.name = 'sample_id',
           value.name = "abd")

    gf$sample_id = factor(gf$sample_id,
                         levels = s_ids)

    if (i == 1){
      chunk_set = gf[, labd := log10(abd)][, .(n_z = sum(abd == 0),
                                               n_nz = sum(abd > 0),
                                               qs = get_qs_and_n_hi_lo(labd[abd > 0]),
                                               quantile = c('q10', 'q50', 'q90', 'n_hi', 'n_lo')), by = sample_id] %>%
        data.table::dcast(sample_id + n_z + n_nz ~ quantile, value.var = 'qs') %>%
        .[, n_diff := (n_lo - n_hi), by = sample_id] %>%
        .[]
    } else {
      chunk_set = rbind(chunk_set,
                        chunk_set = gf[, labd := log10(abd)][, .(n_z = sum(abd == 0),
                                                                 n_nz = sum(abd > 0),
                                                                 qs = get_qs_and_n_hi_lo(labd[abd > 0]),
                                                                 quantile = c('q10', 'q50', 'q90', 'n_hi', 'n_lo')), by = sample_id] %>%
                          data.table::dcast(sample_id + n_z + n_nz ~ quantile, value.var = 'qs') %>%
                          .[, n_diff := (n_lo - n_hi), by = sample_id] %>%
                          .[])
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
  gf[, labd := log10(abd)][, .(n_z = sum(abd == 0),
                               n_nz = sum(abd > 0),
                               qs = get_qs_and_n_hi_lo(labd[abd > 0]),
                               quantile = c('q10', 'q50', 'q90', 'n_hi', 'n_lo')), by = sample_id] %>%
    data.table::dcast(sample_id + n_z + n_nz ~ quantile, value.var = 'qs') %>%
    .[, n_diff := (n_lo - n_hi), by = sample_id] %>%
    .[]
}

filter_with_kmeans = function(gf,
                              samp_stats,
                              save_filter_stats,
                              filter_stats_dir,
                              bug_name,
                              plot_ext = plot_ext) {

  em_input = na.omit(samp_stats[,.(sample_id, n_nz, q50)])
  log_it = dplyr::n_distinct(gf$gene) > 50000
  if (log_it) em_input$n_nz = log1p(em_input$n_nz)
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

  samp_stats$in_right = samp_stats$clust == present_clust
  samp_stats$clust = NULL

  if (save_filter_stats) {
    make_kmeans_dotplot(samp_stats = samp_stats,
                        plot_dir = file.path(filter_stats_dir, "plots"),
                        bug_name,
                        was_logged = log_it,
                        plot_ext = plot_ext)
  }

  filtered_gf = samp_stats[,.(sample_id, in_right)][gf, on = "sample_id"]
  filtered_gf$abd[!filtered_gf$in_right] = 0
  filtered_gf$present = filtered_gf$abd > 0

  filtered_gf
}

#' @export
filter_gf = function(gf,
                     samp_stats,
                     filtering_method = "kmeans",
                     covariates = NULL,
                     outcome = NULL,
                     discard_absent_samples = TRUE,
                     save_filter_stats = FALSE,
                     filter_stats_dir = NULL,
                     plot_ext = "pdf",
                     bug_name = NULL) {

  if (!(filtering_method %in% c("kmeans", "none"))) stop("Specified filtering method not implemented")

  if (filtering_method == 'kmeans') {
    filtered_gf = filter_with_kmeans(gf = gf,
                                     samp_stats = samp_stats,
                                     save_filter_stats = save_filter_stats,
                                     filter_stats_dir = filter_stats_dir,
                                     bug_name = bug_name,
                                     plot_ext = plot_ext)
  } else if (filtering_method == "none") {
    filtered_gf = gf
    filtered_gf$present = filtered_gf$abd > 0
  }

  if (save_filter_stats & filtering_method != "none") {
    make_line_plot(fgf = filtered_gf,
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
                       n_dropped_initial_prevalence_filter = n_start - n_end) # TODO write out to file
  if (!file.exists(initial_prev_filter)) {
    write_tsv_no_progress(drop_df,
                     file = initial_prev_filter)
  } else {
    write_tsv_no_progress(drop_df,
                     file = initial_prev_filter, append = TRUE)
  }

  gf
}

check_table = function(outcome_presence_table,
                       minmax_thresh = 5) {

  if (!all(dim(outcome_presence_table) == c(2, 2))) return(FALSE)

  mins = apply(outcome_presence_table, 2, min)
  maxs = apply(outcome_presence_table, 2, max)
  n = sum(outcome_presence_table)

  ncol(outcome_presence_table) == 2 && all(mins > minmax_thresh) && all(maxs < (n - minmax_thresh))
}

final_prevalence_filter = function(filtered_gf,
                                   outcome,
                                   bn,
                                   minmax_thresh = 5,
                                   filter_stats_dir,
                                   verbose) {

  select_cols = c(outcome, "present", "gene")

  if (dplyr::n_distinct(filtered_gf[[outcome]]) <= 2) {
    to_check = filtered_gf[,..select_cols][,.(sd_table = list(table(.SD))), by = gene]
    to_check$varies_enough = sapply(to_check$sd_table, check_table)
  } else {
    to_check = filtered_gf[,..select_cols][,.(n_pres = sum(present),
                                              n_abs = sum(!present),
                                              .N), by = gene][, varies_enough := (n_pres > minmax_thresh) & (n_abs > minmax_thresh), by = gene][]
  }

  if (any(!to_check$varies_enough)) {
    n_drop = nrow(to_check[!(varies_enough)])
    if (verbose) message(paste0("* Final prevalence filter dropped ", n_drop, " genes."))
    final_filter_file = file.path(filter_stats_dir, "final_prevalence_filter.tsv.gz")
    if (!file.exists(final_filter_file)) {
      write_tsv_no_progress(data.table(bug_name = bn,
                                  n_dropped = n_drop),
                       file = final_filter_file)

    } else {
      write_tsv_no_progress(data.table(bug_name = bn,
                                  n_dropped = n_drop),
                       file = final_filter_file,
                       append = TRUE)
    }
  }

  res = filtered_gf[gene %in% to_check[(varies_enough)]$gene]
  return(res)
}

#' @export
read_and_filter = function(bug_file, metadata, # TODO make metadata optional for this step
                           pivot_wide = TRUE,
                           minmax_thresh = 5, # TODO expose this higher up
                           covariates = NULL,
                           outcome = NULL,
                           filtering_method = "kmeans",
                           discard_absent_samples = TRUE,
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

  gf = read_bug(bug_file, meta = metadata)


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
                          samp_stats = samp_stats,
                          filtering_method = filtering_method,
                          covariates = covariates,
                          outcome = outcome,
                          save_filter_stats = save_filter_stats,
                          filter_stats_dir = filter_stats_dir,
                          plot_ext = plot_ext,
                          bug_name = bug_name)

  if (filtering_method != "none") { # TODO separate these four blocks out to a distinct function
    sample_labels = unique(filtered_gf[,.(sample_id, bug_present = in_right)])
    n_absent = sum(!sample_labels$bug_present)
  }

  if (verbose && filtering_method != "none" && n_absent > 0) {
    message(paste0("* ", n_absent, " samples out of ", nrow(sample_labels), " were determined to not have ", bug_name, " present."))
  }

  if (save_filter_stats && filtering_method != "none") {
    write_tsv_no_progress(sample_labels,
                     file = file.path(filter_stats_dir, 'labels', paste0('sample_labels_', bug_name, '.tsv.gz')))
  }

  if (discard_absent_samples && filtering_method != "none") {
    filtered_gf = filtered_gf[(in_right)]
    if (n_absent != 0 & verbose) message("* Samples with no ", bug_name, " discarded.")
  }

  # third filter: final prevalence filter -----------------------------------

  select_cols = c(outcome, "sample_id")
  filtered_gf = final_prevalence_filter(metadata[, ..select_cols][filtered_gf, on = "sample_id"],
                                        outcome = outcome,
                                        minmax_thresh = minmax_thresh,
                                        filter_stats_dir = filter_stats_dir,
                                        bn = bug_name,
                                        verbose = verbose) %>%
    dplyr::select(-dplyr::all_of(outcome))

  select_cols = c("gene", "present", "sample_id", covariates, outcome)

  if (!is.null(metadata)) {
    joined = filtered_gf[metadata, on = 'sample_id', nomatch = 0] %>%
      dplyr::select(dplyr::all_of(select_cols))
  } else {
    joined = filtered_gf %>%
      dplyr::select(dplyr::all_of(select_cols))
  }

  if (pivot_wide) {
    spread_formula = paste(paste(covariates, collapse = " + "), " + sample_id + ", outcome,  " ~ gene",
                           sep = "") %>%
      as.formula()

    wide_dat = dcast(joined,
                     formula = spread_formula,
                     value.var = 'present')

    return(wide_dat)
  } else {
    return(joined)
  }
}


