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
               select = c(1, (2:nc)[chunks == i])) %>%
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

fit_mixture = function(samp_stats) {
  em_input = na.omit(samp_stats[,.(sample_id, n_z, q50)])
  em_input$n_z = scale(em_input$n_z)
  em_input$q50 = scale(em_input$q50)

  mf = mixtools::mvnormalmixEM(as.matrix(em_input[,-1]),
                               lambda = c(.75, .25),
                               mu = list(c(-1, 1), 1, -1))
  return(list(res = mf, input = em_input))
}

get_component_densities = function(mix_fit, plot = FALSE) {
  input_mat = as.matrix(mix_fit$input[,-1])

  left_i = which.min(sapply(mix_fit$res$mu, function(.x) .x[1]))
  right_i = setdiff(1:2, left_i)

  samp_df = mix_fit$input
  samp_df$left_dens = mixtools::dmvnorm(input_mat,
                                        mu = mix_fit$res$mu[[left_i]],
                                        sigma = mix_fit$res$sigma[[left_i]])
  samp_df$right_dens = mixtools::dmvnorm(input_mat,
                                        mu = mix_fit$res$mu[[right_i]],
                                        sigma = mix_fit$res$sigma[[right_i]])

  samp_df[,in_right := right_dens > left_dens]

  if (plot) {
    samp_df %>%
      ggplot(aes(n_z, q50)) +
      geom_point(aes(color = in_right)) +
      labs(title = 'EM based decision for A. putredinis') +
      theme_light()
  }

  samp_df[,.(sample_id, in_right, left_dens, right_dens)]
}

# TODO implement this
# Check if the components overlap mostly overlap
# Turns out this is pretty hard to do:
# https://math.stackexchange.com/questions/1114879/detect-if-two-ellipses-intersect
check_mix_fit = function(mix_fit) {
  # Check they don't overlap too much
  # Check that one is to the upper left of the other
  TRUE
}

filter_with_mixture = function(gf,
                               samp_stats,
                               save_filter_stats,
                               filter_stats_dir,
                               bug_name){
  mix_fit = fit_mixture(samp_stats)

  mix_checks = check_mix_fit(mix_fit)

  if (!mix_checks) {
    stop("Mixture fitting seems to have failed")
  }

  if (save_filter_stats) {

    make_hex_plot(samp_stats = samp_stats,
                  mix_fit = mix_fit,
                  plot_dir = file.path(filter_stats_dir, "plots"),
                  bug_name = bug_name)

    mix_fit_res = mix_fit$res

    save(samp_stats, mix_fit_res,
         file = file.path(filter_stats_dir, paste0(bug_name, ".RData")))
  }

  mix_labels = get_component_densities(mix_fit)
  label_df = mix_labels[samp_stats, on = "sample_id"] # TODO write out this table
  label_df$in_right[is.na(label_df$in_right)] = TRUE # all zero samples don't have the species

  filtered_gf = label_df[,.(sample_id, in_right)][gf, on = "sample_id"]

  filtered_gf$abd[filtered_gf$in_right] = 0
  filtered_gf$present = filtered_gf$abd > 0

  filtered_gf
}

filter_with_kmeans = function(gf,
                              samp_stats,
                              save_filter_stats,
                              filter_stats_dir,
                              bug_name) {
  em_input = na.omit(samp_stats[,.(sample_id, n_z, q50)])
  if (dplyr::n_distinct(gf$gene) > 70000) em_input$n_z = log1p(em_input$n_z)
  em_input$n_z = scale(em_input$n_z)

  scrunch = 2 # Scrunch the y-axis of the plots to make sure k-means doesn't accidentally produce a horizontal decision boundary. (i.e cut off the top of the U shape)
  # TODO make sure this doesn't mess up anything downstream
  em_input$q50 = scale(em_input$q50) / scrunch

  km_res = kmeans(as.matrix(em_input[,-1]),
                  centers = matrix(c(-1,1,1,-1), nrow = 2))

  samp_stats$in_right = NA
  samp_stats$clust = NA
  low_clust = which.max(km_res$centers[,1])

  samp_stats$clust[!is.na(samp_stats$q50)] = km_res$cluster
  samp_stats$clust[is.na(samp_stats$q50)] = low_clust

  samp_stats$in_right = samp_stats$clust == low_clust
  samp_stats$clust = NULL

  if (save_filter_stats) {
    make_kmeans_dotplot(samp_stats = samp_stats,
                        plot_dir = file.path(filter_stats_dir, "plots"),
                        bug_name)
  }

  filtered_gf = samp_stats[,.(sample_id, in_right)][gf, on = "sample_id"]
  filtered_gf$abd[filtered_gf$in_right] = 0
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
                     bug_name = NULL) {

  if (!(filtering_method %in% c("med_by_nz_components", "kmeans", "none"))) stop("Specified filtering method not implemented")

  if (filtering_method == "med_by_nz_components"){
    filtered_gf = filter_with_mixture(gf,
                                      samp_stats = samp_stats,
                                      save_filter_stats = save_filter_stats,
                                      filter_stats_dir = filter_stats_dir,
                                      bug_name = bug_name)
  } else if (filtering_method == 'kmeans') {
    filtered_gf = filter_with_kmeans(gf = gf,
                                     samp_stats = samp_stats,
                                     save_filter_stats = save_filter_stats,
                                     filter_stats_dir = filter_stats_dir,
                                     bug_name = bug_name)
  } else if (filtering_method == "none") {
    filtered_gf = gf
    filtered_gf$present = filtered_gf$abd > 0
  }

  if (save_filter_stats & filtering_method != "none") {
    make_line_plot(fgf = filtered_gf,
                   bug_name = bug_name,
                   covariates = covariates,
                   outcome = outcome,
                   plot_dir = file.path(filter_stats_dir, "plots"))
  }

  return(filtered_gf)
}

initial_prevalence_filter = function(gf,
                                     bug_name,
                                     minmax_thresh,
                                     filter_stats_dir,
                                     verbose) {



  # Run the initial prevalence filter
  gf = gf[,.(sample_id, abd,
             varies_enough = sum(abd != 0) < (.N - minmax_thresh) & sum(abd != 0) > minmax_thresh),
          by = gene]

  if (!any(gf$varies_enough)) {
    return(data.table::data.table())
  }

  n_start = dplyr::n_distinct(gf$gene)
  gf = gf[(varies_enough)][,.(gene, sample_id, abd)]
  n_end = dplyr::n_distinct(gf$gene)

  if ((n_end != n_start) & verbose) {
    message(paste0("* Initial prevalence filter dropped ", n_start - n_end, " genes out of ", n_start, " present in the input file."))
    initial_prev_filter = file.path(filter_stats_dir, "initial_prevalence_filter.tsv.gz")
    drop_df = data.table(bug = bug_name,
                         n_dropped_initial_prevalence_filter = n_start - n_end) # TODO write out to file
    if (!file.exists(initial_prev_filter)) {
      readr::write_tsv(drop_df,
                       file = initial_prev_filter)
    } else {
      readr::write_tsv(drop_df,
                       file = initial_prev_filter, append = TRUE)
    }
  }

  gf
}

#' @export
read_and_filter = function(bug_file, meta_cov, # TODO make metadata optional for this step
                           pivot_wide = TRUE,
                           minmax_thresh = 5, # TODO expose this higher up
                           covariates = NULL,
                           outcome = NULL,
                           filtering_method = "kmeans",
                           discard_absent_samples = TRUE,
                           save_filter_stats = TRUE,
                           filter_stats_dir = NULL,
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

  gf = read_bug(bug_file, meta = meta_cov)

  gf = initial_prevalence_filter(gf,
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

  if (verbose) message("* Filtering samples based ")
  filtered_gf = filter_gf(gf,
                          samp_stats = samp_stats,
                          filtering_method = filtering_method,
                          covariates = covariates,
                          outcome = outcome,
                          save_filter_stats = save_filter_stats,
                          filter_stats_dir = filter_stats_dir,
                          bug_name = bug_name) # Might need to reapply the minmax_thresh here

  if (filtering_method != "none") { # TODO separate these four blocks out to a distinct function
    sample_labels = unique(filtered_gf[,.(sample_id, bug_present = !in_right)])
    n_absent = sum(!sample_labels$bug_present)
  }

  if (verbose && filtering_method != "none" && n_absent > 0) {
    message(paste0("* ", n_absent, " samples out of ", nrow(sample_labels), " were determined to not have ", bug_name, " present."))
  }

  if (save_filter_stats && filtering_method != "none") {
    readr::write_tsv(sample_labels,
                     file = file.path(filter_stats_dir, 'labels', paste0('sample_labels_', bug_name, '.tsv.gz')))
  }

  if (discard_absent_samples && filtering_method != "none") {
    filtered_gf = filtered_gf[!(in_right)]
    if (n_absent != 0 & verbose) message("* Samples with no ", bug_name, " discarded.")
  }

  select_cols = c("gene", "present", "sample_id", covariates, outcome)

  if (!is.null(meta_cov)) {
    joined = filtered_gf[meta_cov, on = 'sample_id', nomatch = 0] %>%
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


