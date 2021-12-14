# Do outlier detection on number of confidently detected genes
# library(data.table)
# library(tidyverse); theme_set(theme_light())
# library(patchwork)
# library(furrr);

# windoze = grepl("Windows",
#                 Sys.info()['sysname'])
#
# if (windoze) {
#   gf_files = list.files("../../../../strain_stats/data/genefamilies/",
#                         full.names = TRUE)
#   out_dir = "outputs/"
#   meta = read_tsv('data/CRC_analysis_metadata_final_version.tsv')
#   plan(multisession, workers = 6)
# } else {
#   gf_files = list.files("/n/holystore01/LABS/huttenhower_lab/Users/yyan/humann2_bug_gene/franzosa/bb3_version/output/genefamilies",
#                         full.names = TRUE)
#   out_dir = "/n/hutlab12_nobackup/users/aghazi/strain_stats/outputs/"
#   meta = read_tsv('/n/home11/aghazi/strain_stats/data/CRC_analysis_metadata_final_version.tsv')
#   plan(multisession, workers = 32)
#
#   if (!dir.exists('outputs/n_z_by_quantiles')) {
#     dir.create('outputs/n_z_by_quantiles')
#     dir.create('outputs/n_z_by_quantiles/hexs')
#     dir.create('outputs/n_z_by_quantiles/hists')
#   }
# }

get_q_large = function(gf_file) {
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
           variable.name = 'sampleID',
           value.name = "abd")

    gf$sampleID = factor(gf$sampleID,
                         levels = s_ids)

    if (i == 1){
      chunk_set = gf[, labd := log10(abd)][, .(n_z = sum(abd == 0),
                                               n_nz = sum(abd > 0),
                                               qs = get_qs_and_n_hi_lo(labd[abd > 0]),
                                               quantile = c('q10', 'q50', 'q90', 'n_hi', 'n_lo')), by = sampleID] %>%
        data.table::dcast(sampleID + n_z + n_nz ~ quantile, value.var = 'qs') %>%
        .[, n_diff := (n_lo - n_hi), by = sampleID] %>%
        .[]
    } else {
      chunk_set = rbind(chunk_set,
                        chunk_set = gf[, labd := log10(abd)][, .(n_z = sum(abd == 0),
                                                                 n_nz = sum(abd > 0),
                                                                 qs = get_qs_and_n_hi_lo(labd[abd > 0]),
                                                                 quantile = c('q10', 'q50', 'q90', 'n_hi', 'n_lo')), by = sampleID] %>%
                          data.table::dcast(sampleID + n_z + n_nz ~ quantile, value.var = 'qs') %>%
                          .[, n_diff := (n_lo - n_hi), by = sampleID] %>%
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
                               quantile = c('q10', 'q50', 'q90', 'n_hi', 'n_lo')), by = sampleID] %>%
    data.table::dcast(sampleID + n_z + n_nz ~ quantile, value.var = 'qs') %>%
    .[, n_diff := (n_lo - n_hi), by = sampleID] %>%
    .[]
}

fit_mixture = function(samp_stats) {
  em_input = na.omit(samp_stats[,.(sampleID, n_z, q50)])
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

  samp_df[,.(sampleID, in_right)]
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

#' @export
filter_gf = function(gf, filtering_method = "med_by_nz_components",
                     plot_mix = FALSE,
                     save_filter_stats = FALSE,
                     filter_stats_dir = NULL,
                     bug_name = NULL) {

  if (filtering_method != "med_by_nz_components") stop("Only median by zero observation count filtering is currently implemented")

  # TODO allow alternate filtering strategies
  samp_stats = get_samp_stats(gf)

  mix_fit = fit_mixture(samp_stats)

  mix_checks = check_mix_fit(mix_fit)

  if (!mix_checks) {
    stop("Mixture fitting seems to have failed")
  }

  if (save_filter_stats) {
    if (!dir.exists(filter_stats_dir)) {
      stop("filter_stats_dir doesn't exist!")
    }
    mix_fit_res = mix_fit$res
    save(samp_stats, mix_fit_res,
         file = paste0(filter_stats_dir, bug_name, ".RData"))
  }

  if (plot_mix) {
    make_hex_plot(samp_stats = samp_stats,
                  mix_fit = mix_fit)
    stop("Plotting the estimated mixture components isn't implemented yet")
  }

  mix_labels = get_component_densities(mix_fit)
  label_df = mix_labels[samp_stats, on = "sampleID"]
  label_df$in_right[is.na(label_df$in_right)] = TRUE # all zero samples don't have the species

  filtered_gf = label_df[,.(sampleID, in_right)][gf, on = "sampleID"]
  filtered_gf$abd[filtered_gf$in_right] = 0
  filtered_gf$present = filtered_gf$abd > 0

  return(filtered_gf)
}

read_and_filter = function(bug_file, meta_cov,
                           pivot_wide = TRUE,
                           minmax_thresh = 5,
                           save_filter_stats = FALSE,
                           filter_stats_dir = NULL) {

  n_lines = R.utils::countLines(bug_file)

  if (n_lines > 161000) {
    stop("This gene family file is huge. Probably E coli. Auto-handling large files isn't implemented yet")
  }

  if (save_filter_stats & is.null(filter_stats_dir)){
    stop("You said save the filter statistics but didn't provide filter_stats_dir")
  }

  bug_name = get_bug_name(bug_file)

  gf = read_bug(bug_file, meta = meta_cov)[,.(gene, sampleID, abd, varies_enough = sum(abd != 0) < (.N - minmax_thresh) & sum(abd != 0) > minmax_thresh), by = gene
  ][(varies_enough)
  ][,.(gene, sampleID, abd)]

  filtered_gf = filter_gf(gf, filtering_method = "med_by_nz_components",
                          save_filter_stats = save_filter_stats,
                          filter_stats_dir = filter_stats_dir,
                          bug_name = bug_name) # Might need to reapply the minmax_thresh here

  joined = filtered_gf[meta_cov, on = 'sampleID', nomatch = 0][,.(gene, present, sampleID, age, gender, crc)]

  if (pivot_wide) {
    wide_dat = dcast(joined,
                     age + gender + sampleID + crc ~ gene, # TODO generalize covariates
                     value.var = 'present')

    return(wide_dat)
  } else {
    return(joined)
  }
}


