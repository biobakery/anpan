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
                               lambda = c(.75, .25))
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
filter_gf = function(gf, filtering_method = "med_by_nz_components", plot_mix = FALSE,
                     return_labels = FALSE) {

  if (filtering_method != "med_by_nz_components") stop("Only median by zero observation count filtering is currently implemented")

  # TODO allow alternate filtering strategies
  samp_stats = get_samp_stats(gf)

  mix_fit = fit_mixture(samp_stats)

  mix_checks = check_mix_fit(mix_fit)

  if (!mix_checks) {
    stop("Mixture fitting seems to have failed")
  }

  if (plot_mix) {
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

read_and_label = function(bug_file, meta_cov,
                          pivot_wide = TRUE,
                          minmax_thresh = 5) {
  # Would be better to make this function an argument of read_and_filter

  n_lines = R.utils::countLines(bug_file)

  if (n_lines > 161000) {
    stop("This gene family file is huge. Probably E coli. Auto-handling large files isn't implemented yet")
  }

  gf = read_bug(bug_file, meta = meta_cov)[,.(gene, sampleID, abd, varies_enough = sum(abd != 0) < (.N - minmax_thresh) & sum(abd != 0) > minmax_thresh), by = gene
  ][(varies_enough)
  ][,.(gene, sampleID, abd)]

  filtered_gf = filter_gf(gf,
                          filtering_method = "med_by_nz_components",
                          return_labels = TRUE) # Might need to reapply the minmax_thresh here
  filtered_gf
}


read_and_filter = function(bug_file, meta_cov,
                           pivot_wide = TRUE,
                           minmax_thresh = 5) {

  n_lines = R.utils::countLines(bug_file)

  if (n_lines > 161000) {
    stop("This gene family file is huge. Probably E coli. Auto-handling large files isn't implemented yet")
  }

  gf = read_bug(bug_file, meta = meta_cov)[,.(gene, sampleID, abd, varies_enough = sum(abd != 0) < (.N - minmax_thresh) & sum(abd != 0) > minmax_thresh), by = gene
  ][(varies_enough)
  ][,.(gene, sampleID, abd)]

  filtered_gf = filter_gf(gf, filtering_method = "med_by_nz_components") # Might need to reapply the minmax_thresh here

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



# fit mixtures ------------------------------------------------------------

# library(mixtools)
#
# q_df %>%
#   ggplot(aes(n_z, qs)) +
#   geom_hex(aes(color = ..count..)) +
#   scale_fill_viridis_c() +
#   scale_color_viridis_c() +
#   guides(color = guide_none())
#
# circ_df = tibble(t = seq(0, 2*pi, by = .01),
#                  x = cos(t),
#                  y = sin(t))
# circ = as.matrix(circ_df[,c('x', 'y')])
#
# mvnormalmixEM(as.matrix(na.omit(q_df[,.(n_z, qs)])),
#               lambda = c(.2, .8),
#               mu = list(c(500, -.5),
#                         c(3000, .5))) # Too hard to figure out scales
#
# em_input = na.omit(q_df[,.(n_z, qs, quantile)] %>% filter(quantile == ".5") %>% select(-quantile))
#
# mixfit = mvnormalmixEM(as.matrix(mutate_all(em_input, scale)),
#                        lambda = c(.8, .2))
#
# input_scales = em_input %>%
#   summarise(mnz = mean(n_z),
#             m_q = mean(qs),
#             snz = sd(n_z),
#             sq = sd(qs))

get_ell_df = function(mu, sigma){
  path = (circ %*% sigma) %*% matrix(c(input_scales$snz, 0, 0, input_scales$sq), nrow = 2)
  path[,1] = path[,1] + mu[1]*input_scales$mnz
  path[,2] = path[,2] + mu[2]*input_scales$m_q

  path %>%
    as_tibble
}

# q_df %>% filter(quantile == ".5") %>%
#   ggplot(aes(n_z, qs)) +
#   geom_hex(aes(color = ..count..)) +
#   scale_fill_viridis_c() +
#   scale_color_viridis_c() +
#   geom_path(color = 'red',
#             data = get_ell_df(mixfit$mu[[1]], mixfit$sigma[[1]]) ,
#             aes(V1, V2)) +
#   guides(color = guide_none())

get_ell = function(mu, sigma) {
  m = (circ %*% sigma) %*% matrix(c(qnorm(.975), 0, 0, qnorm(.975)), nrow = 2)
  m[,1] = m[,1] + mu[1]
  m[,2] = m[,2] + mu[2]
  as_tibble(m)
}

# mutate_all(em_input, scale) %>%
#   ggplot(aes(n_z, qs)) +
#   geom_hex(aes(color = ..count..),
#            lwd = .1) +
#   scale_fill_viridis_c() +
#   scale_color_viridis_c() +
#   geom_path(color = 'red',
#             data = get_ell(mixfit$mu[[2]], mixfit$sigma[[2]]),
#             aes(V1, V2)) +
#   geom_path(color = 'red',
#             data = get_ell(mixfit$mu[[1]], mixfit$sigma[[1]]),
#             aes(V1, V2)) +
#   guides(color = guide_none()) +
#   labs(title = "Akkermansia result for median with mixtools results",
#        caption = "Scaled axes to make the EM easier")
#
# ggsave('outputs/med_by_z_mixture.png', width = 8, height = 6)
