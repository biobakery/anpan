
get_ell_df = function(mu, sigma){
  circ_df = tibble(t = seq(0, 2*pi, by = .01),
                   x = cos(t),
                   y = sin(t))
  circ = as.matrix(circ_df[,c('x', 'y')])
  path = (circ %*% sigma) %*% matrix(c(input_scales$snz, 0, 0, input_scales$sq), nrow = 2)
  path[,1] = path[,1] + mu[1]*input_scales$mnz
  path[,2] = path[,2] + mu[2]*input_scales$m_q

  path %>%
    as_tibble
}

get_ellipse = function(mu, sigma) {
  m = (circ %*% sigma) %*% matrix(c(qnorm(.975), 0, 0, qnorm(.975)), nrow = 2)
  m[,1] = m[,1] + mu[1]
  m[,2] = m[,2] + mu[2]
  as_tibble(m)
}

make_q_plot = function(gf_file) {
  bug_name = gsub(".genefamilies.tsv", "", basename(gf_file))
  n_lines = R.utils::countLines(gf_file)

  if (n_lines > 161000){
    q_df = get_q_large(gf_file)
  } else {
    gf = read_bug(gf_file)
    q_df = get_samp_stats(gf)
  }

  int_width = q_df[,.(q80 = qs[quantile == ".9"] - qs[quantile == ".1"], n_nz = n_nz[1], n_z = n_z[1]), by = sampleID]

  q80_hist = int_width %>% ggplot(aes(q80)) + geom_histogram() + theme_light()
  q80_hex = int_width %>%
    ggplot(aes(n_z, q80)) +
    geom_hex(aes(color = ..count..)) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    guides(color = guide_none())+ theme_light()

  n80_hist = q_df[n_diff != 0][quantile == 'n_hi'] %>% ggplot(aes(n_diff)) + geom_histogram() + labs(x = 'n80') + theme_light()
  # n80_hex = q_df[n_diff != 0][quantile == 'n_hi'] %>% # OBVIOUSLY it's just a linear relationship
  #   dplyr::rename(n80 = n_diff) %>%
  #   ggplot(aes(n_z, n80)) +
  #   geom_hex(aes(color = ..count..)) +
  #   scale_fill_viridis_c() +
  #   scale_color_viridis_c() +
  #   guides(color = guide_none())+ theme_light()

  q_hists = q_df %>% filter(!grepl("n", quantile)) %>%
    ggplot(aes(qs)) +
    geom_histogram() +
    facet_wrap('quantile', labeller = label_both) +
    theme_light()

  # q_df %>%
  #   ggplot(aes(n_z, qs)) +
  #   geom_point(size = 1) +
  #   facet_wrap('quantile', labeller = label_both)

  scale_fun = ifelse(n_lines > 161000,
                     scale_x_log10,
                     scale_x_continuous)

  title_str = NULL
  if (n_lines > 161000) title_str = 'log scale on x-axis!'

  q_by_z_hex = q_df %>% filter(!grepl("n", quantile)) %>%
    ggplot(aes(n_z, qs)) +
    geom_hex(aes(color = ..count..)) +
    scale_fun() +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    facet_wrap('quantile', labeller = label_both) +
    guides(color = guide_none()) +
    theme_light() +
    labs(title = title_str, y = "quantile_value", x = "number of zero observations")



  layout_str = "
  1##
  222
  "
  ggsave(plot =  q80_hex + q_by_z_hex +  plot_annotation(title = bug_name) + plot_layout(design = layout_str),
         filename = paste0("outputs/n_z_by_quantiles/hexs/", bug_name, ".png"),
         width = 12,
         height = 7)
  ggsave(plot = (q80_hist + n80_hist) / q_hists + plot_annotation(title = bug_name),
         filename = paste0("outputs/n_z_by_quantiles/hists/", bug_name, ".png"),
         width = 12,
         height = 7)

  NULL
}

make_mix_plot = function(gf_file){
  bug_name = gsub(".genefamilies.tsv", "", basename(gf_file))
  n_lines = R.utils::countLines(gf_file)

  if (n_lines > 161000){
    q_df = get_q_large(gf_file)
  } else {
    gf = read_bug(gf_file)
    q_df = get_samp_quantiles(gf)
  }

  em_input = na.omit(q_df[,.(n_z, qs, quantile)] %>% filter(quantile == ".5") %>% select(-quantile))

  mixfit = mvnormalmixEM(as.matrix(mutate_all(em_input, scale)),
                         lambda = c(.8, .2))

  input_scales = em_input %>%
    summarise(mnz = mean(n_z),
              m_q = mean(qs),
              snz = sd(n_z),
              sq = sd(qs))

  p = mutate_all(em_input, scale) %>%
    ggplot(aes(n_z, qs)) +
    geom_hex(aes(color = ..count..),
             lwd = .1) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    geom_path(color = 'red',
              data = get_ell(mixfit$mu[[2]], mixfit$sigma[[2]]),
              aes(V1, V2)) +
    geom_path(color = 'red',
              data = get_ell(mixfit$mu[[1]], mixfit$sigma[[1]]),
              aes(V1, V2)) +
    guides(color = guide_none()) +
    labs(title = paste0(bug_name, " - 2 component mixture of median by n_zero observations"),
         caption = "Scaled axes to make the EM easier")

  ggsave(p,
         filename = paste0('outputs/2_comp_mix_plots/', bug_name, '.png'),
         width = 8, height = 6)
  NULL
}

safely_mix = purrr::safely(make_mix_plot)
