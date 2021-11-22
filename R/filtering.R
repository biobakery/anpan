# Do outlier detection on number of confidently detected genes
library(data.table)
library(tidyverse); theme_set(theme_light())
library(patchwork)
library(furrr);

windoze = grepl("Windows",
                Sys.info()['sysname'])

if (windoze) {
  gf_files = list.files("../../../../strain_stats/data/genefamilies/",
                        full.names = TRUE)
  out_dir = "outputs/"
  meta = read_tsv('data/CRC_analysis_metadata_final_version.tsv')
  plan(multisession, workers = 6)
} else {
  gf_files = list.files("/n/holystore01/LABS/huttenhower_lab/Users/yyan/humann2_bug_gene/franzosa/bb3_version/output/genefamilies",
                        full.names = TRUE)
  out_dir = "/n/hutlab12_nobackup/users/aghazi/strain_stats/outputs/"
  meta = read_tsv('/n/home11/aghazi/strain_stats/data/CRC_analysis_metadata_final_version.tsv')
  plan(multisession, workers = 32)

  if (!dir.exists('outputs/n_z_by_quantiles')) {
    dir.create('outputs/n_z_by_quantiles')
    dir.create('outputs/n_z_by_quantiles/hexs')
    dir.create('outputs/n_z_by_quantiles/hists')
  }
}

read_bug = function(bug_file) {

  nc = readLines(bug_file,
                 n = 1) %>%
    strsplit('\t') %>%
    .[[1]] %>%
    length

  gf = fread(bug_file,
             colClasses = list(character = 1, numeric = 2:nc)) %>%
    dplyr::select_all(~gsub("_Abundance-RPKs", "", .))

  names(gf)[1] = "gene_spec"
  gf = gf %>%
    tidyr::separate(gene_spec,
                    into = c("u", "s"),
                    sep = "\\|")

  s_ids = names(gf)[-c(1, 2)]

  gf = gf %>%
    melt(id.vars = c("u", "s"),
         variable.name = 'sampleID',
         value.name = "abd")

  gf$sampleID = factor(gf$sampleID,
                       levels = s_ids)
  gf
}

# i = 6
gf = read_bug(gf_files[i])

q_df = gf[, labd := log10(abd)][, .(n_z = sum(abd == 0),
                                    n_nz = sum(abd > 0),
                                    qs = quantile(labd[abd > 0], probs = c(.1, .5, .9)),
                                    quantile = c('.1', '.5', '.9')), by = sampleID]

# int_width = q_df[,.(q80 = qs[quantile == ".9"] - qs[quantile == ".1"], n_nz = n_nz[1]), by = sampleID]

# q_hists = q_df %>%
#   ggplot(aes(qs)) +
#   geom_histogram() +
#   facet_wrap('quantile', labeller = label_both)

# q_df %>%
#   ggplot(aes(n_nz, qs)) +
#   geom_point(size = 1) +
#   facet_wrap('quantile', labeller = label_both)

# q_by_nz_hex = q_df %>%
#   ggplot(aes(n_nz, qs)) +
#   geom_hex() +
#   scale_fill_viridis_c() +
#   facet_wrap('quantile', labeller = label_both)

get_q_large = function(gf_file){
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

    names(gf)[1] = "gene_spec"
    gf = gf %>%
      tidyr::separate(gene_spec,
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
                                               quantile = c('.1', '.5', '.9', 'n_hi', 'n_lo')), by = sampleID] %>%
        .[, n_diff := (qs[quantile == "n_lo"] - qs[quantile == "n_hi"]), by = sampleID] %>%
        .[]
    } else {
      chunk_set = rbind(chunk_set,
                        chunk_set = gf[, labd := log10(abd)][, .(n_z = sum(abd == 0),
                                                                 n_nz = sum(abd > 0),
                                                                 qs = get_qs_and_n_hi_lo(labd[abd > 0]),
                                                                 quantile = c('.1', '.5', '.9', 'n_hi', 'n_lo')), by = sampleID] %>%
                          .[, n_diff := (qs[quantile == "n_lo"] - qs[quantile == "n_hi"]), by = sampleID] %>% .[])
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

get_samp_quantiles = function(gf){
  gf[, labd := log10(abd)][, .(n_z = sum(abd == 0),
                               n_nz = sum(abd > 0),
                               qs = get_qs_and_n_hi_lo(labd[abd > 0]),
                               quantile = c('.1', '.5', '.9', 'n_hi', 'n_lo')), by = sampleID] %>%
    .[, n_diff := (qs[quantile == "n_lo"] - qs[quantile == "n_hi"]), by = sampleID] %>%
    .[]
}

make_q_plot = function(gf_file) {
  bug_name = gsub(".genefamilies.tsv", "", basename(gf_file))
  n_lines = R.utils::countLines(gf_file)

  if (n_lines > 161000){
    q_df = get_q_large(gf_file)
  } else {
    gf = read_bug(gf_file)
    q_df = get_samp_quantiles(gf)
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

safely_plot = purrr::safely(make_q_plot)

# pl = future_map(gf_files %>% head,
#                 safely_plot,
#                 .progress = TRUE)


# fit mixtures ------------------------------------------------------------

library(mixtools)

q_df %>%
  ggplot(aes(n_z, qs)) +
  geom_hex(aes(color = ..count..)) +
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  guides(color = guide_none())

circ_df = tibble(t = seq(0, 2*pi, by = .01),
                 x = cos(t),
                 y = sin(t))
circ = as.matrix(circ_df[,c('x', 'y')])

mvnormalmixEM(as.matrix(na.omit(q_df[,.(n_z, qs)])),
              lambda = c(.2, .8),
              mu = list(c(500, -.5),
                        c(3000, .5))) # Too hard to figure out scales

em_input = na.omit(q_df[,.(n_z, qs, quantile)] %>% filter(quantile == ".5") %>% select(-quantile))

mixfit = mvnormalmixEM(as.matrix(mutate_all(em_input, scale)),
                       lambda = c(.8, .2))

input_scales = em_input %>%
  summarise(mnz = mean(n_z),
            m_q = mean(qs),
            snz = sd(n_z),
            sq = sd(qs))

get_ell_df = function(mu, sigma){
  path = (circ %*% sigma) %*% matrix(c(input_scales$snz, 0, 0, input_scales$sq), nrow = 2)
  path[,1] = path[,1] + mu[1]*input_scales$mnz
  path[,2] = path[,2] + mu[2]*input_scales$m_q

  path %>%
    as_tibble
}

q_df %>% filter(quantile == ".5") %>%
  ggplot(aes(n_z, qs)) +
  geom_hex(aes(color = ..count..)) +
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  geom_path(color = 'red',
            data = get_ell_df(mixfit$mu[[1]], mixfit$sigma[[1]]) ,
            aes(V1, V2)) +
  guides(color = guide_none())

get_ell = function(mu, sigma) {
  m = (circ %*% sigma) %*% matrix(c(qnorm(.975), 0, 0, qnorm(.975)), nrow = 2)
  m[,1] = m[,1] + mu[1]
  m[,2] = m[,2] + mu[2]
  as_tibble(m)
}

mutate_all(em_input, scale) %>%
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
  labs(title = "Akkermansia result for median with mixtools results",
       caption = "Scaled axes to make the EM easier")

ggsave('outputs/med_by_z_mixture.png', width = 8, height = 6)
