get_ellipse = function(mu, sigma) {
  circ_df = tibble(t = seq(0, 2*pi, by = .01),
                   x = cos(t),
                   y = sin(t))
  circ = as.matrix(circ_df[,c('x', 'y')])

  m = (circ %*% sigma) %*% matrix(c(qnorm(.975), 0, 0, qnorm(.975)), nrow = 2)
  m[,1] = m[,1] + mu[1]
  m[,2] = m[,2] + mu[2]
  # as_tibble(m)
  t(m)
}

get_transformed_ellipse = function(id, mix_fit, A, input_scales){
  mat = A %*% get_ellipse(mix_fit$res$mu[[id]], mix_fit$res$sigma[[id]])
  mat = t(mat)
  mat[,1] = mat[,1] + input_scales$m_n_z
  mat[,2] = mat[,2] + input_scales$m_med
  as_tibble(mat) %>% purrr::set_names(c('n_z', 'q50'))
}

#' Make a filter-labelled line plot
#' @param bug_file a path to a gene family file
#' @param meta_file a path to the corresponding metadata file
#' @param fgf a filtered gene family data frame
#' @param bug_name name of the bug
#' @details The required input is either \itemize{ \item{the gene family file
#'   and the metadata file} \item{OR a pre-filtered gene family file}}
#' @inheritParams anpan
#' @export
make_line_plot = function(bug_file = NULL,
                          meta_file = NULL,
                          covariates,
                          outcome,
                          fgf = NULL,
                          bug_name = NULL,
                          plot_dir = NULL) {

  precomputed = !is.null(fgf)
  to_compute = !is.null(bug_file) & !is.null(meta_file)

  if (!(precomputed|to_compute)) {
    # filtered gene families
    fgf = read_and_filter(bug_file,
                          covariates = covariates,
                          outcome = outcome,
                          read_meta(meta_file), pivot_wide = FALSE)

    # bug_name = gsub(".genefamilies.tsv", "", basename(bug_file))
  }

  p = fgf[is.finite(labd)][order(-labd)][, i := 1:(nrow(.SD)), by = sample_id][] %>%
    dplyr::mutate(labelled_as = factor(c('present', 'absent')[in_right + 1],
                                       levels = c("present", "absent"))) %>%
    ggplot(aes(i, labd)) +
    geom_line(aes(group = sample_id,
                  color = labelled_as),
              alpha = .30) +
    labs(x = NULL,
         y = 'log abundance',
         title = bug_name) +
    scale_color_brewer(palette = "Set1") +
    theme_light()

  if (!is.null(plot_dir)) {
    ggsave(p,
           filename = file.path(plot_dir, paste0(bug_name, "_labelled_lines.png")),
           width = 8, height = 7)
  }

  p
}

make_kmeans_dotplot = function(samp_stats,
                               plot_dir = NULL,
                               bug_name = NULL) {

  p = samp_stats %>%
    dplyr::mutate(labelled_as = factor(c('present', 'absent')[in_right + 1],
                                       levels = c("present", "absent"))) %>%
    ggplot(aes(n_z, q50)) +
    geom_point(aes(color = labelled_as),
               alpha = .5) +
    scale_color_brewer(palette = "Set1")+
    guides(color = guide_none()) +
    labs(title = paste0(bug_name, " - labelled by kmeans"),
         x = "Number of zero observations",
         y = 'Median log abundance') +
    theme_light()

  if (!is.null(plot_dir)) {
    ggsave(p,
           filename = file.path(plot_dir, paste0(bug_name, "_hex.png")),
           width = 6, height = 4)
  }

  p
}

#' Make a hex plot of the abundance by zeroness of a gene family dataset
#' @param bug_file path to a gene family file
#' @param meta_file path to a metadata file
#' @param samp_stats data frame of sample statistics
#' @param mix_fit mixture fit object
#' @details The required input is either \itemize{ \item{}{the path to the gene
#'   family file and the metadata file} \item{}{OR a set of sample statistics
#'   and a mixture fit} } The first option is simpler but slower and (very
#'   slightly) stochastic. Defaults to using samp_stats and mix_fit if all three
#'   are supplied (but don't do that anyway).
#' @param bug_name the name of the bug
#' @export
make_hex_plot = function(bug_file = NULL,
                         meta_file = NULL,
                         minmax_thresh = 5,
                         samp_stats = NULL,
                         mix_fit = NULL,
                         bug_name = NULL,
                         plot_dir = NULL) {

  precomputed = (!is.null(samp_stats) & !is.null(mix_fit))
  to_compute = (!is.null(bug_file) & !is.null(meta_file))

  if (to_compute) {
    gf = read_bug(bug_file, meta = read_meta(meta_file))[,.(gene, sample_id, abd, varies_enough = sum(abd != 0) < (.N - minmax_thresh) & sum(abd != 0) > minmax_thresh), by = gene
    ][(varies_enough)
    ][,.(gene, sample_id, abd)]

    samp_stats = get_samp_stats(gf)

    mix_fit = fit_mixture(samp_stats)
  } else if (precomputed) {
    # Do nothing
  } else {
    stop("You somehow misspecified your inputs. Specify either a gene family file and metadata file OR a sample statistics data frame and a mixture fit object")
  }

  em_input = na.omit(samp_stats[,.(sample_id, n_z, q50)])
  em_input$sn_z = scale(em_input$n_z)
  em_input$sq50 = scale(em_input$q50)
  input_scales = em_input[,.(m_n_z = mean(n_z),
                             m_med = mean(q50),
                             s_n_z = sd(n_z),
                             s_med = sd(q50))]
  em_input[,`:=`(centered_n_z = n_z - input_scales$m_n_z,
                 centered_med = q50 - input_scales$m_med)]

  s = t(as.matrix(em_input[1:2, .(sn_z, sq50)]))
  C = t(as.matrix(em_input[1:2, .(centered_n_z, centered_med)]))
  s_inv = solve(s)
  A = C %*% s_inv # This is the transformation between standardized and centered scales

  p = em_input %>%
    ggplot(aes(n_z, q50)) +
    geom_hex(aes(color = ..count..),
             lwd = .1) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    geom_path(color = 'red',
              data = get_transformed_ellipse(id = 1, mix_fit, A, input_scales)) +
    geom_path(color = 'red',
              data = get_transformed_ellipse(id = 2, mix_fit, A, input_scales)) +
    guides(color = guide_none()) +
    labs(title = paste0(bug_name, " - 2 component mixture of median by n_zero observations"),
         x = "Number of zero observations",
         y = 'Median log abundance') +
    theme_light()

  if (!is.null(plot_dir)) {
    ggsave(p,
           filename = file.path(plot_dir, paste0(bug_name, "_hex.png")),
           width = 6, height = 4)
  }

  p
}

make_data_plot = function(res, covariates, outcome, model_input, plot_dir, bug_name,
                          annotation_file = NULL,
                          plot_ext = ".pdf",
                          n_top = 50) {
  if (!is.null(annotation_file)) {
    anno = fread(annotation_file) # must have two columns: gene and annotation
  }

  gene_levels = res[1:n_top,]$gene

  select_cols = c("sample_id", covariates, outcome)
  order_cols = c(outcome, rev(covariates))

  color_bars = model_input %>%
    dplyr::select(dplyr::all_of(select_cols)) %>%
    unique %>%
    setorderv(cols = order_cols, na.last = TRUE)

  color_bars$sample_id = factor(color_bars$sample_id,
                               levels = unique(color_bars$sample_id))

  plot_data = model_input[gene %in% res[1:n_top,]$gene] %>%
    mutate(gene = factor(gene, levels = gene_levels),
           sample_id = factor(sample_id,
                             levels = levels(color_bars$sample_id)))

  n_healthy = sum(color_bars[[outcome]] == 0)

  if (length(covariates) > 2) {
    stop("data plots can't handle more than two covariates right now")
  }

  anno_plot = color_bars %>%
    ggplot(aes(x = sample_id)) +
    geom_tile(aes_string(y = 1, fill = covariates[2])) +
    scale_fill_manual(values = c("male" = "cornflowerblue",
                                 "female" = "sienna1")) +

    ggnewscale::new_scale("fill") +
    geom_tile(aes_string(x = "sample_id", y = 2, fill = covariates[1])) +
    scale_fill_viridis_c() +

    ggnewscale::new_scale("fill") +
    geom_tile(aes_string(y = 3, fill = outcome)) +
    scale_fill_manual(values = c("Healthy" = 'antiquewhite3', 'CRC' = 'indianred2')) +

    coord_cartesian(expand = FALSE) +
    labs(y = NULL,
         x = NULL) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank())

  plot_data = as.data.table(res)[anno[plot_data, on = 'gene'], on = 'gene']

  plot_data$sample_id = factor(plot_data$sample_id,
                              levels = unique(color_bars$sample_id))
  plot_data$gene = factor(plot_data$gene,
                          levels = gene_levels)


  plot_data = plot_data %>%
    mutate(g_lab = map2_chr(gene, annotation, \(.x, .y) if_else(is.na(.y),
                                                                as.character(.x),
                                                                paste(.x, ": ", .y, sep = ""))))

  lab_df = plot_data[,.(gene, g_lab)] %>%
    unique

  ns = n_distinct(plot_data$sample_id)

  pres_plot = plot_data %>%
    mutate(present = as.logical(present)) %>%
    ggplot(aes(y = gene, x = sample_id)) +
    geom_tile(aes(fill = present)) +
    geom_vline(lwd = .5,
               color = 'black',
               xintercept = n_healthy + .5) +
    scale_fill_manual(values = c("FALSE"  = "dodgerblue4", "TRUE"  = "chartreuse")) +
    scale_y_discrete(breaks = lab_df$gene,
                     labels = lab_df$g_lab) +
    labs(x = "samples",
         y = NULL) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank()) +
    coord_cartesian(expand = FALSE)

  p = patchwork::wrap_plots(anno_plot, pres_plot,
                            ncol = 1,
                            heights = c(1, 11),
                            guides = 'collect') +
    patchwork::plot_annotation(title = paste(bug_name, " (n = ", ns, ")", sep = "", collapse = ""))

  w = ifelse(max(nchar(lab_df$g_lab)) > 50,
             12, 8)
  ggsave(plot = p,
         width = w,
         height = 10,
         filename = file.path(plot_dir, paste0(bug_name, "_data.", plot_ext)))

  p
}

make_interval_plot = function() {

}

make_composite_plot = function(bug_file,
                               model_results,
                               covariates,
                               outcome,
                               return_components = FALSE) {

  lp = make_line_plot(bug_file,
                      covariates = covariates,
                      outcome = outcome) #
  hp = make_hex_plot(bug_file)

  dp = make_data_plot(bug_file,
                      covariates = covariates,
                      outcome = outcome,
                      model_results, annotation_file)
  ip = make_int_plot(model_results)

  layout_str = "
  AACD
  BBCD
  "

  p = lp + hp + dp + ip + plot_layout(design = layout_str)
  print(p)
  if (return_components) {
    return(list(lp,hp,dp,ip))
  } else {
    return(p)
  }

}
