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
#' @param bug_name
#' @details The required input is either \itemize{ \item{the gene family file
#'   and the metadata file} \item{OR a pre-filtered gene family file}}
#' @export
make_line_plot = function(bug_file = NULL,
                          meta_file = NULL,
                          fgf = NULL,
                          bug_name = NULL,
                          plot_dir = NULL) {

  precomputed = !is.null(fgf)
  to_compute = !is.null(bug_file) & !is.null(meta_file)
  if (!(precomputed|to_compute)) {
    # filtered gene families
    fgf = read_and_filter(bug_file,
                          read_meta(meta_file), pivot_wide = FALSE)

    # bug_name = gsub(".genefamilies.tsv", "", basename(bug_file))
  }

  p = fgf[is.finite(labd)][order(-labd)][, i := 1:(nrow(.SD)), by = sampleID][] %>%
    dplyr::mutate(labelled_as = c('present', 'absent')[in_right + 1]) %>%
    ggplot(aes(i, labd)) +
    geom_line(aes(group = sampleID,
                  color = labelled_as),
              alpha = .33) +
    labs(x = NULL,
         y = 'log abundance',
         title = bug_name) +
    scale_color_brewer(palette = "Set1") +
    theme_light()

  if (!is.null(plot_dir)) {
    ggsave(p,
           filename = file.path(plot_dir, paste0(bug_name, "_labelled_lines.png")),
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
    gf = read_bug(bug_file, meta = read_meta(meta_file))[,.(gene, sampleID, abd, varies_enough = sum(abd != 0) < (.N - minmax_thresh) & sum(abd != 0) > minmax_thresh), by = gene
    ][(varies_enough)
    ][,.(gene, sampleID, abd)]

    samp_stats = get_samp_stats(gf)

    mix_fit = fit_mixture(samp_stats)
  } else if (precomputed) {
    # Do nothing
  } else {
    stop("You somehow misspecified your inputs. Specify either a gene family file and metadata file OR a sample statistics data frame and a mixture fit object")
  }

  em_input = na.omit(samp_stats[,.(sampleID, n_z, q50)])
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

make_data_plot = function() {

}

make_interval_plot = function() {

}

make_composite_plot = function(bug_file,
                               model_results,
                               return_components = FALSE) {

  lp = make_line_plot(bug_file) #
  hp = make_hex_plot(bug_file)

  dp = make_data_plot(bug_file, model_results, annotation_file)
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
