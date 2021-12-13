get_ellipse = function(mu, sigma) {
  circ_df = tibble(t = seq(0, 2*pi, by = .01),
                   x = cos(t),
                   y = sin(t))
  circ = as.matrix(circ_df[,c('x', 'y')])

  m = (circ %*% sigma) %*% matrix(c(qnorm(.975), 0, 0, qnorm(.975)), nrow = 2)
  m[,1] = m[,1] + mu[1]
  m[,2] = m[,2] + mu[2]
  as_tibble(m)
}

#' @export
make_line_plot = function(bug_file,
                          meta_file,
                          model_type = "glm") {
  # filtered gene families
  fgf = read_and_label(bug_file,
                       read_meta(meta_file),
                       pivot_wide = model_type == "scone")

  bug_name = gsub(".genefamilies.tsv", "", basename(bug_file))

  fgf[is.finite(labd)][order(-labd)][, i := 1:(nrow(.SD)), by = sampleID][] %>%
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

}

#' Make a hex plot of the abundance by zeroness of a gene family dataset
#' @param bug_file path to a gene family file
#' @param samp_stats data frame of sample statistics
#' @param mix_fit mixture fit object
#' @details The required input is either \itemize{ \item{}{the path to
#'   the gene family file} \item{}{OR a set of sample statistics and
#'   a mixture fit} } The first option is simpler but slower and stochastic.
#'   Defaults to using samp_stats and mix_fit if all three are supplied (but
#'   don't do that anyway).
#' @export
make_hex_plot = function(bug_file = NULL,
                         samp_stats = NULL,
                         mix_fit = NULL,
                         bug_name) {

  if ((is.null(samp_stats) & is.null(mix_fit)) & !is.null(bug_file)) {
    samp_stats = get_samp_stats(gf)

    mix_fit = fit_mixture(samp_stats)
  } else {
    stop("You somehow misspecified your inputs. Specify either a gene family file OR a sample statistics data frame and a mixture fit object")
  }

  em_input = na.omit(samp_stats[,.(sampleID, n_z, q50)])
  em_input$sn_z = scale(em_input$n_z)
  em_input$sq50 = scale(em_input$q50)

  em_input %>%
    ggplot(aes(sn_z, sq50)) +
    geom_hex(aes(color = ..count..),
             lwd = .1) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    geom_path(color = 'red',
              data = get_ellipse(mix_fit$res$mu[[2]], mix_fit$res$sigma[[2]]),
              aes(V1, V2)) +
    geom_path(color = 'red',
              data = get_ellipse(mix_fit$res$mu[[1]], mix_fit$res$sigma[[1]]),
              aes(V1, V2)) +
    guides(color = guide_none()) +
    labs(title = paste0(bug_name, " - 2 component mixture of median by n_zero observations"),
         caption = "Scaled axes to make the EM easier")


}

make_data_plot = function() {

}

make_interval_plot = function() {

}

make_composite_plot = function(bug_file,
                               model_results) {

  lp = make_line_plot(bug_file) #
  hp = make_hex_plot(bug_file)

  dp = make_data_plot(model_results)
  ip = make_int_plot(model_results)

  layout_str = "
  AACD
  BBCD
  "

  lp + hp + dp + ip + plot_layout(design = layout_str)

}
