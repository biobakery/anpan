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
                          plot_ext = "pdf",
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
    dplyr::mutate(labelled_as = factor(c('absent', 'present')[in_right + 1],
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
           filename = file.path(plot_dir, paste0(bug_name, "_labelled_lines.", plot_ext)),
           width = 8, height = 7)
  }

  p
}

make_kmeans_dotplot = function(samp_stats,
                               plot_dir = NULL,
                               bug_name = NULL,
                               was_logged = FALSE,
                               plot_ext = "pdf") {

  if (was_logged) {
    scale_x = scale_x_continuous(trans = "log1p")
  } else {
    scale_x = scale_x_continuous()
  }
  p = samp_stats %>%
    dplyr::mutate(labelled_as = factor(c('absent', 'present')[in_right + 1],
                                       levels = c("present", "absent"))) %>%
    ggplot(aes(n_nz, q50)) +
    geom_point(aes(color = labelled_as),
               alpha = .5) +
    scale_color_brewer(palette = "Set1") +
    scale_x +
    labs(title = paste0(bug_name, " - labelled by kmeans"),
         x = "Number of non-zero observations",
         y = 'Median log abundance') +
    theme_light()

  if (!is.null(plot_dir)) {
    ggsave(p,
           filename = file.path(plot_dir, paste0(bug_name, "_kmeans.", plot_ext)),
           width = 6, height = 4)
  }

  p
}

#' Plot the data for top results
#'
#' @description This funciton makes a tile plot of the top results of a fit
#'   alongside another tile plot showing the covariates included. Optional
#'   annotations can be included.
#' @param res a data frame of model results (from \code{anpan} or
#'   \code{anpan_batch})
#' @param covariates character string of the covariates to show
#' @param outcome character string of the outcome variable
#' @param model_input data frame of the model input
#' @param plot_dir directory to write output to
#' @param bug_name character string giving the name to use in the title/output
#'   file
#' @param annotation_file optional path file giving annotations
#' @param plot_ext extension to use for plots
#' @param n_top number of top elements to show from the results
#' @param q_threshold FDR threshold to use for inclusion in the plot.
#' @param cluster axis to cluster. either "none", "samples", "genes", or "both"
#' @details If included, \code{annotation_file} must be a tsv with two columns:
#'   "gene" and "annotation".
#'
#'   \code{n_top} is ignored if \code{q_threshold} is specified.
#'
#'   When \code{cluster = "none"}, the samples are ordered by metadata and the
#'   genes are ordered by statistical significance.
#' @export
make_results_plot = function(res, covariates, outcome, model_input, plot_dir, bug_name,
                             annotation_file = NULL,
                             cluster = 'none',
                             n_top = 50,
                             q_threshold = NULL,
                             show_intervals = TRUE,
                             plot_ext = "pdf") {

  if (!is.null(annotation_file)) {
    # TODO allow annotations to get passed from higher up so you only have to read the (potentially large) annotation file once)
    anno = fread(annotation_file) # must have two columns: gene and annotation
  } else {
    anno = NULL
  }

  if (!is.null(q_threshold)) {
    gene_levels = res[q_global < q_threshold]$gene
  } else {
    gene_levels = res[1:n_top,]$gene
  }

  input_mat = model_input %>% dplyr::select('sample_id', all_of(gene_levels)) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()
  input_mat = 1*input_mat # convert to numeric

  if (cluster %in% c('genes', 'both')) {
    g_clust = hclust(dist(t(input_mat)))

    gene_levels = colnames(input_mat)[g_clust$order]
  }

  select_cols = c("sample_id", covariates, outcome)

  if (cluster %in% c('samples', 'both')) {
    color_bars = model_input %>%
      dplyr::select(dplyr::all_of(select_cols)) %>%
      unique

    if (dplyr::n_distinct(model_input[[outcome]]) == 2) {
      ctls = unique(model_input$sample_id[model_input[[outcome]] == sort(unique(model_input[[outcome]]))[1]])
      ctl_clust = hclust(dist(input_mat[rownames(input_mat) %in% ctls,],
                              method = "binary"))

      cases = unique(model_input$sample_id[model_input[[outcome]] == sort(unique(model_input[[outcome]]))[2]])
      case_clust = hclust(dist(input_mat[rownames(input_mat) %in% cases,],
                               method = "binary"))

      s_levels = c(rownames(input_mat)[rownames(input_mat) %in% ctls][ctl_clust$order],
                   rownames(input_mat)[rownames(input_mat) %in% cases][case_clust$order])
    } else {
      s_clust = hclust(dist(input_mat))
      s_levels = rownames(input_mat)[s_clust$order]
    }

    color_bars$sample_id = factor(color_bars$sample_id,
                                  levels = s_levels)
  } else {
    order_cols = c(outcome, rev(covariates))

    color_bars = model_input %>%
      dplyr::select(dplyr::all_of(select_cols)) %>%
      unique %>%
      setorderv(cols = order_cols, na.last = TRUE)

    color_bars$sample_id = factor(color_bars$sample_id,
                                  levels = unique(color_bars$sample_id))
  }

  model_input = data.table::melt(model_input %>% dplyr::select(all_of(select_cols), all_of(gene_levels)),
                                 id.vars = c(covariates, outcome, "sample_id"),
                                 variable.name = "gene",
                                 value.name = "present")

  plot_data = model_input %>%
    mutate(gene = factor(gene, levels = gene_levels),
           sample_id = factor(sample_id,
                             levels = levels(color_bars$sample_id)))

  n_healthy = sum(color_bars[[outcome]] == 0) # TODO adapt to continuous outcome

  if (length(covariates) > 2) {
    stop("data plots can't handle more than two covariates right now")
  }

  # TODO adapt to continuous outcome
  outcome_fill_values = c("FALSE" = '#abd9e9', 'TRUE' = '#d73027')

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
    scale_fill_manual(values = outcome_fill_values) +

    coord_cartesian(expand = FALSE) +
    labs(y = NULL,
         x = NULL) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank())

  if (!is.null(annotation_file)) {
    plot_data = as.data.table(res)[anno[plot_data, on = 'gene'], on = 'gene']
  } else {
    plot_data = as.data.table(res)[plot_data, on = 'gene']
  }

  plot_data$sample_id = factor(plot_data$sample_id,
                              levels = levels(color_bars$sample_id))
  plot_data$gene = factor(plot_data$gene,
                          levels = rev(gene_levels))


  if (!is.null(annotation_file)) {
    plot_data[, g_lab := paste(gene, annotation, sep = ": ")]
    no_annotation = is.na(plot_data$annotation)
    plot_data$g_lab[no_annotation] = plot_data$gene[no_annotation]

    lab_df = plot_data[,.(gene, g_lab)] %>%
      unique

    y_scale = scale_y_discrete(breaks = lab_df$gene,
                               labels = lab_df$g_lab,
                               position = 'right')
    w = ifelse(max(nchar(lab_df$g_lab)) > 50,
               12, 8)
  } else {
    lab_df = plot_data[,.(gene)] %>% unique
    y_scale = scale_y_discrete(breaks = lab_df$gene,
                               position = 'right')
    w = 8
  }

  ns = n_distinct(plot_data$sample_id)

  pres_plot = plot_data %>%
    mutate(present = as.logical(present)) %>%
    ggplot(aes(y = gene, x = sample_id)) +
    geom_tile(aes(fill = present)) +
    geom_vline(lwd = .5,
               color = 'black',
               xintercept = n_healthy + .5) +
    scale_fill_manual(values = c("FALSE"  = "dodgerblue4", "TRUE"  = "chartreuse")) +
    y_scale +
    labs(x = "samples",
         y = NULL) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank()) +
    coord_cartesian(expand = FALSE)

  if (show_intervals) {
    int_plot = plot_data[,.(estimate, gene, std.error)] %>% unique %>%
      ggplot(aes(estimate, gene)) +
      geom_segment(aes(y = gene,
                       yend = gene,
                       x = estimate - 1.96*std.error,
                       xend = estimate + 1.96*std.error),
                   color = 'grey20') +
      geom_point(color = "grey10") +
      geom_vline(xintercept = 0,
                 lty = 2,
                 color = 'grey30') +
      theme(panel.background = element_rect(fill = "white",
                                            colour = NA),
            panel.border = element_rect(fill = NA,
                                        colour = "grey70", size = rel(1)),
            panel.grid = element_line(colour = "grey87"),
            panel.grid.major.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank())

    design_str = "
    #AAAAA
    CBBBBB
    CBBBBB
    CBBBBB
    CBBBBB
    CBBBBB
    CBBBBB
    CBBBBB
    CBBBBB
    CBBBBB
    CBBBBB
    " # lol

    p = patchwork::wrap_plots(anno_plot, pres_plot, int_plot,
                              ncol = 2,
                              guides = 'collect',
                              design = design_str) +
      patchwork::plot_annotation(title = paste(bug_name, " (n = ", ns, ")", sep = "", collapse = ""),
                                 theme = theme(legend.position = "left"))
  } else {
    p = patchwork::wrap_plots(anno_plot, pres_plot,
                              ncol = 1,
                              heights = c(1, 11),
                              guides = 'collect') +
      patchwork::plot_annotation(title = paste(bug_name, " (n = ", ns, ")", sep = "", collapse = ""),
                                 theme = theme(legend.position = "left"))
  }

  ggsave(plot = p,
         width = w,
         height = 10,
         filename = file.path(plot_dir, paste0(bug_name, "_data.", plot_ext)))

  p
}

make_interval_plot = function(res,
                              covariates, outcome, model_input, plot_dir, bug_name,
                              annotation_file = NULL,
                              plot_ext = "pdf",
                              n_top = 50,
                              q_threshold = NULL) {
  if (!is.null(annotation_file)) {
    # TODO allow annotations to get passed from higher up so you only have to read the (potentially large) annotation file once)
    anno = fread(annotation_file) # must have two columns: gene and annotation
  }

  if (!is.null(q_threshold)) {
    plot_data = res[q_global < q_threshold]
    gene_levels = plot_data$gene
  } else {
    plot_data = res[1:n_top,]
    gene_levels = plot_data$gene
  }
  plot_data$gene = factor(plot_data$gene,
                          levels = rev(gene_levels))



}

#' Make a p-value histogram
#'
#' @description This function makes a p-value histogram from a collection of
#'   bug:gene glm fits.
#'
#' @param all_bug_terms a data frame of bug:gene glm fits
#' @param out_dir string giving the output directory
#' @param plot_ext string giving the extension to use
#' @details The plot will be written out to \code{aaa_p_value_histogram.<ext>}
#'   in the specified output directory. The "aaa" is there for alphabetical
#'   superiority.
#'
#'   If you don't understand the purpose of this type of plot,
#'   \href{http://varianceexplained.org/statistics/interpreting-pvalue-histogram/}{this
#'   blog post by David Robinson} has a lot of helpful information.
#' @export
make_p_value_histogram = function(all_bug_terms,
                                  out_dir = NULL,
                                  plot_ext = "pdf",
                                  n_bins = 50) {
  p = all_bug_terms %>%
    ggplot(aes(`p.value`)) +
    geom_histogram(breaks = seq(0, 1, length.out = n_bins)) +
    labs(title = "p-value histogram for all bug:gene glm fits",
         subtitle = paste0("There are ",
                           dplyr::n_distinct(all_bug_terms$bug_name), " unique bugs, ",
                           dplyr::n_distinct(all_bug_terms$gene), " unique genes, and ",
                           nrow(unique(all_bug_terms[,.(bug_name, gene)])), " bug:gene combinations.")) +
    theme_light()

  if (!is.null(out_dir)) {
    ggsave(plot = p,
           filename = file.path(out_dir, paste0("aaa_p_value_histogram.", plot_ext)),
           width = 8, height = 6)
  }

  p
}

make_composite_plot = function(bug_file,
                               model_results,
                               covariates,
                               outcome,
                               return_components = FALSE) {

  lp = make_line_plot(bug_file,
                      covariates = covariates,
                      outcome = outcome) #

  dp = make_data_plot(bug_file,
                      covariates = covariates,
                      outcome = outcome,
                      model_results, annotation_file)
  ip = make_interval_plot(model_results)

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

#' @export
make_cov_mat_plot = function(cov_mat,
                             bug_name = NULL) {

  if (!is.null(bug_name)) {
    title_str = paste0(bug_name, " tree\nas a covariance matrix")
  } else {
    title_str = NULL
  }

  cov_mat %>%
    as.data.frame %>%
    tibble::rownames_to_column('rn') %>%
    tibble::as_tibble() %>%
    dplyr::mutate(rn = factor(rn,
                       levels = unique(rn))) %>%
    tidyr::pivot_longer(-rn,
                 names_to = 'V2',
                 values_to = 'cov') %>%
    dplyr::mutate(V2 = factor(V2,
                       levels = rev(levels(rn)))) %>%
    ggplot(aes(rn, V2)) +
    geom_tile(aes(fill = cov)) +
    labs(x = "sample1",
         y = "sample2",
         title = title_str) +
    scale_fill_viridis_c() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) +
    coord_equal()
}
