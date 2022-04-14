#' Make a filter-labelled line plot
#' @param bug_file a path to a gene family file
#' @param meta_file a path to the corresponding metadata file
#' @param fgf a filtered gene family data frame
#' @param bug_name name of the bug
#' @param subset_line integer gives the number of points to take along the lines
#' @details The required input is either \itemize{ \item{the gene family file
#'   and the metadata file} \item{OR a pre-filtered gene family file}}
#'   \code{subset_line} is used to make saving plots faster. Plotting thousands
#'   of lines each with tens of thousands of points along them is too much
#'   visual detail and makes saving the plot very slow. Set \code{subset_line}
#'   to 0 to turn off subsetting.
#' @inheritParams anpan
#' @export
plot_lines = function(bug_file = NULL,
                          meta_file = NULL,
                          covariates,
                          outcome,
                          omit_na = FALSE,
                          fgf = NULL,
                          bug_name = NULL,
                          plot_ext = "pdf",
                          plot_dir = NULL,
                          subset_line = 200) {

  precomputed = !is.null(fgf)
  to_compute = !is.null(bug_file) & !is.null(meta_file)

  if (!(precomputed|to_compute)) {
    # filtered gene families
    fgf = read_and_filter(bug_file,
                          covariates = covariates,
                          outcome = outcome,
                          read_meta(meta_file,
                                    select_cols = c("sample_id", outcome, covariates),
                                    omit_na = omit_na), pivot_wide = FALSE)

    # bug_name = gsub(".genefamilies.tsv", "", basename(bug_file))
  }

  plot_df = fgf[is.finite(labd)][order(-labd)][, i := 1:(nrow(.SD)), by = sample_id][] %>%
    dplyr::mutate(labelled_as = factor(c('absent', 'present')[in_right + 1],
                                       levels = c("present", "absent")))

  if (subset_line != 0) {
    plot_df = plot_df[,.SD[floor(seq(1, nrow(.SD), length.out = subset_line))], by = sample_id]
  }

  p = plot_df %>%
    ggplot(aes(i, labd)) +
    geom_line(aes(group = sample_id,
                  color = labelled_as),
              alpha = .30) +
    labs(x = "rank",
         y = 'log abundance',
         title = bug_name,
         color = NULL) +
    scale_color_brewer(palette = "Set1") +
    theme_light()

  if (!is.null(plot_dir)) {
    ggsave(p,
           filename = file.path(plot_dir, paste0(bug_name, "_labelled_lines.", plot_ext)),
           width = 8, height = 7)
  }

  p
}

plot_kmeans_dots = function(samp_stats,
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
    na.omit %>%
    dplyr::mutate(labelled_as = factor(c('absent', 'present')[in_right + 1],
                                       levels = c("present", "absent"))) %>%
    ggplot(aes(n_nz, q50)) +
    geom_point(aes(color = labelled_as),
               alpha = .5) +
    scale_color_brewer(palette = "Set1") +
    scale_x +
    labs(title = paste0(bug_name, " - labelled by kmeans"),
         x = "Number of non-zero observations",
         y = 'Median log abundance',
         color = NULL) +
    theme_light()

  if (!is.null(plot_dir)) {
    ggsave(p,
           filename = file.path(plot_dir, paste0(bug_name, "_kmeans.", plot_ext)),
           width = 6, height = 4)
  }

  p
}

blank_tree = function(clust) {
  ggdendro::ggdendrogram(clust, labels = FALSE, leaf_labels = FALSE) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'cm')) +
    scale_x_continuous(expand = c(0,0))
}

get_cov_color_map = function(unique_covs) {

  disc_scales = list(scale_fill_brewer(palette = "Set1"),
                     scale_fill_brewer(palette = "Set2"))

  cont_scales = list(scale_fill_viridis_c(),
                     scale_fill_viridis_c(option = "magma"))

  col_scales = list(discrete = disc_scales,
                    continuous = cont_scales)

  covs = names(unique_covs)

  cov_types = unique_covs %>%
    imap_dfr(function(.x, .y){tibble(covariate = .y,
                                     is_num = is.numeric(.x),
                                     n_uniq = dplyr::n_distinct(.x))}) %>%
    mutate(cov_type = c('discrete', 'continuous')[((n_uniq >= 5 & (is_num)) + 1)]) %>%
    group_by(cov_type) %>%
    mutate(cov_i = 1:n()) %>%
    ungroup %>%
    mutate(color_scales = map2(cov_type, cov_i,
                               function(.x, .y){col_scales[[.x]][[.y]]}),
           y = 2:(n()+1))

  return(cov_types)
}

plot_color_bars = function(color_bars, model_input,
                          covariates, outcome, binary_outcome) {

  if (binary_outcome) {
    n_healthy = sum(color_bars[[outcome]] == 0)
    n_case = sum(color_bars[[outcome]] == 1)
    outcome_fill_values = c("FALSE" = '#abd9e9', 'TRUE' = '#d73027')
    outcome_fill_scale = scale_fill_manual(values = outcome_fill_values)
    # TODO add color scales too to avoid grey outlines around tiles
  } else{
    outcome_fill_scale = scale_fill_viridis_c(option = "cividis")
  }

  coords = coord_cartesian(expand = FALSE)
  labs_obj = labs(y = NULL,
              x = NULL)
  theme_obj = theme(axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    panel.border = element_blank())

  if (length(covariates) > 0) {
    unique_covs = model_input %>%
      dplyr::select(dplyr::all_of(covariates)) %>%
      unique
    covariate_color_map = get_cov_color_map(unique_covs)
  }

  base_plot = color_bars %>%
    ggplot(aes(x = sample_id)) +
    geom_tile(aes_string(y = 1, fill = outcome)) +
    outcome_fill_scale +
    coords + labs_obj + theme_obj

  p = base_plot

  if (length(covariates) > 0) {
    p = p +

      ggnewscale::new_scale("fill") +
      geom_tile(aes_string(x = "sample_id", y = covariate_color_map$y[1],
                           fill = covariate_color_map$covariate[1])) +
      covariate_color_map$color_scales[[1]]
  }

  if (length(covariates) > 1) {
    p = p +

      ggnewscale::new_scale("fill") +
      geom_tile(aes_string(x = "sample_id",
                           y = covariate_color_map$y[2],
                           fill = covariate_color_map$covariate[2])) +
      covariate_color_map$color_scales[[2]]
  }

  return(p)
}

#' Plot the data for top results
#'
#' @description This funciton makes a tile plot of the top results of a fit
#'   alongside another tile plot showing the covariates included. Optional
#'   annotations can be included.
#' @param res a data frame of model results (from \code{anpan})
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
#' @param show_trees logical to show the trees for the samples (if clustered)
#' @param width width of saved plot in inches
#' @param height height of saved plot in inches
#' @details If included, \code{annotation_file} must be a tsv with two columns:
#'   "gene" and "annotation".
#'
#'   \code{n_top} is ignored if \code{q_threshold} is specified.
#'
#'   When \code{cluster = "none"}, the samples are ordered by metadata and the
#'   genes are ordered by statistical significance.
#' @export
plot_results = function(res, covariates, outcome, model_input, plot_dir = NULL, bug_name,
                             annotation_file = NULL,
                             cluster = 'none',
                             show_trees = FALSE,
                             n_top = 50,
                             q_threshold = NULL,
                             show_intervals = TRUE,
                             width = NULL,
                             height = NULL,
                             plot_ext = "pdf") {

  if (cluster %in% c("none", "genes")) {
    show_trees = FALSE
  }

  binary_outcome = dplyr::n_distinct(model_input[[outcome]]) == 2

  if (!is.null(annotation_file)) {
    # TODO allow annotations to get passed from higher up so you only have to read the (potentially large) annotation file once)
    anno = fread(annotation_file) # must have two columns: gene and annotation
  } else {
    anno = NULL
  }

  n_top = min(n_top, dplyr::n_distinct(res$gene))

  if (!is.null(q_threshold)) {
    gene_levels = res[q_global < q_threshold]$gene

    if (length(gene_levels) == 0) {
      threshold_warning_string = paste0("Note: no genes passed the specified q-value threshold. Displaying the top ", n_top, " genes instead.")
      gene_levels = res[1:n_top,]$gene
      subtitle_str = paste0("Top ", n_top, " hits")
    } else {
      threshold_warning_string = NULL
      subtitle_str = paste0(length(gene_levels), " genes with Q below ", q_threshold)
    }
  } else {
    gene_levels = res[1:n_top,]$gene
    threshold_warning_string = NULL
    subtitle_str = paste0("Top ", n_top, " hits")
  }

  input_mat = model_input %>% dplyr::select('sample_id', all_of(gene_levels)) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()
  input_mat = 1*input_mat # convert to numeric

  if (cluster %in% c('genes', 'both') && ncol(input_mat) > 2) {
    g_clust = hclust(dist(t(input_mat)))

    gene_levels = colnames(input_mat)[g_clust$order]
  }

  select_cols = c("sample_id", covariates, outcome)

  if (cluster %in% c('samples', 'both') && nrow(input_mat) > 2) {
    color_bars = model_input %>%
      dplyr::select(dplyr::all_of(select_cols)) %>%
      unique

    if (binary_outcome) {
      ctls = unique(model_input$sample_id[model_input[[outcome]] == sort(unique(model_input[[outcome]]))[1]])
      ctl_clust = hclust(dist(input_mat[rownames(input_mat) %in% ctls,],
                              method = "binary"))
      ctl_tree = blank_tree(ctl_clust)

      cases = unique(model_input$sample_id[model_input[[outcome]] == sort(unique(model_input[[outcome]]))[2]])
      case_clust = hclust(dist(input_mat[rownames(input_mat) %in% cases,],
                               method = "binary"))

      s_levels = c(rownames(input_mat)[rownames(input_mat) %in% ctls][ctl_clust$order],
                   rownames(input_mat)[rownames(input_mat) %in% cases][case_clust$order])

      case_tree = blank_tree(case_clust)
    } else {
      s_clust = hclust(dist(input_mat,
                            method = 'binary'))
      s_levels = rownames(input_mat)[s_clust$order]
      s_tree = blank_tree(s_clust)
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

  if (binary_outcome) {
    n_healthy = sum(color_bars[[outcome]] == 0)
    n_case = sum(color_bars[[outcome]] == 1)
  }

  if (length(covariates) > 2) {
    stop("data plots can't handle more than two covariates right now")
  }

  anno_plot = plot_color_bars(color_bars, model_input,
                             covariates, outcome, binary_outcome)

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
    plot_data$g_lab = paste(plot_data$gene, plot_data$annotation, sep = ": ")
    no_annotation = is.na(plot_data$annotation)
    plot_data$g_lab[no_annotation] = as.character(plot_data$gene[no_annotation])

    lab_df = plot_data[,.(gene, g_lab)] %>%
      unique

    y_scale = scale_y_discrete(breaks = lab_df$gene,
                               labels = lab_df$g_lab,
                               position = 'right')
    if (is.null(width)) {
      width = ifelse(max(nchar(lab_df$g_lab)) > 50,
                     12, 8)
    }
  } else {
    lab_df = plot_data[,.(gene)] %>% unique
    y_scale = scale_y_discrete(breaks = lab_df$gene,
                               position = 'right')
    if (is.null(width)) {
      width = 8
    }
  }

  ns = dplyr::n_distinct(plot_data$sample_id)
  ng = length(gene_levels)

  glab_frac = ifelse(ng > 50,
                     50^(1/2) / (ng^(1/2)),
                     1)

  if (binary_outcome) {
    black_vline = geom_vline(lwd = .5,
               color = 'black',
               xintercept = n_healthy + .5)
  } else {
    black_vline = NULL
  }

  heatmap_tile = plot_data %>%
    mutate(present = as.logical(present)) %>%
    ggplot(aes(y = gene, x = sample_id)) +
    geom_tile(aes(fill = present)) +
    black_vline +
    scale_fill_manual(values = c("FALSE"  = "dodgerblue4", "TRUE"  = "chartreuse")) +
    y_scale +
    labs(x = "samples",
         y = NULL) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_text(size = ggplot2::rel(glab_frac))) +
    coord_cartesian(expand = FALSE)

  int_plot_df = plot_data[,.(estimate, gene, std.error, `p.value`)] %>% unique %>%
    dplyr::mutate(max_val = estimate + 1.96*std.error,
                  min_val = estimate - 1.96*std.error,
                  p_group = dplyr::case_when(p.value < .001 ~ "***",
                                             p.value < .01  ~ "**",
                                             p.value < .05  ~ "*",
                                             p.value < .1   ~ ".",
                                             p.value < 1    ~ " "))

  est_range = max(int_plot_df$max_val) - min(int_plot_df$min_val)

  star_loc = min(int_plot_df$min_val) - .25*est_range

  int_plot = int_plot_df %>%
    ggplot(aes(estimate, gene)) +
    geom_segment(aes(y = gene,
                     yend = gene,
                     x = min_val,
                     xend = max_val),
                 color = 'grey20') +
    geom_point(color = "grey10") +
    geom_vline(xintercept = 0,
               lty = 2,
               color = 'grey70') +
    geom_text(aes(x = star_loc,
                  y = gene,
                  label = p_group), hjust = 0, vjust = .7) +
    xlim(c(min(0, min(int_plot_df$min_val) - .3*est_range),
           max(0, max(int_plot_df$max_val)))) +
    theme(panel.background = element_rect(fill = "white",
                                          colour = NA),
          panel.border = element_rect(fill = NA,
                                      colour = "grey70", size = rel(1)),
          panel.grid = element_line(colour = "grey92"),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())

  if (binary_outcome && show_trees) {
    n = n_healthy + n_case
    ctl_width =  5 * n_healthy/n
    case_width = 5 * n_case/n
    tree_plot = patchwork::wrap_plots(ctl_tree, case_tree) +
      patchwork::plot_layout(nrow = 1, widths = c(ctl_width, case_width))
  } else if (!binary_outcome && show_trees) {
    n = nrow(input_mat)
    tree_plot = patchwork::wrap_plots(s_tree) +
      patchwork::plot_layout(nrow = 1, widths = 5)
  }

  if (show_intervals && !show_trees) {

    design_str = "
    #AAAAAA
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    " # lol

    p = patchwork::wrap_plots(anno_plot, heatmap_tile, int_plot,
                              ncol = 2,
                              guides = 'collect',
                              design = design_str) +
      patchwork::plot_annotation(title = paste(bug_name, " (n = ", ns, ")", sep = "", collapse = ""),
                                 theme = theme(legend.position = "left"),
                                 caption = threshold_warning_string,
                                 subtitle = subtitle_str)
  } else if (show_intervals && show_trees) {

    design_str = "
    #DDDDDD
    #AAAAAA
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    CBBBBBB
    " # lol

    p = patchwork::wrap_plots( anno_plot, heatmap_tile, int_plot, tree_plot,
                               guides = 'collect',
                               design = design_str) +
      patchwork::plot_annotation(title = paste(bug_name, " (n = ", ns, ")", sep = "", collapse = ""),
                                 theme = theme(legend.position = "left"),
                                 caption = threshold_warning_string,
                                 subtitle = subtitle_str)


  } else if (!show_intervals && show_trees) {
    design_str = "
    CCCCCC
    AAAAAA
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    " # lol

    p = patchwork::wrap_plots( anno_plot, heatmap_tile, tree_plot,
                               guides = 'collect',
                               design = design_str) +
      patchwork::plot_annotation(title = paste(bug_name, " (n = ", ns, ")", sep = "", collapse = ""),
                                 theme = theme(legend.position = "left"),
                                 caption = threshold_warning_string,
                                 subtitle = subtitle_str)

  } else {
    p = patchwork::wrap_plots(anno_plot, heatmap_tile,
                              ncol = 1,
                              heights = c(1, 11),
                              guides = 'collect') +
      patchwork::plot_annotation(title = paste(bug_name, " (n = ", ns, ")", sep = "", collapse = ""),
                                 theme = theme(legend.position = "left"),
                                 caption = threshold_warning_string,
                                 subtitle = subtitle_str)
  }

  if (!is.null(plot_dir)) {
    if (is.null(height)) height = 10
    ggsave(plot = p,
           width = width,
           height = height,
           filename = file.path(plot_dir, paste0(bug_name, "_data.", plot_ext)))
  }

  p
}

safely_plot_results = purrr::safely(plot_results)

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
plot_p_value_histogram = function(all_bug_terms,
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

#' @export
plot_cor_mat = function(cor_mat,
                             bug_name = NULL) {

  if (!is.null(bug_name)) {
    title_str = paste0(bug_name, " tree\ncorrelation matrix")
  } else {
    title_str = NULL
  }

  cor_mat %>%
    as.data.frame %>%
    tibble::rownames_to_column('rn') %>%
    tibble::as_tibble() %>%
    dplyr::mutate(rn = factor(rn,
                       levels = unique(rn))) %>%
    tidyr::pivot_longer(-rn,
                 names_to = 'V2',
                 values_to = 'cor') %>%
    dplyr::mutate(V2 = factor(V2,
                       levels = rev(levels(rn)))) %>%
    ggplot(aes(rn, V2)) +
    geom_tile(aes(fill = cor,
                  color = cor)) +
    labs(x = "sample1",
         y = "sample2",
         title = title_str) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) +
    coord_equal()
}

#' Plot a tree file showing the outcome variable
#' @description Plot a tree file, and show the outcome variable as a colored dot
#'   on the end of each tip.
#' @details Showing the covariates as color bar annotations isn't supported yet.
#' @param return_tree_df if true, return a list containing 1) the plot, 2) the
#'   segment data frame, and 3) the labelled terminal segment data frame.
#'   Otherwise, just return the plot.
#' @inheritParams anpan_pglmm
#' @export
plot_tree = function(tree_file,
                     meta_file,
                     covariates = c("age", "gender"),
                     outcome = 'crc',
                     omit_na = FALSE,
                     verbose = TRUE,
                     return_tree_df = FALSE) {

  if (length(covariates) > 2) {
    stop("more than two covariates is currently not supported")
  }

  olap_list = olap_tree_and_meta(tree_file,
                                 meta_file,
                                 covariates,
                                 outcome,
                                 omit_na,
                                 verbose)

  bug_tree = olap_list[[1]]
  model_input = olap_list[[2]]

  if (dplyr::n_distinct(model_input[[outcome]]) == 2) {
    outcome_color_values = c('#abd9e9', '#d73027')
    names(outcome_color_values) = sort(unique(model_input[[outcome]]))
    outcome_color_scale = scale_color_manual(values = outcome_color_values)
  } else {
    outcome_color_scale = scale_color_viridis_c()
  }

  dend_df = ggdendro::dendro_data(bug_tree %>% phylogram::as.dendrogram())

  seg_df = dend_df$segments %>%
    as_tibble

  tip_df = dend_df$labels %>%
    as_tibble() %>%
    dplyr::select("x", "label")

  terminal_seg_df = seg_df %>%
    filter(x == xend & (x %% 1) == 0) %>%
    group_by(x) %>%
    filter(yend == min(yend)) %>%
    ungroup %>%
    left_join(tip_df, by = "x") %>% # join on tip labels
    left_join(model_input, by = c("label" = "sample_id")) # join on metadata

  p = ggplot(seg_df, aes(x = x, y = yend)) +
    geom_segment(aes(x = x, xend = xend,
                     y = y, yend = yend)) +
    geom_point(data = terminal_seg_df,
               aes_string(color = outcome)) +
    outcome_color_scale +
    scale_x_continuous(breaks = 1:nrow(terminal_seg_df),
                       labels = terminal_seg_df$label) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1,
                                     size = 3),
          panel.background = element_blank(),
          panel.grid = element_blank())

  if (return_tree_df) {
    return(list(tree_plot = p,
                seg_df = seg_df,
                terminal_seg_df = terminal_seg_df))
  } else {
    return(p)
  }
}

#' Plot a tree and the PGLMM posterior predictive
#' @param fit a pglmm fit from \code{anpan_pglmm()}
#' @param labels the ordered tip labels from the tree
#' @inheritParams plot_tree
#' @export
plot_tree_with_post_pred = function(tree_file,
                                    meta_file,
                                    covariates = c("age", "gender"),
                                    outcome = 'crc',
                                    omit_na = FALSE,
                                    verbose = TRUE,
                                    fit = NULL,
                                    labels) {

  if (is.null(fit)) {
    stop("You must provide a pglmm fit to plot the posterior predictive")
  }

  tree_plot = plot_tree(tree_file,
                        meta_file,
                        covariates = covariates,
                        outcome = outcome,
                        omit_na = omit_na,
                        verbose = verbose,
                        return_tree_df = TRUE)

  if (!all(labels == tree_plot$terminal_seg_df$label)) {
    stop('Mismatch between yrep ordering and tree label ordering. This should never happen.')
  }

  yrep_draws = fit$draws(format = "draws_df") %>%
    tibble::as_tibble() %>%
    dplyr::select(matches("yrep"))

  if (dplyr::n_distinct(tree_plot$terminal_seg_df[[outcome]]) == 2) {
    yrep_df = yrep_draws %>%
      dplyr::summarise_all(mean) %>%
      tidyr::pivot_longer(dplyr::everything(), 'param', values_to = "prop_1") %>%
      dplyr::mutate(prop_0 = 1 - prop_1,
                    one = 1.1, # lol
                    zero = -0.1) %>%
      bind_cols(tree_plot$terminal_seg_df) %>%
      mutate(param = factor(param,
                            levels = param))

    yrep_df$outcome = as.numeric(yrep_df[[outcome]])

    yrep_df = yrep_df %>%
      dplyr::select(param, mean_yrep = prop_1, y = outcome) %>%
      tidyr::pivot_longer(-param,
                          names_to = "y_type",
                          values_to = "value")

    yrep_plot = ggplot(yrep_df,
                       aes(x = param)) +
      geom_point(aes(y = value,
                     alpha = y_type)) + # might be better to change these to little bars
      scale_x_discrete(labels = tree_plot$terminal_seg_df$label) +
      scale_alpha_discrete(range = c(.25, 1)) +
      theme(axis.title = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1,
                                       size = 3.5),
            panel.grid = element_blank(),
            panel.background = element_blank())

  } else {
    # TODO handle continuous outcome
    yrep_df = yrep_draws %>%
      posterior::summarise_draws(posterior::default_summary_measures(),
                                 q = ~quantile(.x, probs = c(.25, .75))) %>%
      bind_cols(tree_plot$terminal_seg_df) %>%
      mutate(variable = factor(variable,
                               levels = variable))

    yrep_plot = ggplot(yrep_df, aes(x = variable)) +
      geom_boxplot(aes(ymin = q5,
                       lower = `25%`,
                       middle = median,
                       upper = `75%`,
                       ymax = q95),
                   stat = 'identity') +
      geom_hline(lty = 2,
                 color = 'grey80',
                 yintercept = mean(yrep_df[[outcome]])) +
      geom_point(aes_string(y = outcome,
                            color = outcome)) +
      scale_color_viridis_c() +
      scale_x_discrete(labels = tree_plot$terminal_seg_df$label) +
      labs(y = paste0(outcome, "\n posterior predictive")) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1,
                                       size = 3),
            panel.grid = element_blank(),
            panel.background = element_blank())

  }

  tree_no_labels = tree_plot$tree_plot +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    coord_cartesian(expand = FALSE)

  tree_with_post_pred = tree_no_labels / yrep_plot + plot_layout(heights = c(3,1),
                                                                 guides = "collect")
  return(tree_with_post_pred)


}
