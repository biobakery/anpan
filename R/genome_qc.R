#' Get genome intersection counts
#' @export
get_genome_intersect_counts = function(unistrain_file) {
  genome_df = fread(unistrain_file,
                    header = TRUE)
  names(genome_df)[1] = "gene"

  starting_count = ncol(genome_df) - 1

  if (ncol(genome_df) < 3) {
    message("Only one genome present in unistrains file, intersections can't be computed.")
    return(NULL)
  }

  if (ncol(genome_df) > 201) {
    message("This unistrains file has many genomes, selecting 200 at random.")
    genome_df = genome_df[,c(1, sample(2:ncol(genome_df),
                                       size = 200)), with = FALSE]
  }

  genome_cols = names(genome_df)[-1]

  genome_df[, (genome_cols) := lapply(.SD, function(.x) .x > 0), .SDcols = genome_cols]

  clust = data.table::transpose(genome_df,
                                keep.names = 'genome',
                                make.names = 'gene') |>
    tibble::column_to_rownames('genome') |>
    dist() |>
    hclust()

  genome_order = clust$labels[clust$order]

  n_gene_df = genome_df[, ..genome_cols] |>
    lapply(sum) |>
    as.data.table() |>
    data.table::melt(measure.vars = genome_cols)

  genome_combn_df = expand.grid(g1 = genome_cols,
                                g2 = genome_cols) |>
    as.data.table() |>
    mutate(g1 = factor(g1, levels = genome_order),
           g2 = factor(g2, levels = genome_order))

  genome_combn_df = genome_combn_df[order(g1, g2)][, combn := mapply(function(.x, .y) paste(c(.x, .y),
                                                                                            collapse = ":"),
                                                                     g1, g2)][!duplicated(combn)][,!"combn"]

  genome_combn_df = genome_combn_df |> # man forget data.table syntax
    dplyr::mutate(dplyr::across(c(g1, g2),
                                as.character))

  int_df = genome_combn_df[g1 != g2][, n_intersect := mapply(function(.x, .y) sum(genome_df[[.x]] & genome_df[[.y]]),
                                                             g1, g2)][]

  off_diag_df = int_df |>
    mutate(g1 = factor(g1, levels = genome_order),
           g2 = factor(g2, levels = genome_order))

  n_gene_df = n_gene_df |>
    mutate(g1 = factor(variable, levels = genome_order),
           g2 = g1)

  return(list(intersection_df = off_diag_df,
              total_gene_count_df = n_gene_df,
              starting_count = starting_count))
}

#' Plot genome intersections
#'
#' @export
plot_genome_intersections = function(unistrain_file,
                                     out_dir = NULL,
                                     file_ext = "pdf",
                                     width = 10,
                                     height = 8) {

  file_id = gsub(".tsv$|.tsv.gz$", "", x = basename(unistrain_file))

  genome_intersect_counts = get_genome_intersect_counts(unistrain_file)

  if (is.null(genome_intersect_counts)) return(NULL)

  off_diag_df = genome_intersect_counts$intersection_df
  n_gene_df   = genome_intersect_counts$total_gene_count_df

  p = ggplot(off_diag_df,
         aes(g1, g2)) +
    geom_tile(aes(fill  = n_intersect,
                  color = n_intersect)) +
    scale_color_viridis_c() +
    scale_fill_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "genome 1",
         y = "genome 2",
         color = "# intersecting\ngenes",
         fill  = "# intersecting\ngenes",
         title = file_id) +
    coord_cartesian(expand = FALSE)

  if (nrow(n_gene_df) < 20) {
    # Label the cell values if there aren't too many
    p = p +
      shadowtext::geom_shadowtext(aes(label = n_intersect)) +
      shadowtext::geom_shadowtext(data = n_gene_df,
                aes(label = value))
  }

  if (!is.null(out_dir)) {

    ggsave(p,
           filename = file.path(out_dir, paste0(file_id, "_genome_intersections.", file_ext)),
           width = width, height = height)
  }

  return(p)

}

#' Intersection histograms
#'
#' @export
plot_genome_intersection_histograms = function(unistrain_file,
                                               out_dir = NULL,
                                               file_ext = "pdf",
                                               width = 10,
                                               height = 5) {

  file_id = gsub(".tsv$|.tsv.gz$", "", x = basename(unistrain_file))

  genome_intersect_counts = get_genome_intersect_counts(unistrain_file)

  if (is.null(genome_intersect_counts)) return(NULL)

  if (genome_intersect_counts$starting_count > 200) {
    caption_str = paste0("From 200 randomly selected genomes from the initial ",
                         genome_intersect_counts$starting_count,
                         " genomes present in the unistrains file.")
  } else {
    caption_str = NULL
  }

  genome_size_hist = genome_intersect_counts$total_gene_count_df |>
    ggplot(aes(value)) +
    geom_histogram(bins = 30) +
    geom_rug() +
    theme_light() +
    labs(x = "Number of genes in genome",
         title = "Total gene counts by genome")

  intersect_hist = genome_intersect_counts$intersection_df |>
    mutate(genome_comb = purrr::map2_chr(g1, g2,
                                         ~paste(sort(c(.x, .y)), collapse = ":"))) |>
    filter(!duplicated(genome_comb)) |>
    ggplot(aes(n_intersect)) +
    geom_histogram(bins = 30) +
    geom_rug() +
    labs(title = "Genome intersection size for each pair") +
    theme_light()

  hists_plot = genome_size_hist + intersect_hist +
    patchwork::plot_annotation(title = file_id,
                               caption = caption_str)

  if (!is.null(out_dir)) {
    ggsave(hists_plot,
           filename = file.path(out_dir, paste0(file_id, "_intersection_hists.", file_ext)),
           width = width, height = height)
  }

  return(hists_plot)
}
