plot_genome_overlaps = function(unistrain_file,
                                out_dir = NULL,
                                file_ext = "pdf") {

  file_id = basename(unistrain_file) |> gsub(".tsv$|.tsv.gz$", "", x = _)

  genome_df = fread(unistrain_file)
  names(genome_df)[1] = "gene"

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


  p = ggplot(off_diag_df,
         aes(g1, g2)) +
    geom_tile(aes(fill  = n_intersect,
                  color = n_intersect)) +
    scale_color_viridis_c() +
    scale_fill_viridis_c() +
    theme_light() +
    labs(x = "genome 1",
         y = "genome 2")

  if (length(genome_order) < 20) {
    p = p + geom_text(aes(label = n_intersect)) +
      geom_text(data = n_gene_df,
                aes(label = value))
  }

  if (!is.null(out_dir)) {
    ggsave(p,
           filename = file.path(out_dir, paste0(file_id, "_genome_overlaps.", file_ext)),
           width = 10, height = 8)
  }

  return(p)

}
