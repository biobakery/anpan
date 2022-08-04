#' Compute a clade effect
#' @description This function computes the average phylogenetic effect of members of a user-defined
#'   clade against that of all non-members.
#' @param clade_members a character vector listing members of the clade of interest
#' @param anpan_pglmm_result a result list from \code{anpan_pglmm()}
#' @param plot_difference logical indicating whether to plot 90% posterior intervals for the three
#'   computed variables
#' @returns a list of two tibbles
#' @export
compute_clade_effects = function(clade_members,
                                 anpan_pglmm_result,
                                 plot_difference = TRUE) {

  n = nrow(anpan_pglmm_result$model_input)

  clade_df = anpan_pglmm_result$model_input |>
    mutate(clade = c("clade_nonmember", "clade_member")[sample_id %in% clade_members + 1],
           param_name = paste("phylo_effect[", 1:n, "]", sep = "")) |>
    tibble::as_tibble()

  phy_effect_df = anpan_pglmm_result$pglmm_fit$draws(format = 'draws_df',
                                         variables = "phylo_effect") |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(matches("phylo_effect"),
                        names_to = 'param_name',
                        values_to = 'value')

  clade_draws = clade_df |>
    select(clade, param_name) |>
    right_join(phy_effect_df, by = 'param_name') |>
    group_by(`.chain`, `.iteration`, `.draw`, clade) |>
    summarise(mean_val = mean(value),
              .groups = "drop_last") |>
    tidyr::pivot_wider(names_from = "clade",
                       values_from = "mean_val") |>
    ungroup() |>
    mutate(clade_difference = clade_member - clade_nonmember)

  clade_summary = clade_draws |>
    posterior::as_draws_df() |>
    posterior::summarise_draws()

  if (plot_difference) {
    p = clade_summary |>
      ggplot(aes(mean, variable)) +
      geom_vline(xintercept = 0,
                 lty = 2,
                 color = 'grey50') +
      geom_point() +
      geom_segment(aes(x = q5,
                       xend = q95,
                       yend = variable)) +
      theme_light() +
      labs(x = "phylogenetic effect",
           y = NULL)

    print(p)
  }

  return(list(clade_draws   = clade_draws,
              clade_summary = clade_summary))
}
