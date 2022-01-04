# library(tidyverse); theme_set(theme_light())
# library(brms)
# library(patchwork)
# library(shadowtext)
# library(ggbeeswarm)
#
# meta_adult = read_csv("data/current_cross_sectional_metadata.csv") %>%
#   mutate(outcome = c("case", "control")[grepl("HC|NIJP", x = diagnosis) + 1]) %>%
#   select(sample_id = X1, age, gender, outcome)
# meta_jia = read_csv("data/current_jia_metadata.csv") %>%
#   mutate(outcome = c("case", "control")[grepl("HC", x = subtype) + 1]) %>%
#   select(sample_id = SampleID, age, gender, outcome)
#
# meta = rbind(meta_jia, meta_adult)
#
# trees = list.files("../../../../strain_stats/data/age_trees",
#                    pattern = ".tre$",
#                    full.names = TRUE)
#
# bug_file = trees[1]
#
# bug_tree = ape::read.tree(bug_file)
# bug_tree$tip.label = bug_tree$tip.label %>%
#   gsub("_bowtie2", "", .)
# cov_mat = ape::vcv.phylo(bug_tree)
#
# meta_gender_n = meta %>% count(gender, outcome)
#
# model_input %>%
#   ggplot(aes(gender, outcome)) +
#   geom_jitter(width = .2, height = .2, color = 'lightblue3') +
#   geom_shadowtext(data = meta_gender_n,
#                   aes(label = n),
#                   size = 10)
# ggsave("outputs/gender_by_outcome.png", width = 6, height = 5)
#
# model_input %>%
#   ggplot(aes(age, outcome)) +
#   geom_beeswarm(groupOnX = FALSE, cex = 1.7)
# ggsave("outputs/age_by_outcome.png", width = 7, height = 5)
#
# model_input = meta %>%
#   filter(sample_id %in% bug_tree$tip.label) %>%
#   mutate(outcome = dplyr::case_when(outcome == 'case' ~ 1,
#                                     outcome == "control" ~ 0))
#
# base_age_model = brm(outcome ~ age,
#                      data = model_input,
#                      family = bernoulli(),
#
#                      prior = c(
#                        prior(normal(0, 50), "Intercept"),
#                        prior(student_t(3, 0, 20), "sd")
#                      ),
#                      cores = 4)
#
# base_tree_model = brm(outcome ~ (1|gr(sample_id, cov = cov_mat)),
#                       data = model_input,
#                       family = bernoulli(),
#                       data2 = list(cov_mat = cov_mat),
#                       prior = c(
#                         prior(normal(0, 50), "Intercept"),
#                         prior(student_t(3, 0, 20), "sd")
#                       ),
#                       cores = 4)
#
# base_tree_model$fit %>% as.data.frame() %>%
#   posterior::summarise_draws() %>% arrange(-abs(mean))

#' @title Run a phylogenetic generalized linear mixed model
#' @param tree_file path to a .tre file
#' @param trim_pattern optional pattern to trim from tip labels of the tree
#' @param ... other arguments to pass to brms::brm
#' @details the tip labels of the tree must be the sample ids from the metadata.
#'   You can use the \code{trim_pattern} argument to automatically trim off any
#'   consistent pattern from the tip labels if necessary. The dots can be used
#'   to pass cores=4 to brms to make the chains run in parallel.
#'
#'   The prior for the intercept is a normal distribution centered on the mean
#'   of the outcome variable with a standard deviation of 3*sd(outcome variable)
#' @inheritParams anpan
#' @export
anpan_pglmm = function(meta_file,
                       tree_file,
                       trim_pattern = NULL,
                       covariates = NULL,
                       plot_cov_mat = TRUE,
                       outcome = "age",
                       verbose = TRUE,
                       ...) {
  meta = read_meta(meta_file,
                   select_cols = c("sample_id", covariates, outcome))

  bug_tree = ape::read.tree(tree_file)

  bug_tree$tip.label = gsub("_bowtie2", "", bug_tree$tip.label)
  cov_mat = ape::vcv.phylo(bug_tree)

  if (plot_cov_mat) {
    p = make_cov_mat_plot(cov_mat,
                          bug_name = basename(tree_file))
  }

  model_input = meta %>%
    dplyr::filter(sample_id %in% bug_tree$tip.label)

  if (nrow(model_input) < nrow(meta) & verbose){
    message(paste0(nrow(model_input), ' samples out of ', nrow(meta), " present in the metadata included in the analysis." ))
  }

  if (!is.null(covariates)) {
    cov_str = paste0(paste(covariates, collapse = " + "), " + ")
  } else {
    cov_str = NULL
  }

  model_formula = as.formula(paste0(outcome, " ~ ", cov_str, "(1|gr(sample_id, cov = cov_mat))"))

  model_fit = brm(formula = model_formula,
                  data = model_input,
                  family = gaussian(), #TODO alternate outcome families
                  data2 = list(cov_mat = cov_mat),
                  prior = c(
                    prior(normal(mean(model_input[[outcome]]), 3*sd(model_input[[outcome]])), "Intercept"),
                    prior(student_t(3, 0, 20), "sd"),
                    prior(student_t(3, 0, 20), "sigma")
                  ),
                  backend = "cmdstanr")

  model_fit

}
