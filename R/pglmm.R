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

anpan_pglmm = function(meta_path,
                       tree_path) {
  meta = read_meta(meta_path,
                   select_cols = c("sample_id", ))

}
