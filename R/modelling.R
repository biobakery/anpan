# Unsure if this syntax is right but it's close I think
# brm(crc ~ 1+ age + gender + (present | uniref) + (1 | bug),
#     family = binomial(),
#     data = whole_species,
#     prior = prior(prior("normal(0,1)",
#                         coef = c("age", "gender", "bug")),
#                   prior("horseshoe(par_ratio = .005)",
#                         coef = "uniref")))


# library(tidyverse); theme_set(theme_light())
# library(data.table)
# library(magrittr)
# library(brms)
# library(cmdstanr)

# windoze = Sys.info()['sysname'] == "Windows"
#
# if (windoze) {
#   # plan(multisession, workers = 5)
#   n_core = 4
#   file_dir = "../../../../strain_stats/data/strains_binary_genes/"
# } else {
#   # plan(multisession, workers = 12)
#   n_core = 4
#   file_dir = 'data/strains_binary_genes'
# }



# minmax_thresh = 5
#
# # priors by predictor, one bug at a time ----------------------------------
# gf = read_bug(bug_files[4]) %>%
#   .[,.(gene_spec, sampleID, bug,
#        varies_enough = sum(bug) < (.N - minmax_thresh) &
#          sum(bug) > minmax_thresh), by = gene_spec] %>%
#   .[(varies_enough)] %>%
#   .[,.(gene_spec, sampleID, bug)]
#
# joined = gf[meta_cov, on = 'sampleID', nomatch = 0][,.(gene_spec, present = bug, sampleID, age, gender, crc)] %>%
#   tidyr::separate(gene_spec,
#                   into = c('uniref', "species"),
#                   sep = "\\|") %>%
#   .[,.(n_crc = sum(crc),
#        tot = .N), by = .(uniref, species, present, age, gender, sampleID)]
#
# wide_dat = dcast(joined[,!"species"],
#                  age + gender + sampleID + n_crc + tot ~ uniref,
#                  value.var = 'present')
# wide_subset = wide_dat[,1:50]


fit_glms = function(model_input, out_dir, non_element_cols = 1:4) {
  # non_element_cols = indices of columns that aren't genes.
  # Select covariates, outcome, and iteratively each additional gene-family column. GLM them each in turn

  element_cols = (1:ncol(model_input))[-non_element_cols] # This is brittle



}

fit_glm = function(model_input, out_dir) {
  # Select covariates, outcome, and iteratively each additional gene-family column. GLM them each in turn
  if (n_distinct(gene_spec_dat$crc) != 2 |
      n_distinct(gene_spec_dat$bug) == 1 |
      min(table(gene_spec_dat$bug) / nrow(gene_spec_dat)) < .05) {
    return(data.table(term = character(),
                      estimate = numeric(),
                      std.error = numeric(),
                      statistic = numeric(),
                      p.value = numeric(),
                      gs = character()))
  }

  glm(crc ~ age + gender + bug,
      data = gene_spec_dat,
      family = 'binomial') %>%
    broom::tidy() %>%
    as.data.table() %>%
    .[, gs := gs]
}

fit_scone = function(model_input, bug_name, tpc = 4, ncore = 4, out_path,
                    save_summ = FALSE) {
  form = brms::bf(n_crc ~ stdCovariates + unirefs, nl = TRUE) +
    brms::lf(stdCovariates  ~ 1 + age + gender , center = TRUE) +
    brms::lf(as.formula(paste0("unirefs ~ 0 + ",
                               paste(names(model_input)[-(1:5)],
                                     collapse = " + "))), cmc = FALSE)

  p = brms::set_prior("horseshoe(par_ratio = .005)",
                      class = 'b', nlpar = "unirefs") +
    brms::prior(normal(0,2),
                class = "b", nlpar = 'stdCovariates')

  model_fit = brms::brm(form,
            prior = p,
            data = model_input,
            family = brms::bernoulli(),
            refresh = 50,
            control = list(adapt_delta = .9),
            backend = 'cmdstanr',
            chains = 4,
            cores = ncore,
            threads = brms::threading(tpc))

  summ = summarise_fit(model_fit)

  if (save_summ){
    save(summ,
         file = paste0(out_path, "/", bug_name, ".RData")) # TODO make the parts of this work
  }

  model_fit

}

#' Run scone
#'
#' @description Run the scone model on a
#' @param bug_file path to a gene family file (usually probably from HUMAnN)
#' @param meta_file path to a metadata tsv. Must contain the specify covariates
#' @param model_type either "scone" or "glm"
#' @param covariates covariates to account for
#' @export
scone = function(bug_file,
                 meta_file,
                 out_dir,
                 model_type = "scone",
                 covariates = c("age", "sex")) {

  if (!(model_type %in% c("scone", "glm"))) stop('model_type must be either "scone" or "glm"')

  bug_name = gsub(".genefamilies.tsv", "", basename(bug_file))
  n_lines = R.utils::countLines(bug_file)

  #
  model_input = read_and_filter(bug_file, read_meta(meta_file),
                                pivot_wide = model_type == "scone")

  res = switch(model_type,
               scone = fit_scone(model_input),
               glm = fit_glms(model_input))
}

#' @export
summarise_scone = function(model_fit) {
  post_sum = brms::posterior_summary(model_fit)

  post_df = post_sum %>%
    as.data.frame %>%
    tibble::rownames_to_column('param') %>%
    as.data.table()

  post_df
}
