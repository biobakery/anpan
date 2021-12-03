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


fit_glm = function(wide, bug_file){

}

fit_scone = function(wide, bug_name, tpc = 4, ncore = 4, out_path,
                    save_summ = FALSE) {
  form = brms::bf(n_crc ~ stdCovariates + unirefs, nl = TRUE) +
    brms::lf(stdCovariates  ~ 1 + age + gender , center = TRUE) +
    brms::lf(as.formula(paste0("unirefs ~ 0 + ",
                               paste(names(wide)[-(1:5)],
                                     collapse = " + "))), cmc = FALSE)

  p = brms::set_prior("horseshoe(par_ratio = .005)",
                      class = 'b', nlpar = "unirefs") +
    brms::prior(normal(0,2),
                class = "b", nlpar = 'stdCovariates')

  model_fit = brms::brm(form,
            prior = p,
            data = wide,
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
                 model_type = "scone",
                 covariates = c("age", "gender")) {

  if (!(model_type %in% c("scone", "glm"))) stop('model_type must be either "scone" or "glm"')

  bug_name = gsub(".genefamilies.tsv", "", basename(gf_file))
  n_lines = R.utils::countLines(gf_file)

  #
  input_data = read_and_filter(bug_file, read_meta(meta_file))

  res = switch(model_type,
               scone = fit_scone(input_data),
               glm = fit_glm(input_data))
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
