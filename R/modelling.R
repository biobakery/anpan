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
#   .[,.(gene, sampleID, bug,
#        varies_enough = sum(bug) < (.N - minmax_thresh) &
#          sum(bug) > minmax_thresh), by = gene] %>%
#   .[(varies_enough)] %>%
#   .[,.(gene, sampleID, bug)]
#
# joined = gf[meta_cov, on = 'sampleID', nomatch = 0][,.(gene, present = bug, sampleID, age, gender, crc)] %>%
#   tidyr::separate(gene,
#                   into = c('uniref', "species"),
#                   sep = "\\|") %>%
#   .[,.(n_crc = sum(crc),
#        tot = .N), by = .(uniref, species, present, age, gender, sampleID)]
#
# wide_dat = dcast(joined[,!"species"],
#                  age + gender + sampleID + n_crc + tot ~ uniref,
#                  value.var = 'present')
# wide_subset = wide_dat[,1:50]

get_bug_name = function(bug_file,
                        remove_pattern = ".genefamilies.tsv") {
  gsub(remove_pattern, "", basename(bug_file))
}


fit_glms = function(model_input, out_dir) {
  glm_fits = model_input[,.(data_subset = list(.SD)), by = gene]

  # Progress won't be that hard: https://furrr.futureverse.org/articles/articles/progress.html#package-developers-1
  # p <- progressor(steps = dplyr::n_distinct(model_input$gene))
  # ^ That goes right here, then activate the p() calls commented out in fit_glm()
  glm_fits$glm_res = furrr::future_map(glm_fits$data_subset,
                                       safely_fit_glm) # TODO progress bar with progressr


  failed = glm_fits[sapply(glm_fits$glm_res,
                           \(.x) !is.null(.x$error))]
  if (nrow(failed) > 0){
    # TODO Write out the failures to a warning file with a message
  }

  worked = glm_fits[sapply(glm_fits$glm_res,
                           \(.x) is.null(.x$error))]

  all_terms = tidyr::unnest(worked[,.(gene, glm_res = lapply(glm_res, \(.x) .x$result))],
                            cols = c(glm_res))

  readr::write_tsv(all_terms,
                   file = paste0(out_dir, "all_terms.tsv"))

  bug_terms = all_terms %>%
    dplyr::filter(term == "presentTRUE") %>%
    dplyr::arrange(p.value) %>%
    dplyr::mutate(q = p.adjust(p.value, method = 'fdr')) %>%
    dplyr::select(-term)

  readr::write_tsv(bug_terms,
                   file = paste0(out_dir, "bug_terms.tsv"))

  return(bug_terms)
}

check_prevalence_okay = function(gene_dat, prevalence_filter) {

  min_prev_by_outcome = table(gene_dat[,.(present, crc)]) %>%
    apply(2, \(.x) .x / sum(.x)) %>%
    apply(2, min) # It has to be variable above the prevalence filter in BOTH groups

  if (n_distinct(gene_dat$crc) != 2 |                         # If only one outcome shows up
      n_distinct(gene_dat$present) == 1 |                     # or if the gene is omni-present or omni-absent
      all(prev_by_outcome < prevalence_filter)) {             # or if it doesn't get through the prevalence filter
    return(FALSE) # Then it's not okay
  } else {
    return(TRUE)
  }
}

#' Fit a GLM to one gene
fit_glm = function(gene_dat, out_dir, prevalence_filter = .05) {

  if (!check_prevalence_okay(gene_dat, prevalence_filter)) {
    # p()
    return(data.table(term = character(),
                      estimate = numeric(),
                      std.error = numeric(),
                      statistic = numeric(),
                      p.value = numeric()))
  }

  res = glm(crc ~ age + gender + present, # TODO adjustable covariates
            data = gene_dat,
            family = 'binomial') %>%
    broom::tidy() %>%
    as.data.table()
  # p()
  return(res)
}

safely_fit_glm = purrr::safely(fit_glm)

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
#' @param save_filter_stats logical indicating whether to save filter statistics
#' @export
scone = function(bug_file,
                 meta_file,
                 out_dir,
                 model_type = "scone",
                 covariates = c("age", "sex"),
                 save_filter_stats = FALSE) {


# Checks ------------------------------------------------------------------
  message("Preparing the mise en place (checking inputs)...")

  if (!(model_type %in% c("scone", "glm"))) stop('model_type must be either "scone" or "glm"')

  bug_name = get_bug_name(bug_file)
  n_lines = R.utils::countLines(bug_file)

  if (!grepl('/$', out_dir)){
    out_dir = paste0(out_dir, "/")
  }

  if (!dir.exists(out_dir)){
    message("* Output directory doesn't exist. Creating it.")
    dir.create(out_dir)
  }

  if (save_filter_stats) {
    message("* Creating the filter stats directory in the output directory.")
    filter_stats_dir = paste0(out_dir, "filter_stats/")
    dir.create(filter_stats_dir)
  }


# Filtering ---------------------------------------------------------------

  message("Reading and filtering the data.")
  model_input = read_and_filter(bug_file, read_meta(meta_file),
                                pivot_wide = model_type == "scone",
                                save_filter_stats = save_filter_stats)


# Fitting -----------------------------------------------------------------


  res = switch(model_type,
               scone = fit_scone(model_input, out_dir),
               glm = fit_glms(model_input, out_dir))


# Summarizing -------------------------------------------------------------


# Write output ------------------------------------------------------------


  return(res)
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
