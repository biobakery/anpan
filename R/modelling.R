# Unsure if this syntax is right but it's close I think
# brm(crc ~ 1+ age + gender + (present | uniref) + (1 | bug),
#     family = binomial(),
#     data = whole_species,
#     prior = prior(prior("normal(0,1)",
#                         coef = c("age", "gender", "bug")),
#                   prior("horseshoe(par_ratio = .005)",
#                         coef = "uniref")))


library(tidyverse); theme_set(theme_light())
library(data.table)
library(magrittr)
library(brms)
library(cmdstanr)

windoze = Sys.info()['sysname'] == "Windows"

if (windoze) {
  # plan(multisession, workers = 5)
  n_core = 4
  file_dir = "../../../../strain_stats/data/strains_binary_genes/"
} else {
  # plan(multisession, workers = 12)
  n_core = 4
  file_dir = 'data/strains_binary_genes'
}

meta = fread('data/CRC_analysis_metadata_final_version.tsv')
meta$crc = c(CRC = TRUE, control = FALSE)[meta$study_condition]
meta_cov = meta %>% select(sampleID, age, gender, crc) %>% as.data.table

read_bug = function(bug_file) {
  nc = readLines(bug_file,
                 n = 1) %>%
    strsplit('\t') %>%
    .[[1]] %>%
    length

  gf = fread(bug_file,
             colClasses = list(character = 1, numeric = 2:nc)) %>%
    dplyr::select_all(~gsub("_Abundance-RPKs", "", .))

  names(gf)[1] = "gene_spec"
  gf = gf %>% select(gene_spec, any_of(unique(meta$sampleID)))

  melt(gf, id.vars = "gene_spec", variable.name = 'sampleID', value.name = "bug")
}


bug_files = list.files(file_dir, full.names = TRUE, pattern = "*.binary.tsv")

minmax_thresh = 5

# priors by predictor, one bug at a time ----------------------------------
gf = read_bug(bug_files[4]) %>%
  .[,.(gene_spec, sampleID, bug,
       varies_enough = sum(bug) < (.N - minmax_thresh) &
         sum(bug) > minmax_thresh), by = gene_spec] %>%
  .[(varies_enough)] %>%
  .[,.(gene_spec, sampleID, bug)]

joined = gf[meta_cov, on = 'sampleID', nomatch = 0][,.(gene_spec, present = bug, sampleID, age, gender, crc)] %>%
  tidyr::separate(gene_spec,
                  into = c('uniref', "species"),
                  sep = "\\|") %>%
  .[,.(n_crc = sum(crc),
       tot = .N), by = .(uniref, species, present, age, gender, sampleID)]

wide_dat = dcast(joined[,!"species"],
                 age + gender + sampleID + n_crc + tot ~ uniref,
                 value.var = 'present')
wide_subset = wide_dat[,1:50]

fit_wide = function(wide) {
  form = bf(n_crc ~ stdCovariates + unirefs, nl = TRUE) +
    lf(stdCovariates  ~ 1 + age + gender , center = TRUE) +
    lf(as.formula(paste0("unirefs ~ 0 + ",
                         paste(names(wide)[-(1:5)],
                               collapse = " + "))), cmc = FALSE)

  p = set_prior("horseshoe(par_ratio = .005)",
                class = 'b', nlpar = "unirefs") +
    prior(normal(0,2),
          class = "b", nlpar = 'stdCovariates')

  brm(form,
      prior = p,
      data = wide,
      family = bernoulli(),
      refresh = 50,
      control = list(adapt_delta = .9),
      backend = 'cmdstanr',
      chains = 4,
      cores = 4,
      threads = threading(4))

}

sub_fit = fit_wide(wide_subset)
whole_fit = fit_wide(wide_dat)

post_sum = posterior_summary(whole_fit)

post_df = post_sum %>%
  as.data.frame %>%
  rownames_to_column('param') %>%
  as.data.table()

save(whole_fit,
     post_df,
     file = 'outputs/A_muciniphilia_fit.RData')
