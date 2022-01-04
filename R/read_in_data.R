read_meta = function(meta_file,
                     select_cols = c("sample_id", "age", "gender", "crc")) {

  meta = fread(meta_file)

  # Handle missing or alternative sample_id column names
  if (!("sample_id" %in% names(meta))) {
    if (sum(c("sampleID", "SampleID", "sampID", "samp_id", "sample_ID")  %in% names(meta)) == 1) {
      sid_i = which(names(meta) %in% c("sampleID", "SampleID", "sampID", "samp_id", "sample_ID"))
      names(meta)[sid_i] = 'sample_id'
    } else {
      stop("Couldn't find the sample_id column in the metadata file.")
    }
  }

  meta_cov = meta %>% dplyr::select(dplyr::all_of(select_cols))
  meta_cov
}

read_bug = function(bug_file, meta = NULL,
                    remove_pattern = "_Abundance-RPKs") {
  nc = readLines(bug_file,
                 n = 1) %>%
    strsplit('\t') %>%
    .[[1]] %>%
    length

  gf = fread(bug_file,
             colClasses = list(character = 1, numeric = 2:nc)) %>%
    dplyr::select_all(~gsub(remove_pattern, "", .)) # This is why we need to import dplyr

  names(gf)[1] = "gene"

  gf$gene = gsub("\\|(.*)", "", gf$gene)
  # ^ This removes the |species_id part of the identifier to make it easier to read

  if (!is.null(meta)){
    gf = gf %>% dplyr::select(gene, any_of(unique(meta$sample_id)))
  }

  melt(gf, id.vars = "gene", variable.name = 'sample_id', value.name = "abd")
}

get_file_list = function(file_dir) {
  list.files(file_dir, full.names = TRUE, pattern = "*.genefamilies.tsv")
}


