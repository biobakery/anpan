read_meta = function(meta_path,
                     select_cols = c("sampleID", "age", "gender", "crc")) {
  meta = fread(meta_path) # 'data/CRC_analysis_metadata_final_version.tsv'
  meta$crc = c(CRC = TRUE, control = FALSE)[meta$study_condition]
  meta_cov = meta[, ..select_cols]
  meta_cov
}

read_bug = function(bug_file, meta,
                    remove_pattern = "_Abundance-RPKs") {
  nc = readLines(bug_file,
                 n = 1) %>%
    strsplit('\t') %>%
    .[[1]] %>%
    length

  gf = fread(bug_file,
             colClasses = list(character = 1, numeric = 2:nc)) %>%
    dplyr::select_all(~gsub(remove_pattern, "", .)) # This is why we need dplyr

  names(gf)[1] = "gene"

  gf$gene = gsub("\\|(.*)", "", gf$gene)
  # ^ This removes the |species_id part of the identifier to make it easier to read

  gf = gf %>% dplyr::select(gene, any_of(unique(meta$sampleID)))

  melt(gf, id.vars = "gene", variable.name = 'sampleID', value.name = "abd")
}

get_file_list = function(file_dir) {
  list.files(file_dir, full.names = TRUE, pattern = "*.binary.tsv")
}


