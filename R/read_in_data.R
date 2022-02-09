read_meta = function(meta_file,
                     select_cols = c("sample_id", "age", "gender", "crc")) {

  meta = fread(meta_file,
               showProgress = FALSE)

  # Handle missing or alternative sample_id column names
  if (!("sample_id" %in% names(meta))) {
    alternate_name = grepl(pattern = "samp[l]?[e]?[:punct:]?[i]?[d]?",
                           x = names(meta),
                           ignore.case = TRUE)
    unique_name_found = sum(alternate_name == 1)
    if (unique_name_found) {
      sid_i = which(alternate_name)
      names(meta)[sid_i] = 'sample_id'
    } else {
      stop("Couldn't find a unique sample_id column in the metadata file.")
    }
  }

  metadata = meta %>% dplyr::select(dplyr::all_of(select_cols))
  metadata
}

#' @export
read_bug = function(bug_file, meta = NULL,
                    remove_pattern = "_Abundance-RPKs") {
  nc = readLines(bug_file,
                 n = 1) %>%
    strsplit('\t') %>%
    .[[1]] %>%
    length

  gf = fread(bug_file,
             colClasses = list(character = 1, numeric = 2:nc),
             showProgress = FALSE) %>%
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

# The readr progress messes up the progressr bar
write_tsv_no_progress = function(x, file) {

  old = options("readr.show_progress")$readr.show_progress

  options(readr.show_progress = FALSE)

  readr::write_tsv(x = x,
                   file = file)

  options(readr.show_progress = old)

  return(invisible(x))
}


