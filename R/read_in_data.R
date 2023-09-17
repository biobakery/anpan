read_meta = function(meta_file,
                     select_cols = c("sample_id", "age", "gender", "crc"),
                     omit_na = FALSE) {

  if (is.character(meta_file) && file.exists(meta_file)) {
    meta = fread(meta_file,
                 showProgress = FALSE,
                 header = TRUE)
  } else if (!is.data.frame(meta_file)) {
    stop("meta_file doesn't seem to be a path to a file nor a data frame.")
  } else {
    meta = as.data.table(meta_file)
  }

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

  metadata = meta |> dplyr::select(dplyr::all_of(select_cols))

  if (omit_na) {
    metadata = na.omit(metadata)
  } else if (any(!complete.cases(metadata))) {
    stop("Incomplete cases detected in metadata. Set omit_na = TRUE to omit NAs.")
  }

  metadata
}

#' Read a bug's genefamily file
#' @param bug_file path to a bug's genefamily file
#' @param meta a data frame of metadata
#' @param remove_pattern pattern to remove from the column names of the
#'   genefamily file
#' @details
#' The input bug_file needs to be readable by data.table::fread()
#'
#' If metadata is provided, the genefamily file is subset to only those samples
#' present in the sample_id column of the metadata.
#' @returns The genefamily file of the bug as a data.table in TALL format
#' @export
read_bug = function(bug_file, meta = NULL,
                    remove_pattern = "_Abundance-RPKs") {

  gf = fread(bug_file,
             showProgress = FALSE,
             header = TRUE) |>
    dplyr::select_all(~gsub(remove_pattern, "", .)) # This is why we need to import dplyr

  names(gf)[1] = "gene"

  col_classes = gf |> sapply(class)

  if (any(col_classes[-1] != "numeric")) {
    for (i in (which(col_classes[-1] != "numeric") + 1)) {
      gf[[i]] = as.numeric(gf[[i]])
    }
  }

  gf$gene = gsub("\\|(.*)", "", gf$gene)
  # ^ This removes the |species_id part of the identifier to make it easier to read

  if (!is.null(meta)){
    gf = gf |> dplyr::select(gene, any_of(unique(meta$sample_id)))
  }

  melt(gf, id.vars = "gene", variable.name = 'sample_id', value.name = "abd")
}

get_file_list = function(file_dir) {
  list.files(file_dir, full.names = TRUE)
}

write_tsv_no_progress = function(x, file, append = FALSE) {

  data.table::fwrite(x = x,
                     file = file,
                     append = append,
                     sep = "\t",
                     showProgress = FALSE)

  return(invisible(x))
}


