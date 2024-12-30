.onAttach = function(libname, pkgname) {
  ver = utils::packageVersion("anpan")

  vignette_built = nrow(vignette(package = 'anpan')$results) != 0

  if (vignette_built) {
    rebuild_advice = ""
  } else {
    rebuild_advice = "Reinstall with {.code build_vignettes = TRUE}, then "
  }

  vign_emph = cli::combine_ansi_styles("bold", "red2")

  vignette_msg = paste0('{cli::symbol$bullet} {.strong Read the guide}: ',
                        rebuild_advice,
                        'run {vign_emph("anpan::anpan_vignette()")}')

  cli::cli_inform(paste0("{cli::symbol$bullet} This is anpan version {.strong ", ver, "}"), class = "packageStartupMessage")
  cli::cli_inform(vignette_msg, class = "packageStartupMessage")
  cli::cli_inform("{cli::symbol$bullet} {.strong Get help}: Visit the biobakery help forum at {.url https://forum.biobakery.org/}", class = "packageStartupMessage")
  cli::cli_inform("{cli::symbol$bullet} {.strong Parallelize}: Run {.fn future::plan} as appropriate for your system.", class = "packageStartupMessage")
  cli::cli_inform("{cli::symbol$bullet} {.strong Activate progress bars}: {.code library(progressr); handlers(global=TRUE)}", class = "packageStartupMessage")
}

#' Open the anpan vignette
#' @description This function checks if you have the vignette built, then opens
#' it with utils::RShowDoc() if you do. This function exists just to make it a
#' little easier for users to find their way to the vignette.
#'
#' @export
anpan_vignette = function() {
  vignette_built = nrow(vignette(package = 'anpan')$results) != 0

  if (vignette_built) {
    message("Opening the vignette in your browser...")
    utils::RShowDoc("anpan_tutorial", type = "html", package = "anpan")
  } else {
    stop("You don't have the vignette built! Reinstall with:\n\nremotes::install_github('biobakery/anpan', force = TRUE, build_vignettes = TRUE)\n\nthen try again!")
  }
}
