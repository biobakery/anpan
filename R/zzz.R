.onAttach = function(libname, pkgname) {
  ver = utils::packageVersion("anpan")

  vignette_built = nrow(vignette(package = 'anpan')$results) != 0

  if (vignette_built) {
    rebuild_advice = ""
  } else {
    rebuild_advice = "Reinstall with build_vignettes = TRUE, then "
  }

  vignette_msg = paste0('- Read the vignette: ',
                        rebuild_advice,
                        'anpan::anpan_vignette()')

  packageStartupMessage("This is anpan version ", ver)
  packageStartupMessage(vignette_msg)
  packageStartupMessage("- Get help: Visit the biobakery help forum at https://forum.biobakery.org/")
  packageStartupMessage("- Parallelize: Before calling anpan, run future::plan() in a way that's appropriate for your system.")
  packageStartupMessage("- Show progress: Before calling anpan, run library(progressr); handlers(global=TRUE)")
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
