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
                        'browseVignettes("anpan")')

  packageStartupMessage("This is anpan version ", ver)
  packageStartupMessage(vignette_msg)
  packageStartupMessage("- Get help: Visit the biobakery help forum at https://forum.biobakery.org/")
  packageStartupMessage("- Parallelize: Before calling anpan, run future::plan() in a way that's appropriate for your system.")
  packageStartupMessage("- Show progress: Before calling anpan, run library(progressr); handlers(global=TRUE)")
}
