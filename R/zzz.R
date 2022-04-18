.onAttach = function(libname, pkgname) {
  ver = utils::packageVersion("anpan")
  packageStartupMessage("This is anpan version ", ver)
  packageStartupMessage("- Get help: Visit the biobakery help forum at https://forum.biobakery.org/")
  packageStartupMessage("- Parallelize: Before calling anpan, run future::plan() in a way that's appropriate for your system.")
  packageStartupMessage("- Show progress: Before calling anpan, run library(progressr); handlers(global=TRUE)")
}
