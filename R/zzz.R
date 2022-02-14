.onAttach = function(libname, pkgname) {
  packageStartupMessage("* Get help: Visit the biobakery help forum at https://forum.biobakery.org/\n* Parallelize: Before calling anpan, run future::plan() in a way that's appropriate for your system.\n* Show progress: Run library(progressr); handlers(global=TRUE)")
}
