.onAttach = function(libname, pkgname) {
  packageStartupMessage("anpan uses the furrr and progressr packages.\n * Run future::plan() (before calling anpan) in a way that's appropriate for your system to parallelize.\n * Run library(progressr); handlers(global=TRUE) to show progress.")
}
