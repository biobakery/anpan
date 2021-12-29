.onAttach = function(libname, pkgname) {
  message("anpan uses the furrr and progressr packages.\n * Run future::plan() (before calling anpan) in a way that's appropriate for your system to parallelize.\n * Run library(progressr); handlers(global=TRUE) to show progress.")
}
