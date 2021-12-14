.onLoad = function(libname, pkgname) {
  message("This package uses furrr. Run future::plan() in a way that's appropriate for your system to parallelize.")
}
