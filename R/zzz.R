.onLoad <- function(libname, pkgname) {
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)  
}
