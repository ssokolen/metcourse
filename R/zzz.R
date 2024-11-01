.onLoad <- function(libname, pkgname) {
  if (requireNamespace("rstan", quietly = TRUE)) {
    rstan::rstan_options(auto_write = TRUE)
    options(mc.cores = min(parallel::detectCores(), 1))
  } else {
    warning("rstan not found, some functionality will not work.")
  }
}
