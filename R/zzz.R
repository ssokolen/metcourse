.onLoad <- function(libname, pkgname) {
  if (requireNamespace("rstan", quietly = TRUE)) {
    rstan::rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
  } else {
    warning("rstan not found, some functionality will not work.")
  }
}
