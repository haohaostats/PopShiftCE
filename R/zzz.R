
NULL

.onLoad <- function(libname, pkgname) {
  # Reasonable, reproducible defaults; users can override
  try({
    set.seed(12345)
    RNGkind("L'Ecuyer-CMRG")
  }, silent = TRUE)
}
