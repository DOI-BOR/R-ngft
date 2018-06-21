#' Time-series utilities for seismic signal processing
#'
#' ngft provides the generalized S transform.
#'
#' @details lots of details...
#'
"_PACKAGE"
#> [1] "_PACKAGE"

# these functions should not be exported

# .onLoad is called by loadNamespace, usually via library()
.onLoad <- function(libname, pkgname) {
  # message(".onLoad: libname=",libname,", pkgname=", pkgname)
  # load fft3, unless it's already loaded
  fft_dll <- "libfftw3-3"
  # library.dynam("seismic.ts", pkgname, libname)
  dlls <- getLoadedDLLs()
  dll_list <- NULL
  for ( ii in 1:length(dlls) )
    dll_list <- c(dll_list, dlls[[ii]][["name"]])
  if ( ! fft_dll %in% dll_list )
    library.dynam(fft_dll, pkgname, lib.loc=libname)
}

# .onAttach is called by attachNamespace, usually via library()
.onAttach <- function(libname, pkgname) {
  # lib.loc <- file.path(libname,pkgname,"libs",.Platform$r_arch)
  # message(".onAttach: libname=",libname,", pkgname=", pkgname)
}

# .onUnLoad is called by unloadNamespace, usually via detach()
.onUnload <- function(libpath) {
  # message(".onUnload: libpath=",libpath)
  library.dynam.unload("ngft", libpath)
  library.dynam.unload("libfftw3-3", libpath)
}

# .onDetach is called by detach
.onDetach <- function(libpath) {
  # message(".onDetach: libpath=",libpath)
}
