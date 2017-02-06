# these functions should not be exported

# .onLoad is called by loadNamespace, usually via library()
.onLoad <- function(libname, pkgname) {
  # message(".onLoad: libname=",libname,", pkgname=", pkgname)
  # don't load ngft if useDynLib directive is used in NAMESPACE file
  # library.dynam("ngft", pkgname, libname)
  library.dynam("libfftw3-3", pkgname, libname)
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
