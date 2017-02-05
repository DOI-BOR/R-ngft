# these functions should not be exported

# .onLoad is called by loadNamespace, after loading
.onLoad <- function(libname, pkgname) {
  # don't load ngft if useDynLib directive is used in NAMESPACE file
  # library.dynam("ngft", pkgname, libname)
  library.dynam("libfftw3-3", pkgname, libname)
}

# .onAttach is called by attachNamespace (ususlly done by library)
.onAttach <- function(libname, pkgname) {
  # lib.loc <- file.path(libname,pkgname,"libs",.Platform$r_arch)
  # message(".onLoad: libname=",libname,", pkgname=", pkgname)
}

# .onUnLoad is called by unloadNamespace, after unloading
.onUnload <- function(libpath) {
  library.dynam.unload("ngft", libpath)
  library.dynam.unload("libfftw3-3", libpath)
}

# .onDetach is called by detach, or unloadNamespace if the namespace is attached
.onDetach <- function(libpath) {
}
