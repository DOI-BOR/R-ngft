#' Fast Inverse S-Transform
#'
#' \code{ifst} is used to compute the fast inverse S-transform of a complex,
#' univariate S-transform.
#'
#' @param dst Fast S-transform, as returned by \code{\link{fst}}. Must convert to complex vector.
#' @param ts_len Length of the original time series input to \code{\link{fst}}.
#' @param dt Sample interval, in seconds. Default is 0.01.
#' @param eps Fractional window overlap, dimensionless. Default is 0.0.
#' @param part.type \describe{
#' \item{Dyadic}{\code{c("dy","Dyadic", "d")}}
#' \item{Equal Division of the Octave}{\code{c("eq","EDO", "e")}}
#' }
#' Characters are case-insensitive, and only need to uniquely
#' define the option. Default is \code{Dyadic}.
#' @param win.type \describe{
#' \item{Gaussian}{\code{c("ga","gauss")}}
#' \item{Box}{\code{c("b","box")}}
#' }
#' Characters are case-insensitive, and only need to uniquely
#' define the option. Default is \code{Gaussian}.
#' @return A list including the complex time series of the data, and
#' the window type used.
#' @details If Gaussian windows are selected, the results are discrete S-transforms.
#' If Box windows are selected, the results are discrete orthonormal S-transforms,
#' with critical sampling in both frequency and time/space dimensions.
#'
#' This function is a wrapper for the New Generalized Fourier Family Transform (NGFT)
#' library, a reimplementation of the Brown et al (2009) GFT library.
#' @seealso \itemize{
#' \item \href{https://en.wikipedia.org/wiki/S_transform}{Wikipedia}
#' \item \href{http://ieeexplore.ieee.org/document/492555/}{Stockwell et al, 1996}
#' Localization of the complex spectrum: the S transform.
#' \item \href{http://www.sciencedirect.com/science/article/pii/S1051200406000546}{Stockwell, 2007}
#' A basis for efficient representation of the S-transform.
#' \item \href{http://ieeexplore.ieee.org/document/5184926/}{Brown et al, 2009}
#' A General Description of Linear Time-Frequency Transforms
#' and Formulation of a Fast, Invertible Transform That Samples the Continuous
#' S-Transform Spectrum Nonredundantly.
#' \item \href{https://sourceforge.net/projects/fst-uofc/}{Original GFT library at SourceForge}
#' }
#' @keywords ts
ifst <- function(dst, ts_len, dt=0.01, eps=0.0, win.type=NA, part.type=NA) {

  xt <- as.double(xt[!is.na(xt)])
  len <- length(xt)
  if ( len < 3 )
    stop("input time series must have at least 3 valid points")

  # call C function
  out <- .Call("CALLngft_1dComplex64Inv",
               as.complex(gft), as.integer(ts_len), as.double(dt),
               as.double(eps), as.character(part.type), as.character(win.type))
  return(out)
}
