#' Fast Inverse S-Transform
#'
#' \code{ifst} is used to compute the fast inverse S-transform of a complex,
#' univariate S-transform.
#'
#' @param fst Fast S-transform. Must convert to complex vector.
#' @param dt Sample interval, in seconds. Default is 0.01.
#' @param gauss TRUE for Gaussian windows, FALSE for Box windows. Default is TRUE.
#' @return A list including the complex time series of the data, and
#' the window type used.
#' @details If Gaussian windows are selected, the results are discrete S-transforms.
#' If Box windows are selected, the results are discrete orthonormal S-transforms,
#' with critical sampling in both frequency and time/space dimensions.
#'
#' This function is a wrapper for the Generalized Fourier Family Transform (GFT)
#' library \url{https://sourceforge.net/projects/fst-uofc/}.
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
#' \item \href{https://sourceforge.net/projects/fst-uofc/}{GFT library at SourceForge}
#' }
#' @keywords ts
ifst <- function(dst, dt=0.01, gauss=NA) {

  xt <- as.double(xt[!is.na(xt)])
  len <- length(xt)
  if ( len < 3 )
    stop("input time series must have at least 3 valid points")

  # call C function
  out <- .Call("CALLngft_1dComplex64Inv",
               as.complex(dst), as.double(dt), as.logical(gauss))
  return(out)
}
