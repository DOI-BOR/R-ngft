#' Fast S-Transform.
#'
#' \code{fst} is used to compute the fast S-transform of a real,
#' univariate time series.
#'
#' @param xt Equally-sampled input time series. Must convert to numeric vector.
#' @param dt Sample interval, in seconds. Default is 0.01.
#' @param gauss TRUE for Gaussian windows, FALSE for Box windows. Default is TRUE.
#' @return The complex S-transform of the data.
#' @details If Gaussian windows are selected, the results are discrete S-transforms.
#' If Box windows are selected, the results are discrete orthonormal S-transforms,
#' with critical sampling in both frequency and time/space dimensions.
#'
#' This function is a wrapper for the Generalized Fourier Family Transform (GFT)
#' library \url{https://sourceforge.net/projects/fst-uofc/}.
#' @seealso \url{https://en.wikipedia.org/wiki/S_transform}
#' \url{http://ieeexplore.ieee.org/document/492555/}
#' \url{http://www.sciencedirect.com/science/article/pii/S1051200406000546}
#' \url{http://ieeexplore.ieee.org/document/5184926/}
#' \url{https://sourceforge.net/projects/fst-uofc/}
#' @keywords ts
fst <- function(xt, dt=0.01, gauss=NA) {

	xt <- as.double(xt[!is.na(xt)])
	len <- length(xt)
	if ( len < 3 )
		stop("input time series must have at least 3 valid points")

	# call C function
	out <- .Call("CALLgft_1dComplex64",
							 as.double(xt), as.double(dt), as.logical(gauss))
	return(out)
}
