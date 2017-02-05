#' Fast S-Transform
#'
#' \code{fst} is used to compute the fast S-transform of a real,
#' univariate time series.
#'
#' @param xt Equally-sampled input time series. Must convert to numeric vector.
#' @param dt Sample interval, in seconds. Default is 0.01.
#' @param gauss TRUE for Gaussian windows, FALSE for Box windows. Default is TRUE.
#' @param img.dim Image dimension of transform, in pixels. Default is no downsampling.
#' @param by.part Image is by time and frequency partition rather than by the implicit
#' N time and frequency values. Default is TRUE
#' @param all.freqs Image contains negative as well as non-negative frequencies.
#' Default is FALSE.
#' @param ind.map Image contains the index of the S-Transform value, rather than
#' the value. Useful for understanding partitions. Default is FALSE.
#' @return A list including the complex S-transform of the data, an image sampled on
#' the time and freq partition centers (usually small in dimension), the height and width
#' of the image, the frequency-partition centers, the time-partition centers,
#' the window type used, and other info.
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
fst <- function(xt, dt=0.01, gauss=TRUE, img.dim=NA, by.part=TRUE,
                all.freqs=FALSE, ind.map=FALSE) {

	xt <- as.double(xt[!is.na(xt)])
	len <- length(xt)
	if ( len < 3 )
		stop("input time series must have at least 3 valid points")

	# call C function
	out <- .Call("CALLngft_1dComplex64",
							 as.double(xt), as.double(dt), as.logical(gauss), as.integer(img.dim),
							 as.logical(by.part), as.logical(all.freqs), as.logical(ind.map))
	return(out)
}
