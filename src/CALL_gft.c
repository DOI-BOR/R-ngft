#include "gft.h"
#include "gft_proto.h"

/* be sure to include these last so they can override some windows defs */
#include <R.h>
#include <Rinternals.h>

DllExport SEXP CALLgft_1dComplex64(SEXP ts_d, SEXP dt_d, SEXP gauss_i)
{
	int ii, pcnt = 0;
	Rcomplex *cts, *ctp;
	SEXP cts_d;
	unsigned int n_samples;
	int stride, p_size;
	int *partitions;
	double *window_sets;
	windowFunction *wf;

	double *ts = REAL(ts_d);
	int len = length(ts_d);
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];
	BOOL gaussian_window = asLogical(gauss_i) == NA_LOGICAL ? TRUE : LOGICAL(gauss_i)[0];

	if ( len < 3 )
		error("time series must have at least 3 points");
	if ( dt <= 0 )
		error("dt must be positive");

	// Utility Functions
	DllExport void gft_1d_shift(double *signal, unsigned int N, unsigned int shiftBy);

	// Interpolation functions
	DllExport double *gft_1d_interpolateNN(double *signal, unsigned int N, unsigned int M);

	/* initialize */
	n_samples = (unsigned)len;
	partitions = gft_1dPartitions(n_samples);
	wf = gaussian_window ? gaussian : box;
	window_sets = windowsFromPars(n_samples, wf, partitions);
	p_size = gft_1dSizeOfPartitions(n_samples);
	stride = 1;

	if ( (cts = calloc( len, sizeof( *cts ) )) == NULL )
		error( "CALLgft_1dComplex64: can't get space" );
	for ( ii = 0 ; ii < len ; ii++ )
		cts[ii].r = ts[ii];

	// Call 1D GFT Function
	gft_1dComplex64((double *)cts, n_samples, window_sets, partitions, stride );

	/* allocate space for R structures for complex transform, and copy gft output */
	cts_d = PROTECT(allocVector(CPLXSXP, len)); pcnt++;
	ctp = COMPLEX(cts_d);
	for ( ii = 0 ; ii < len ; ii++ )
		ctp[ii] = cts[ii];

	free(cts);
	free( partitions );
	free( window_sets );

	UNPROTECT(pcnt);

	return cts_d;
}
