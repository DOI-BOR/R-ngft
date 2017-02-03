#include "ngft.h"
#include "ngft_proto.h"

/* be sure to include these last so they can override some windows defs */
#include <R.h>
#include <Rinternals.h>

DllExport SEXP CALLngft_1dComplex64(SEXP ts_d, SEXP dt_d, SEXP gauss_i, SEXP image_dim_i)
{
	int ii, pcnt = 0;
	DCMPLX *cts;
	DIMAGE *image;
	Rcomplex *ctp;
	double *dp;
	SEXP cts_st;
	SEXP cts_img;
	SEXP img_wd;
	SEXP img_ht;
	SEXP freq_c;
	SEXP time_c;
	SEXP win_name;
	SEXP ret_l, names_s;
	int stride;
	windowFunction *window_fn;
	FPCOL *partitions;
	TPCOL *tpcol;
	ILIST *f_centers;
	ILIST *t_centers;

	double *ts = REAL(ts_d);
	int n_samples = length(ts_d);
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];
	BOOL gaussian_window = asLogical(gauss_i) == NA_LOGICAL ? TRUE : LOGICAL(gauss_i)[0];
	int image_dim = asInteger(image_dim_i) == NA_INTEGER ? -1 : INTEGER(image_dim_i)[0];

	if ( n_samples < 3 )
		error("time series must have at least 3 points");
	if ( dt <= 0 )
		error("dt must be positive");

	/* initialize */
	window_fn = gaussian_window ? gaussian : box;
	partitions = ngft_DyadicPartitions(n_samples);
	ngft_AddWindowsToParts(partitions, window_fn);
	stride = 1;

	if ( (cts = calloc( n_samples, sizeof( *cts ) )) == NULL )
		oops( "CALLngft_1dComplex64","can't get space" );
	for ( ii = 0 ; ii < n_samples ; ii++ )
		cts[ii].r = ts[ii];

	// Call 1D GFT Function
	ngft_1dComplex64(cts, n_samples, partitions, stride);

	// get complex image
	tpcol = ngft_TimePartitions( partitions );
	image = ngft_1d_logfInterpolateNN(cts, partitions, tpcol, image_dim);

	// get time and frequency centers
	f_centers = getFreqCenters(partitions);
	t_centers = getTimeCenters(tpcol);

	/* allocate space for R structures for complex transform, and copy gft output */
	cts_st = PROTECT(allocVector(CPLXSXP, n_samples)); pcnt++;
	ctp = COMPLEX(cts_st);
	for ( ii = 0 ; ii < n_samples ; ii++ ) {
		ctp[ii].r = cts[ii].r;
		ctp[ii].i = cts[ii].i;
	}
	cts_img = PROTECT(allocVector(CPLXSXP, image->img->count)); pcnt++;
	ctp = COMPLEX(cts_img);
	for ( ii = 0 ; ii < image->img->count ; ii++ ) {
		ctp[ii].r = image->img->values[ii].r;
		ctp[ii].i = image->img->values[ii].i;
	}
	img_wd = PROTECT(allocVector(INTSXP, 1)); pcnt++; INTEGER(img_wd)[0] = image->wd;
	img_ht = PROTECT(allocVector(INTSXP, 1)); pcnt++; INTEGER(img_ht)[0] = image->ht;
	freq_c = PROTECT(allocVector(REALSXP, f_centers->count)); pcnt++;
	dp = REAL( freq_c );
	for ( ii = 0 ; ii < f_centers->count ; ii++ )
		dp[ii] = f_centers->values[ii];
	time_c = PROTECT(allocVector(REALSXP, t_centers->count)); pcnt++;
	dp = REAL( time_c );
	for ( ii = 0 ; ii < t_centers->count ; ii++ )
		dp[ii] = t_centers->values[ii];
	win_name = PROTECT(allocVector(STRSXP, 1)); pcnt++;
	SET_STRING_ELT(win_name, 0, mkChar(gaussian_window ? "Gaussian" : "Box"));

	/* put the return values into a list */
	ret_l = PROTECT(allocVector(VECSXP, 7)); pcnt++;
	names_s = PROTECT(allocVector(VECSXP, 7)); pcnt++;
	SET_VECTOR_ELT(ret_l, 0, cts_st); SET_VECTOR_ELT(names_s, 0, mkChar("DST"));
	SET_VECTOR_ELT(ret_l, 1, cts_img); SET_VECTOR_ELT(names_s, 1, mkChar("image"));
	SET_VECTOR_ELT(ret_l, 2, img_wd); SET_VECTOR_ELT(names_s, 2, mkChar("wd"));
	SET_VECTOR_ELT(ret_l, 3, img_ht); SET_VECTOR_ELT(names_s, 3, mkChar("ht"));
	SET_VECTOR_ELT(ret_l, 4, freq_c); SET_VECTOR_ELT(names_s, 4, mkChar("f_centers"));
	SET_VECTOR_ELT(ret_l, 5, time_c); SET_VECTOR_ELT(names_s, 5, mkChar("t_centers"));
	SET_VECTOR_ELT(ret_l, 6, win_name); SET_VECTOR_ELT(names_s, 6, mkChar("win"));
	setAttrib(ret_l, R_NamesSymbol, names_s);

	UNPROTECT(pcnt);

	free(cts);
	freeDImage(image);
	ngft_FreeFreqPartitions(partitions);
	ngft_FreeTimePartitions(tpcol);
	freeIlist(f_centers);
	freeIlist(t_centers);

	return ret_l;
}


DllExport SEXP CALLngft_1dComplex64Inv(SEXP dst_c, SEXP dt_d, SEXP gauss_i)
{
	int ii, pcnt = 0;
	DCMPLX *cts;
	Rcomplex *ctp;
	SEXP cts_ts;
	SEXP win_name;
	SEXP ret_l, names_s;
	int stride;
	windowFunction *window_fn;
	FPCOL *partitions;

	Rcomplex *dst = COMPLEX(dst_c);
	int n_samples = length(dst_c);
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];
	BOOL gaussian_window = asLogical(gauss_i) == NA_LOGICAL ? TRUE : LOGICAL(gauss_i)[0];

	if ( n_samples < 3 )
		error("S-transform must have at least 3 points");
	if ( dt <= 0 )
		error("dt must be positive");

	/* initialize */
	window_fn = gaussian_window ? gaussian : box;
	partitions = ngft_DyadicPartitions(n_samples);
	ngft_AddWindowsToParts(partitions, window_fn);
	stride = 1;

	if ( (cts = calloc( n_samples, sizeof( *cts ) )) == NULL )
		oops( "CALLngft_1dComplex64Inv","can't get space" );
	for ( ii = 0 ; ii < n_samples ; ii++ ) {
		cts[ii].r = dst[ii].r;
		cts[ii].i = dst[ii].i;
	}

	// Call 1D GFT Inverse Function
	ngft_1dComplex64Inv(cts, partitions, stride);

	/* allocate space for R structures for complex time series, and copy inv_gft output */
	cts_ts = PROTECT(allocVector(CPLXSXP, n_samples)); pcnt++;
	ctp = COMPLEX(cts_ts);
	for ( ii = 0 ; ii < n_samples ; ii++ ) {
		ctp[ii].r = cts[ii].r;
		ctp[ii].i = cts[ii].i;
	}
	win_name = PROTECT(allocVector(STRSXP, 1)); pcnt++;
	SET_STRING_ELT(win_name, 0, mkChar(gaussian_window ? "Gaussian" : "Box"));

	/* put the return values into a list */
	ret_l = PROTECT(allocVector(VECSXP, 2)); pcnt++;
	names_s = PROTECT(allocVector(VECSXP, 2)); pcnt++;
	SET_VECTOR_ELT(ret_l, 0, cts_ts); SET_VECTOR_ELT(names_s, 0, mkChar("TS"));
	SET_VECTOR_ELT(ret_l, 1, win_name); SET_VECTOR_ELT(names_s, 1, mkChar("win"));
	setAttrib(ret_l, R_NamesSymbol, names_s);

	UNPROTECT(pcnt);

	free(cts);
	ngft_FreeFreqPartitions(partitions);

	return ret_l;
}
