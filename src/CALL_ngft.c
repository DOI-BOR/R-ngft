#include "ngft.h"
#include "ngft_proto.h"

/* be sure to include these last so they can override some windows defs */
#include <R.h>
#include <Rinternals.h>

DllExport SEXP CALLngft_1dComplex64(SEXP ts_d, SEXP dt_d, SEXP gauss_l, SEXP image_dim_i,
																		 SEXP by_part_l, SEXP all_freqs_l, SEXP ind_map_l)
{
	int ii, jj, pcnt = 0, nf, nt;
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
	double epsilon = -1;

	double *ts = REAL(ts_d);
	int n_samples = length(ts_d);
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];
	BOOL gaussian_window = asLogical(gauss_l) == NA_LOGICAL ? TRUE : LOGICAL(gauss_l)[0];
	int image_dim = asInteger(image_dim_i) == NA_INTEGER ? -1 : INTEGER(image_dim_i)[0];
	BOOL by_part = asLogical( by_part_l ) == NA_LOGICAL ? TRUE : LOGICAL( by_part_l )[0];
	BOOL all_freqs = asLogical( all_freqs_l ) == NA_LOGICAL ? FALSE : LOGICAL( all_freqs_l )[0];
	BOOL ind_map = asLogical( ind_map_l ) == NA_LOGICAL ? FALSE : LOGICAL( ind_map_l )[0];

	if ( n_samples < 3 )
		error("time series must have at least 3 points");
	if ( dt <= 0 )
		error("dt must be positive");

	/* create partitions */
	partitions = ngft_FrequencyPartitions( n_samples, epsilon );
	window_fn = gaussian_window ? gaussian : box;

	// copy the data
	stride = 1;
	if ( (cts = calloc( n_samples, sizeof( *cts ) )) == NULL )
		oops( "CALLngft_1dComplex64","can't get space" );
	for ( ii = 0 ; ii < n_samples ; ii++ )
		cts[ii].r = ts[ii];

	// Call 1D GFT Function
	ngft_1dComplex64(&cts, &n_samples, &partitions, window_fn, stride);

	// get complex image
	image = ngft_1d_InterpolateNN(cts, partitions, image_dim, by_part, all_freqs, ind_map);
	nf = image->y_centers->count;
	nt = image->x_centers->count;

	/* allocate space for R structures for complex transform, and copy ngft output */

	// S-Transform
	cts_st = PROTECT(allocVector(CPLXSXP, n_samples)); pcnt++;
	ctp = COMPLEX(cts_st);
	for ( ii = 0 ; ii < n_samples ; ii++ ) {
		ctp[ii].r = cts[ii].r;
		ctp[ii].i = cts[ii].i;
	}

	// Image, with frequencies on the vertical axis and times on the horizontal axis
	cts_img = PROTECT(allocMatrix(CPLXSXP, nf, nt)); pcnt++;
	ctp = COMPLEX(cts_img);
	for ( ii = 0 ; ii < nf ; ii++ ) {
		for ( jj = 0 ; jj < nt ; jj++ ) {
			int kk = ii * nt + jj;	// NGFT stores images by row
			int ll = jj * nf + ii;	// R stores images by column
			ctp[ll].r = image->img->values[kk].r;
			ctp[ll].i = image->img->values[kk].i;
		}
	}

	// image width and height
	img_wd = PROTECT(allocVector(INTSXP, 1)); pcnt++; INTEGER(img_wd)[0] = nt; // time
	img_ht = PROTECT(allocVector(INTSXP, 1)); pcnt++; INTEGER(img_ht)[0] = nf; // frequency

	// frequency centers
	freq_c = PROTECT(allocVector(REALSXP, nf)); pcnt++;
	dp = REAL( freq_c );
	for ( ii = 0 ; ii < nf ; ii++ )
		dp[ii] = image->y_centers->values[ii];

	// time centers
	time_c = PROTECT(allocVector(REALSXP, nt)); pcnt++;
	dp = REAL( time_c );
	for ( ii = 0 ; ii < nt ; ii++ )
		dp[ii] = image->x_centers->values[ii];

	// other returned items
	win_name = PROTECT(allocVector(STRSXP, 1)); pcnt++;
	SET_STRING_ELT(win_name, 0, mkChar(gaussian_window ? "Gaussian" : "Box"));

	/* put the return values into a list */
	ret_l = PROTECT(allocVector(VECSXP, 10)); pcnt++;
	names_s = PROTECT(allocVector(VECSXP, 10)); pcnt++;
	SET_VECTOR_ELT(ret_l, 0, cts_st); SET_VECTOR_ELT(names_s, 0, mkChar("DST"));
	SET_VECTOR_ELT(ret_l, 1, cts_img); SET_VECTOR_ELT(names_s, 1, mkChar("image"));
	SET_VECTOR_ELT(ret_l, 2, img_wd); SET_VECTOR_ELT(names_s, 2, mkChar("wd"));
	SET_VECTOR_ELT(ret_l, 3, img_ht); SET_VECTOR_ELT(names_s, 3, mkChar("ht"));
	SET_VECTOR_ELT(ret_l, 4, freq_c); SET_VECTOR_ELT(names_s, 4, mkChar("f.centers"));
	SET_VECTOR_ELT(ret_l, 5, time_c); SET_VECTOR_ELT(names_s, 5, mkChar("t.centers"));
	SET_VECTOR_ELT(ret_l, 6, win_name); SET_VECTOR_ELT(names_s, 6, mkChar("win"));
	SET_VECTOR_ELT(ret_l, 7, by_part_l); SET_VECTOR_ELT(names_s, 7, mkChar("by.partition"));
	SET_VECTOR_ELT(ret_l, 8, all_freqs_l); SET_VECTOR_ELT(names_s, 8, mkChar("all.freqs"));
	SET_VECTOR_ELT(ret_l, 9, ind_map_l); SET_VECTOR_ELT(names_s, 9, mkChar("ind.map"));
	setAttrib(ret_l, R_NamesSymbol, names_s);

	UNPROTECT(pcnt);

	free(cts);
	freeDImage(image);
	ngft_FreeFreqPartitions(partitions);

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
	double epsilon = -1;

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
	partitions = ngft_FrequencyPartitions( n_samples, epsilon );
	stride = 1;

	if ( (cts = calloc( n_samples, sizeof( *cts ) )) == NULL )
		oops( "CALLngft_1dComplex64Inv","can't get space" );
	for ( ii = 0 ; ii < n_samples ; ii++ ) {
		cts[ii].r = dst[ii].r;
		cts[ii].i = dst[ii].i;
	}

	// Call 1D GFT Inverse Function
	ngft_1dComplex64Inv(&cts, &partitions, window_fn, stride);

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
